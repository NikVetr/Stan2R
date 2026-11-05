# Stan-to-Stan utilities built on the front end: block contents and flatten_model().

#' Declarations or sampling statements of one program block.
#'
#' @param stan_code Stan source.
#' @param block block name, e.g. "data", "parameters", "model".
#' @param type "declaration" (declared variables) or "sampling" (outcomes of
#'   sampling statements and of `target += *_lpdf(y | ...)` terms).
#' @return a data frame with `var_name` and, for declarations, `declaration_type`
#'   (int, real, vector, row_vector, matrix), `bounds`, `dimensions` and the full
#'   Stan `type`; for sampling statements, `distribution` and `line`.
#' @export
retrieve_block_objects <- function(stan_code, block = c("data", "transformed data", "parameters",
                                                        "transformed parameters", "model",
                                                        "generated quantities")[1],
                                   type = "declaration") {
  program <- if (inherits(stan_code, "stan2r_program")) stan_code else stan_program(stan_code)
  stmts <- program$blocks[[block]]
  if (type == "declaration") {
    rows <- list()
    for (s in stmts) {
      if (s$node != "Decl") next
      for (v in s$vars) {
        tr <- s$transform
        lower <- if (tr$kind == "bounds" && !is.null(tr$lower)) stan_expr_text(tr$lower) else "-Inf"
        upper <- if (tr$kind == "bounds" && !is.null(tr$upper)) stan_expr_text(tr$upper) else "Inf"
        dims <- vapply(c(s$stype$array_dims, s$stype$dims), stan_expr_text, "")
        rows[[length(rows) + 1]] <- data.frame(
          var_name = v$name, declaration_type = s$stype$base,
          bounds = paste0("lower=", lower, ", upper=", upper),
          dimensions = if (length(dims)) paste(dims, collapse = ", ") else "1",
          type = sub(paste0(" ", v$name, "$"), "", stan_decl_text(s, v$name)),
          stringsAsFactors = FALSE)
      }
    }
    if (!length(rows)) return(data.frame(var_name = character(0), declaration_type = character(0),
                                         bounds = character(0), dimensions = character(0), type = character(0)))
    return(do.call(rbind, rows))
  }
  if (type == "sampling") {
    rows <- list()
    note <- function(outcome, dist, line) {
      b <- outcome_base(outcome)
      rows[[length(rows) + 1]] <<- data.frame(var_name = b %||% stan_expr_text(outcome), distribution = dist,
                                              line = line, stringsAsFactors = FALSE)
    }
    walk_stmts(stmts, function(s) {
      if (s$node == "Tilde") note(s$arg, s$dist, s$line)
      if (s$node %in% c("TargetPE", "Assign")) {
        for (cl in ir_find(s, is_lp_call)) note(cl$args[[1]], lp_call_dist(cl), s$line)
      }
    })
    if (!length(rows)) return(data.frame(var_name = character(0), distribution = character(0), line = integer(0)))
    return(do.call(rbind, rows))
  }
  stop("Stan2R: type must be 'declaration' or 'sampling'", call. = FALSE)
}

#' Fix or move a multilevel model's scale parameters to flatten it.
#'
#' Parameters whose names match `scale_identifier` are removed from the
#' parameters block, their sampling statements are commented out, and they are
#' declared either in transformed data with the value `set_as`, or (with
#' `put_params_in_data = TRUE`) in the data block. Everything else, including
#' comments and layout, is left as written.
#'
#' @param stan_code Stan source.
#' @param scale_identifier regular expression matched against parameter names.
#' @param set_as value for the fixed parameters: one number, or a named vector
#'   with one value per matched parameter.
#' @param put_params_in_data move the parameters to the data block instead.
#' @return the new Stan program as a single string (checked with stanc).
#' @export
flatten_model <- function(stan_code, scale_identifier = "sd_", set_as = 1, put_params_in_data = FALSE) {
  program <- stan_program(stan_code)
  text <- program$source_text
  params <- program$symbols[program$symbols$role == "param", ]
  targets <- params$name[grepl(scale_identifier, params$name)]
  if (!length(targets)) {
    warning("flatten_model: no parameter name matches '", scale_identifier, "'", call. = FALSE)
    return(text)
  }
  if (length(set_as) > 1 && !all(targets %in% names(set_as))) {
    stop("flatten_model: set_as needs a value for each of ", paste(targets, collapse = ", "), call. = FALSE)
  }
  value_of <- function(name) if (length(set_as) > 1) set_as[[name]] else set_as
  edits <- list()
  edit <- function(start, end, replacement) edits[[length(edits) + 1]] <<- list(start = start, end = end, text = replacement)

  # 1. parameter declarations
  decl_text <- character(0)
  for (s in program$blocks$parameters) {
    if (s$node != "Decl") next
    names_here <- vapply(s$vars, `[[`, "", "name")
    hit <- names_here %in% targets
    if (!any(hit)) next
    keep <- vapply(s$vars[!hit], function(v) paste0(stan_decl_text(s, v$name), ";"), "")
    edit(s$span[["start"]], s$span[["end"]], paste(c(commented_out(s, text, standalone_body = FALSE), keep), collapse = "\n  "))
    for (v in s$vars[hit]) {
      if (put_params_in_data) {
        decl_text <- c(decl_text, paste0(stan_decl_text(s, v$name), ";"))
      } else {
        decl_text <- c(decl_text, paste0(stan_decl_text(s, v$name), " = ", fill_value_text(s, value_of(v$name)), ";"))
      }
    }
  }

  # 2. sampling statements on those parameters (Tilde, or target += with only such terms)
  visit <- function(stmts, parent_is_block) {
    for (s in stmts) {
      outcome <- NULL
      if (s$node == "Tilde") outcome <- outcome_base(s$arg)
      if (s$node == "TargetPE") {
        bases <- unique(vapply(ir_find(s$expr, is_lp_call), function(cl) outcome_base(cl$args[[1]]) %||% "", ""))
        if (length(bases) && all(bases %in% targets)) outcome <- bases[1]
        if (any(bases %in% targets) && !all(bases %in% targets)) {
          stop(sprintf("flatten_model: line %d mixes the scale parameters with other terms; edit it by hand", s$line), call. = FALSE)
        }
      }
      if (!is.null(outcome) && outcome %in% targets) {
        edit(s$span[["start"]], s$span[["end"]], commented_out(s, text, standalone_body = !parent_is_block))
        next
      }
      kids <- switch(s$node, Block = , Profile = list(list(stmts = s$stmts, block = TRUE)),
                     For = , ForEach = , While = list(list(stmts = list(s$body), block = FALSE)),
                     If = c(list(list(stmts = list(s$then), block = FALSE)),
                            if (!is.null(s$otherwise)) list(list(stmts = list(s$otherwise), block = FALSE))),
                     list())
      for (k in kids) visit(k$stmts, k$block)
    }
  }
  visit(program$blocks$model, TRUE)

  # 3. new declarations
  target_block <- if (put_params_in_data) "data" else "transformed data"
  info <- program$block_info[[target_block]]
  added <- paste0("  ", decl_text, "  // set by flatten_model()")
  if (!is.null(info)) {
    edit(info$close_brace, info$close_brace - 1L, paste0(paste(added, collapse = "\n"), "\n"))
  } else {
    after <- if (!put_params_in_data && !is.null(program$block_info$data)) program$block_info$data$span[["end"]]
             else if (!is.null(program$block_info$functions)) program$block_info$functions$span[["end"]] else 0L
    block_text <- paste0(if (after > 0) "\n" else "", target_block, " {\n", paste(added, collapse = "\n"), "\n}",
                         if (after == 0) "\n" else "")
    edit(after + 1L, after, block_text)
  }

  # apply edits from the end of the text backwards
  starts <- vapply(edits, `[[`, 0, "start")
  for (e in edits[order(starts, decreasing = TRUE)]) {
    text <- paste0(substr(text, 1, e$start - 1L), e$text, substr(text, e$end + 1L, nchar(text)))
  }
  stan_program(text)  # fails loudly if the rewrite is not valid Stan
  text
}

# Stan text that removes a statement but keeps it visible as a comment
commented_out <- function(s, text, standalone_body) {
  original <- substr(text, s$span[["start"]], s$span[["end"]])
  rest_of_line <- sub("\n.*$", "", substr(text, s$span[["end"]] + 1L, nchar(text)))
  prefix <- if (standalone_body) "; " else ""
  if (grepl("^\\s*(//.*)?$", rest_of_line)) {
    lines <- strsplit(original, "\n", fixed = TRUE)[[1]]
    return(paste0(prefix, paste(paste0("// ", lines), collapse = "\n")))
  }
  if (grepl("*/", original, fixed = TRUE)) {
    stop(sprintf("flatten_model: cannot comment out line %d automatically", s$line), call. = FALSE)
  }
  paste0(prefix, "/* ", original, " */")
}

# Stan expression filling a declared type with one value
fill_value_text <- function(decl, value) {
  st <- decl$stype
  dims <- vapply(st$dims, stan_expr_text, "")
  v <- format(value)
  inner <- switch(st$base, int = , real = v,
                  vector = paste0("rep_vector(", v, ", ", dims[1], ")"),
                  row_vector = paste0("rep_row_vector(", v, ", ", dims[1], ")"),
                  matrix = paste0("rep_matrix(", v, ", ", dims[1], ", ", dims[2], ")"),
                  stop("flatten_model: cannot fill a ", st$base, call. = FALSE))
  if (length(st$array_dims)) {
    inner <- paste0("rep_array(", inner, ", ", paste(vapply(st$array_dims, stan_expr_text, ""), collapse = ", "), ")")
  }
  inner
}
