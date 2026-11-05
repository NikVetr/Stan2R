# Front end: Stan source text -> checked, typed IR with a symbol table.
#   stanc (typed AST dump)  -> read_sexp()  -> ast_program()  -> attach_source_map()

stan2r_cache <- new.env(parent = emptyenv())

as_stan_lines <- function(stan_code) {
  if (!is.character(stan_code)) stop("Stan2R: stan_code must be a character vector", call. = FALSE)
  if (length(stan_code) == 1 && file.exists(stan_code) && grepl("\\.stan$", stan_code)) {
    stan_code <- readLines(stan_code, warn = FALSE)
  }
  unlist(strsplit(stan_code, "\n", fixed = TRUE), use.names = FALSE)
}

#' Parse a Stan program with stanc and return Stan2R's typed IR.
#'
#' @param stan_code Stan source (lines, a single string, or a path to a .stan file).
#' @param include_paths directories searched for `#include` files.
#' @param model_name name used in stanc messages.
#' @return a `stan2r_program`: blocks of IR statements with line numbers and
#'   comments, a symbol table of top-level declarations, and user functions.
#' @export
stan_program <- function(stan_code, include_paths = character(0), model_name = "model") {
  lines <- as_stan_lines(stan_code)
  key <- paste(c(lines, "\u001f", include_paths), collapse = "\n")
  cached <- stan2r_cache$programs[[key]]
  if (!is.null(cached)) return(cached)

  if (any(grepl("^\\s*#include", lines))) {
    lines <- stanc_expand_includes(lines, include_paths, model_name)
  }
  dump <- stanc_typed_ast(lines, model_name)
  blocks <- ast_program(read_sexp(dump$ast))
  mapped <- attach_source_map(blocks, lines)

  program <- structure(list(
    blocks = mapped$blocks,
    block_info = mapped$block_info,
    trailing_comments = mapped$trailing_comments,
    symbols = build_symbols(mapped$blocks),
    functions = collect_functions(mapped$blocks$functions),
    source = lines,
    source_text = mapped$text,
    stanc_warnings = dump$warnings,
    stanc_version = stanc_version_string()
  ), class = "stan2r_program")

  if (is.null(stan2r_cache$programs)) stan2r_cache$programs <- list()
  if (length(stan2r_cache$programs) >= 20) stan2r_cache$programs <- stan2r_cache$programs[-1]
  stan2r_cache$programs[[key]] <- program
  program
}

#' @export
print.stan2r_program <- function(x, ...) {
  present <- names(x$blocks)[!vapply(x$blocks, is.null, logical(1))]
  cat("<stan2r_program> blocks:", paste(present, collapse = ", "), "\n")
  cat("declared:", nrow(x$symbols), "top-level variables;", length(x$functions), "user functions\n")
  invisible(x)
}

# variable role used by the translator
block_roles <- c("data" = "data", "transformed data" = "tdata", "parameters" = "param",
                 "transformed parameters" = "tparam", "model" = "local",
                 "generated quantities" = "gq")

# Stan type of a declaration as text, e.g. "array[] vector", "real"
decl_type_text <- function(stype, transform) {
  base <- if (!is.null(transform) && !(transform$kind %in% c("none", "bounds", "affine", "tuple"))) transform$kind else stype$base
  if (length(stype$array_dims)) paste0("array[", strrep(",", length(stype$array_dims) - 1), "] ", base) else base
}

build_symbols <- function(blocks) {
  rows <- list()
  for (b in names(block_roles)) {
    for (s in blocks[[b]]) {
      if (!identical(s$node, "Decl")) next
      for (v in s$vars) {
        rows[[length(rows) + 1]] <- list(
          name = v$name, block = b, role = unname(block_roles[b]),
          base = s$stype$base, adims = length(s$stype$array_dims),
          rank = length(s$stype$array_dims) + length(s$stype$dims),
          type_text = decl_type_text(s$stype, s$transform),
          transform = s$transform$kind, line = s$line %||% NA_integer_,
          decl = s, has_init = !is.null(v$init))
      }
    }
  }
  if (length(rows) == 0) {
    return(data.frame(name = character(0), block = character(0), role = character(0),
                      base = character(0), adims = integer(0), rank = integer(0),
                      type_text = character(0), transform = character(0), line = integer(0),
                      stringsAsFactors = FALSE))
  }
  df <- data.frame(
    name = vapply(rows, `[[`, "", "name"), block = vapply(rows, `[[`, "", "block"),
    role = vapply(rows, `[[`, "", "role"), base = vapply(rows, `[[`, "", "base"),
    adims = vapply(rows, `[[`, 0L, "adims"), rank = vapply(rows, `[[`, 0L, "rank"),
    type_text = vapply(rows, `[[`, "", "type_text"), transform = vapply(rows, `[[`, "", "transform"),
    line = vapply(rows, function(r) as.integer(r$line), 0L), stringsAsFactors = FALSE)
  df$decl <- lapply(rows, `[[`, "decl")
  df
}

collect_functions <- function(stmts) {
  defs <- Filter(function(s) identical(s$node, "FunDef") && !isTRUE(s$forward), stmts)
  if (!length(defs)) return(list())
  names(defs) <- vapply(defs, `[[`, "", "name")
  dup <- unique(names(defs)[duplicated(names(defs))])
  if (length(dup)) {
    s <- defs[[which(names(defs) == dup[1])[2]]]
    stan2r_unsupported(sprintf("overloaded user-defined function '%s'", dup[1]), s$line,
                       "overloaded functions")
  }
  defs
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# Structured error for constructs Stan2R cannot translate faithfully.
stan2r_unsupported <- function(what, line = NA, construct = what) {
  where <- if (!is.null(line) && !is.na(line)) sprintf(" (line %d)", line) else ""
  cond <- structure(class = c("stan2r_unsupported", "error", "condition"),
                    list(message = paste0("Stan2R: unsupported: ", what, where),
                         call = NULL, construct = construct, line = line))
  stop(cond)
}
