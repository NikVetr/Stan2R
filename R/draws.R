# Posterior draws -> R values with the declared Stan shape.
# Column names may be CmdStan style (`theta[2,1]`) or the data.frame() mangling
# of it (`theta.2.1.`).

draw_columns_pattern <- function(var_name) {
  paste0("^", var_name, "(\\[[0-9, ]+\\]|(\\.[0-9]+)+\\.?)?$")
}

#' Columns of a draws table that belong to one Stan variable.
#' @param var_name Stan variable name.
#' @param samps draws table (data frame or data.table) with CmdStan-style column
#'   names (`theta[1,2]`) or their `data.frame()` form (`theta.1.2.`).
#' @return a data frame with the matching columns.
#' @export
subset_samps <- function(var_name, samps) {
  cols <- grep(draw_columns_pattern(var_name), names(samps), value = TRUE)
  out <- lapply(cols, function(cl) samps[[cl]])
  names(out) <- cols
  as.data.frame(out, check.names = FALSE)
}

# parse "theta[2,1]" / "theta.2.1." into an integer matrix of indices (one row per column)
draw_indices <- function(var_name, cols) {
  idx <- sub(paste0("^", var_name), "", cols)
  idx <- gsub("^\\[|\\]$|^\\.|\\.$", "", idx)
  if (all(idx == "")) return(matrix(integer(0), nrow = length(cols), ncol = 0))
  parts <- strsplit(idx, "[,.]\\s*")
  lens <- lengths(parts)
  if (length(unique(lens)) != 1) stop("Stan2R: inconsistent index depth in draws for '", var_name, "'", call. = FALSE)
  matrix(as.integer(unlist(parts)), nrow = length(cols), byrow = TRUE)
}

# one draw (a named numeric vector of a row) -> array shaped by the index columns
shape_draw <- function(var_name, values, cols) {
  ix <- draw_indices(var_name, cols)
  if (ncol(ix) == 0) return(unname(values[[1]]))
  dims <- apply(ix, 2, max)
  if (length(dims) == 1) {
    out <- rep(NA_real_, dims)
    out[ix[, 1]] <- values
    return(out)
  }
  out <- array(NA_real_, dim = dims)
  out[ix] <- values
  if (length(dims) == 2) out <- matrix(out, dims[1], dims[2])
  out
}

#' Reshape the draws of one Stan variable into R objects (one per row).
#' @param var_name Stan variable name.
#' @param df draws table restricted to that variable (see [subset_samps()]).
#' @return the shaped value for a single row, otherwise a list of them.
#' @export
munge_samps <- function(var_name, df) {
  df <- as.data.frame(df)
  shape_row <- function(i) shape_draw(var_name, unlist(df[i, , drop = TRUE]), names(df))
  if (nrow(df) == 1) return(shape_row(1))
  lapply(seq_len(nrow(df)), shape_row)
}

# R literal for a parameter's value in draw `sample_index`
draw_value_code <- function(samps, sample_index, name, decl, ctx) {
  cols <- grep(draw_columns_pattern(name), names(samps), value = TRUE)
  if (!length(cols)) stop(sprintf("Stan2R: parameter '%s' (line %d) not found in samps", name, decl$line), call. = FALSE)
  if (sample_index < 1 || sample_index > nrow(samps)) stop("Stan2R: sample_index ", sample_index, " is outside samps", call. = FALSE)
  value <- shape_draw(name, vapply(cols, function(cl) as.numeric(samps[[cl]][sample_index]), 0), cols)
  rank <- length(decl$stype$array_dims) + length(decl$stype$dims)
  r_literal(value, rank)
}

r_literal <- function(value, rank) {
  num <- function(x) {
    s <- as.character(x)
    strwrap(paste(s, collapse = ", "), width = 100, exdent = 4)
  }
  join <- function(x) paste(num(x), collapse = "\n")
  if (rank == 0) return(as.character(value))
  if (rank == 1) return(paste0("c(", join(value), ")"))
  d <- dim(value)
  if (rank == 2) return(paste0("matrix(c(", join(as.vector(value)), "), nrow = ", d[1], ", ncol = ", d[2], ")"))
  paste0("array(c(", join(as.vector(value)), "), dim = c(", paste(d, collapse = ", "), "))")
}
