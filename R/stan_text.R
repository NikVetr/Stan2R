# Print IR expressions and declared types back as Stan source (for messages,
# comments, and Stan-to-Stan rewrites such as flatten_model()).

stan_binop_prec <- c("||" = 9, "&&" = 8, "==" = 7, "!=" = 7, "<" = 6, "<=" = 6, ">" = 6, ">=" = 6,
                     "+" = 5, "-" = 5, "*" = 4, "/" = 4, "%" = 4, "%/%" = 4, "\\" = 3, ".*" = 2, "./" = 2,
                     "^" = 0.5, ".^" = 0.5)

stan_expr_text <- function(e) stan_text_prec(e)$code

stan_text_prec <- function(e) {
  atom <- function(code) list(code = code, prec = 0)
  paren_if <- function(x, limit) if (x$prec > limit) paste0("(", x$code, ")") else x$code
  idx_text <- function(idx) {
    switch(idx$kind, single = , multi = stan_expr_text(idx$expr), all = "",
           upfrom = paste0(stan_expr_text(idx$lower), ":"), downfrom = paste0(":", stan_expr_text(idx$upper)),
           between = paste0(stan_expr_text(idx$lower), ":", stan_expr_text(idx$upper)))
  }
  out <- switch(e$node,
    Var = atom(e$name), Int = atom(e$value), Real = atom(e$value), Imag = atom(paste0(e$value, "i")),
    Str = atom(paste0("\"", e$value, "\"")),
    BinOp = {
      p <- stan_binop_prec[[e$op]]
      l <- stan_text_prec(e$lhs)
      r <- stan_text_prec(e$rhs)
      list(code = paste0(paren_if(l, p), " ", e$op, " ", paren_if(r, p - 0.01)), prec = p)
    },
    PrefixOp = list(code = paste0(e$op, paren_if(stan_text_prec(e$arg), 1)), prec = 1),
    PostfixOp = atom(paste0(paren_if(stan_text_prec(e$arg), 0), "'")),
    Ternary = list(code = paste0(stan_expr_text(e$cond), " ? ", stan_expr_text(e$yes), " : ", stan_expr_text(e$no)), prec = 10),
    Call = {
      args <- vapply(e$args, stan_expr_text, "")
      if (isTRUE(e$conditional) && length(args)) {
        atom(paste0(e$name, "(", args[1], if (length(args) > 1) paste0(" | ", paste(args[-1], collapse = ", ")), ")"))
      } else atom(paste0(e$name, "(", paste(args, collapse = ", "), ")"))
    },
    Index = atom(paste0(paren_if(stan_text_prec(e$arg), 0), "[", paste(vapply(e$indices, idx_text, ""), collapse = ", "), "]")),
    ArrayExpr = atom(paste0("{", paste(vapply(e$elems, stan_expr_text, ""), collapse = ", "), "}")),
    RowVectorExpr = atom(paste0("[", paste(vapply(e$elems, stan_expr_text, ""), collapse = ", "), "]")),
    GetTarget = atom("target()"),
    TupleProjection = atom(paste0(stan_expr_text(e$arg), ".", e$index)),
    TupleExpr = atom(paste0("(", paste(vapply(e$elems, stan_expr_text, ""), collapse = ", "), ")")),
    atom(paste0("<", e$node, ">")))
  if (isTRUE(e$paren) && out$prec > 0) out <- atom(paste0("(", out$code, ")"))
  out
}

# Stan declaration text for one variable, e.g. "vector<lower=0>[K] tau"
stan_decl_text <- function(decl, name, strip_constraints = FALSE) {
  st <- decl$stype
  tr <- decl$transform
  dims <- vapply(st$dims, stan_expr_text, "")
  constraint <- ""
  if (!strip_constraints) {
    if (tr$kind == "bounds") {
      parts <- c(if (!is.null(tr$lower)) paste0("lower=", stan_expr_text(tr$lower)),
                 if (!is.null(tr$upper)) paste0("upper=", stan_expr_text(tr$upper)))
      constraint <- paste0("<", paste(parts, collapse = ", "), ">")
    } else if (tr$kind == "affine") {
      parts <- c(if (!is.null(tr$offset)) paste0("offset=", stan_expr_text(tr$offset)),
                 if (!is.null(tr$multiplier)) paste0("multiplier=", stan_expr_text(tr$multiplier)))
      constraint <- paste0("<", paste(parts, collapse = ", "), ">")
    }
  }
  special <- tr$kind %in% c("simplex", "unit_vector", "sum_to_zero_vector", "ordered", "positive_ordered",
                            "cholesky_factor_corr", "corr_matrix", "cov_matrix", "cholesky_factor_cov",
                            "column_stochastic_matrix", "row_stochastic_matrix", "sum_to_zero_matrix")
  base <- if (special && !strip_constraints) tr$kind else st$base
  size <- if (special && !strip_constraints && tr$kind %in% c("cholesky_factor_corr", "corr_matrix", "cov_matrix")) dims[1] else dims
  type <- paste0(base, constraint, if (length(size)) paste0("[", paste(size, collapse = ", "), "]"))
  if (length(st$array_dims)) {
    type <- paste0("array[", paste(vapply(st$array_dims, stan_expr_text, ""), collapse = ", "), "] ", type)
  }
  paste(type, name)
}
