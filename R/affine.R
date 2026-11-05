# Sampling statements whose left-hand side is an affine function of one variable,
# e.g. `(r + 1) / 2 ~ beta(a, b)` or `1 - theta ~ exponential(l)`: the variable is
# drawn by inverting e = s * v + o, and its density gains log|s|.

# Var/Index chains (maximal) inside an expression
var_chains <- function(e) {
  if (!is.list(e) || is.null(e$node)) return(list())
  if (e$node %in% c("Var", "Index") && !is.null(outcome_base(e))) return(list(e))
  kids <- switch(e$node, BinOp = list(e$lhs, e$rhs), PrefixOp = , PostfixOp = list(e$arg), Call = e$args,
                 Ternary = list(e$cond, e$yes, e$no), list())
  unlist(lapply(kids, var_chains), recursive = FALSE)
}

# e = s * v + o for the chain of base `b`; s and o are IR expressions or NULL (1 and 0)
affine_decompose <- function(e, b) {
  contains <- function(x) b %in% expr_vars(x)
  if (e$node %in% c("Var", "Index") && identical(outcome_base(e), b)) return(list(target = e, s = NULL, o = NULL))
  num <- function(v) list(node = "Real", value = v, type = list(base = "real", adims = 0L), ad = "DataOnly")
  bin <- function(l, op, r) if (is.null(l)) r else if (is.null(r)) l else
    list(node = "BinOp", lhs = l, op = op, rhs = r, type = e$type,
         ad = if (identical(l$ad, "DataOnly") && identical(r$ad, "DataOnly")) "DataOnly" else "AutoDiffable")
  neg <- function(x) if (is.null(x)) num("-1") else list(node = "PrefixOp", op = "-", arg = x, type = x$type, ad = x$ad)
  if (e$node == "PrefixOp" && e$op == "-") {
    inner <- affine_decompose(e$arg, b)
    if (is.null(inner)) return(NULL)
    return(list(target = inner$target, s = neg(inner$s), o = if (is.null(inner$o)) NULL else neg(inner$o)))
  }
  if (e$node != "BinOp") return(NULL)
  in_l <- contains(e$lhs)
  in_r <- contains(e$rhs)
  if (in_l == in_r) return(NULL)
  inner <- affine_decompose(if (in_l) e$lhs else e$rhs, b)
  if (is.null(inner)) return(NULL)
  other <- if (in_l) e$rhs else e$lhs
  scalar_other <- type_rank(other$type) == 0
  switch(e$op,
    "+" = list(target = inner$target, s = inner$s, o = bin(inner$o, "+", other)),
    "-" = if (in_l) list(target = inner$target, s = inner$s, o = bin(inner$o, "-", other))
          else list(target = inner$target, s = neg(inner$s), o = bin(other, "-", inner$o)),
    "*" = , ".*" = if (e$op == "*" && !scalar_other && type_rank(e$type) > 0) NULL else
      list(target = inner$target, s = bin(other, ".*", inner$s), o = if (is.null(inner$o)) NULL else bin(other, ".*", inner$o)),
    "/" = , "./" = if (!in_l) NULL else
      list(target = inner$target, s = if (is.null(inner$s)) bin(num("1"), "./", other) else bin(inner$s, "./", other),
           o = if (is.null(inner$o)) NULL else bin(inner$o, "./", other)),
    NULL)
}

# the variable a sampling statement draws or scores, seen through to_vector() and affine maps
sampled_target <- function(e, ctx) {
  t <- unwrap_reshape(e)
  if (t$node %in% c("Var", "Index") && !is.null(outcome_base(t))) return(t)
  aff <- affine_outcome(e, ctx)
  if (is.null(aff)) NULL else aff$target
}

sampled_base <- function(e, ctx) {
  t <- sampled_target(e, ctx)
  if (is.null(t)) NULL else outcome_base(t)
}

affine_outcome <- function(e, ctx) {
  if (e$node %in% c("Var", "Index", "Call")) return(NULL)
  bases <- vapply(var_chains(e), outcome_base, "")
  roles <- vapply(bases, function(b) sym_role(ctx, b), "")
  cand <- unique(bases[roles %in% c("param", "data", "tdata")])
  cand <- cand[vapply(cand, function(b) sum(bases == b) == 1, logical(1))]
  decomp <- Filter(Negate(is.null), lapply(setNames(cand, cand), function(b) affine_decompose(e, b)))
  if (!length(decomp)) return(NULL)
  observed <- names(decomp)[vapply(names(decomp), function(b) sym_role(ctx, b) %in% c("data", "tdata"), logical(1))]
  pick <- if (length(observed) == 1) observed else if (length(decomp) == 1) names(decomp) else NULL
  if (is.null(pick)) return(NULL)
  decomp[[pick]]
}

# R code: v = (draw - o) / s
affine_inverse_code <- function(aff, draw, ctx) {
  out <- draw
  if (!is.null(aff$o)) out <- paste0("(", out, ") - ", wrap(translate_expr(aff$o, ctx), PREC_ADD + 1))
  if (!is.null(aff$s)) out <- paste0("(", out, ") / ", wrap(translate_expr(aff$s, ctx), PREC_MUL + 1))
  out
}

affine_log_jacobian <- function(aff, ctx) {
  if (is.null(aff$s)) return(NULL)
  paste0("log(abs(", translate_expr(aff$s, ctx)$code, "))")
}

# value of R code made only of numeric literals and arithmetic, else NA
constant_value <- function(code) {
  if (!grepl("^[-+*/(). 0-9eE]+$", code)) return(NA_real_)
  as.numeric(eval(parse(text = code), envir = baseenv()))
}

# the sign of s when it is a numeric constant (needed to map bounds of v to bounds of e)
affine_sign <- function(aff, ctx) {
  if (is.null(aff$s)) return(1)
  sign(constant_value(translate_expr(aff$s, ctx)$code))
}

# bounds of e from bounds of v (IR nodes), for a numeric s
affine_bounds <- function(aff, bounds, ctx) {
  sgn <- affine_sign(aff, ctx)
  if (is.na(sgn)) {
    stan2r_unsupported("affine left-hand side with a non-constant scale on a bounded parameter", ctx$line,
                       "expression on left of ~")
  }
  map <- function(b) {
    if (is.null(b)) return(NULL)
    out <- b
    if (!is.null(aff$s)) out <- list(node = "BinOp", lhs = aff$s, op = ".*", rhs = out, type = b$type, ad = b$ad)
    if (!is.null(aff$o)) out <- list(node = "BinOp", lhs = out, op = "+", rhs = aff$o, type = b$type, ad = b$ad)
    out
  }
  lo <- map(bounds$lower)
  hi <- map(bounds$upper)
  if (sgn < 0) list(lower = hi, upper = lo) else list(lower = lo, upper = hi)
}
