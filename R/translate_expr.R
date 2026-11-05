# Stan expressions (typed IR) -> R code. Every result is list(code, prec, type);
# `prec` is the R precedence of the outermost operator so parentheses are added
# only where R needs them (or where the Stan source had them).

PREC_ATOM <- 100
PREC_POW <- 90
PREC_UNARY <- 80
PREC_COLON <- 70
PREC_SPECIAL <- 60
PREC_MUL <- 50
PREC_ADD <- 40
PREC_CMP <- 30
PREC_NOT <- 20
PREC_AND <- 15
PREC_OR <- 10
PREC_IF <- 5

rx <- function(code, prec = PREC_ATOM, type = NULL) list(code = code, prec = prec, type = type)

wrap <- function(x, min_prec) if (x$prec < min_prec) paste0("(", x$code, ")") else x$code

type_rank <- function(t) t$adims + switch(t$base, vector = , row_vector = 1L, matrix = 2L, 0L)
is_scalar_type <- function(t) t$adims == 0 && t$base %in% c("int", "real")
is_int_type <- function(t) t$adims == 0 && t$base == "int"
la_kind <- function(t) if (t$adims > 0) "array" else if (t$base %in% c("int", "real")) "scalar" else t$base

r_reserved <- c("if", "else", "repeat", "while", "function", "for", "next", "break", "TRUE", "FALSE",
                "NULL", "Inf", "NaN", "NA", "NA_integer_", "NA_real_", "NA_character_", "in")

r_name <- function(name) if (name %in% r_reserved) paste0("`", name, "`") else name

translate_expr <- function(e, ctx) {
  out <- switch(e$node,
    Var = rx(r_name(e$name)),
    Int = rx(e$value),
    Real = rx(e$value),
    Str = rx(deparse(e$value)),
    BinOp = translate_binop(e, ctx),
    PrefixOp = translate_prefix(e, ctx),
    PostfixOp = {
      a <- translate_expr(e$arg, ctx)
      if (e$arg$type$base == "matrix" && e$arg$type$adims == 0) rx(paste0("t(", a$code, ")")) else a
    },
    Ternary = {
      cond <- translate_expr(e$cond, ctx)
      yes <- translate_expr(e$yes, ctx)
      no <- translate_expr(e$no, ctx)
      rx(paste0("if (", cond$code, ") ", wrap(yes, PREC_IF + 1), " else ", wrap(no, PREC_IF + 1)), PREC_IF)
    },
    Call = translate_call(e, ctx),
    Index = translate_index(e, ctx),
    ArrayExpr = translate_array_expr(e, ctx),
    RowVectorExpr = {
      elems <- lapply(e$elems, translate_expr, ctx = ctx)
      codes <- vapply(elems, `[[`, "", "code")
      if (e$type$base == "matrix") rx(paste0("rbind(", paste(codes, collapse = ", "), ")"))
      else rx(paste0("c(", paste(codes, collapse = ", "), ")"))
    },
    Imag = stan2r_unsupported("complex numbers", ctx$line, "complex numbers"),
    GetTarget = stan2r_unsupported("target()", ctx$line, "target()"),
    TupleProjection = , TupleExpr = stan2r_unsupported("tuples", ctx$line, "tuples"),
    stop("Stan2R: internal error, unknown expression node ", e$node, call. = FALSE))
  if (isTRUE(e$promote_complex) || grepl("complex", e$type$base)) {
    stan2r_unsupported("complex numbers", ctx$line, "complex numbers")
  }
  out$type <- e$type
  if (isTRUE(e$paren) && out$prec < PREC_ATOM) out <- rx(paste0("(", out$code, ")"), PREC_ATOM, e$type)
  out
}

translate_prefix <- function(e, ctx) {
  a <- translate_expr(e$arg, ctx)
  switch(e$op,
         "-" = rx(paste0("-", wrap(a, PREC_UNARY + 1)), PREC_UNARY),
         "+" = a,
         "!" = rx(paste0("!", wrap(a, PREC_ATOM)), PREC_NOT))
}

translate_binop <- function(e, ctx) {
  l <- translate_expr(e$lhs, ctx)
  r <- translate_expr(e$rhs, ctx)
  lt <- e$lhs$type
  rt <- e$rhs$type
  infix <- function(op, prec, right_assoc = FALSE) {
    left <- if (right_assoc) wrap(l, prec + 1) else wrap(l, prec)
    right <- if (right_assoc) wrap(r, prec) else wrap(r, prec + 1)
    rx(paste0(left, " ", op, " ", right), prec)
  }
  call2 <- function(fn, a = l, b = r) rx(paste0(fn, "(", a$code, ", ", b$code, ")"))
  matmul <- function() paste0(wrap(l, PREC_SPECIAL), " %*% ", wrap(r, PREC_SPECIAL + 1))
  switch(e$op,
    "+" = infix("+", PREC_ADD),
    "-" = infix("-", PREC_ADD),
    ".*" = infix("*", PREC_MUL),
    "./" = infix("/", PREC_MUL),
    "^" = , ".^" = infix("^", PREC_POW, right_assoc = TRUE),
    "*" = {
      kinds <- paste(la_kind(lt), la_kind(rt))
      if (la_kind(lt) == "scalar" || la_kind(rt) == "scalar") infix("*", PREC_MUL)
      else switch(kinds,
        "matrix vector" = , "row_vector matrix" = rx(paste0("as.vector(", matmul(), ")")),
        "matrix matrix" = rx(matmul(), PREC_SPECIAL),
        "row_vector vector" = rx(paste0("sum(", wrap(l, PREC_MUL), " * ", wrap(r, PREC_MUL + 1), ")")),
        "vector row_vector" = call2("outer"),
        stan2r_unsupported(sprintf("multiplication of %s by %s", kinds_text(lt), kinds_text(rt)),
                           ctx$line, "matrix multiplication shape"))
    },
    "/" = {
      if (is_int_type(lt) && is_int_type(rt)) rx(paste0("trunc(", wrap(l, PREC_MUL), " / ", wrap(r, PREC_MUL + 1), ")"))
      else if (la_kind(rt) == "scalar" || la_kind(lt) == "scalar") infix("/", PREC_MUL)
      else if (la_kind(lt) == "matrix" && la_kind(rt) == "matrix") rx(paste0(wrap(l, PREC_SPECIAL), " %*% solve(", r$code, ")"), PREC_SPECIAL)
      else if (la_kind(lt) == "row_vector" && la_kind(rt) == "matrix") rx(paste0("as.vector(", wrap(l, PREC_SPECIAL), " %*% solve(", r$code, "))"))
      else stan2r_unsupported(sprintf("division of %s by %s", kinds_text(lt), kinds_text(rt)), ctx$line, "matrix division shape")
    },
    "%/%" = rx(paste0("trunc(", wrap(l, PREC_MUL), " / ", wrap(r, PREC_MUL + 1), ")")),
    "%" = { ctx_use_helper(ctx, "stan_modulus"); call2("stan_modulus") },
    "\\" = call2("solve"),
    "==" = infix("==", PREC_CMP), "!=" = infix("!=", PREC_CMP),
    "<" = infix("<", PREC_CMP), "<=" = infix("<=", PREC_CMP),
    ">" = infix(">", PREC_CMP), ">=" = infix(">=", PREC_CMP),
    "&&" = infix("&&", PREC_AND), "||" = infix("||", PREC_OR),
    stop("Stan2R: internal error, unknown operator ", e$op, call. = FALSE))
}

kinds_text <- function(t) if (t$adims > 0) paste0(t$adims, "-d array of ", t$base) else t$base

ctx_use_helper <- function(ctx, names) {
  ctx$helpers <- unique(c(ctx$helpers, names))
  invisible(NULL)
}

#### indexing ####

# R code for the extent of dimension `pos` of x (used by open-ended ranges)
dim_extent <- function(x_code, pos, rank) {
  if (rank == 1) paste0("length(", x_code, ")") else paste0("dim(", x_code, ")[", pos, "]")
}

range_code <- function(lower, upper, ctx) {
  if (lower$code == "1") return(paste0("seq_len(", upper$code, ")"))
  if (grepl("^[0-9]+$", lower$code) && grepl("^[0-9]+$", upper$code) &&
      as.numeric(upper$code) >= as.numeric(lower$code)) {
    return(paste0(lower$code, ":", upper$code))
  }
  ctx_use_helper(ctx, "stan_range")
  paste0("stan_range(", lower$code, ", ", upper$code, ")")
}

# translate Stan indices into R subscripts for an object of R rank `rank`
index_codes <- function(indices, x_code, rank, ctx) {
  codes <- character(length(indices))
  singles <- integer(0)
  for (i in seq_along(indices)) {
    idx <- indices[[i]]
    codes[i] <- switch(idx$kind,
      single = { singles <- c(singles, i); translate_expr(idx$expr, ctx)$code },
      multi = translate_expr(idx$expr, ctx)$code,
      all = "",
      upfrom = range_code(translate_expr(idx$lower, ctx), rx(dim_extent(x_code, i, rank)), ctx),
      downfrom = paste0("seq_len(", translate_expr(idx$upper, ctx)$code, ")"),
      between = range_code(translate_expr(idx$lower, ctx), translate_expr(idx$upper, ctx), ctx))
  }
  if (rank > length(indices)) codes <- c(codes, rep("", rank - length(indices)))
  list(codes = codes, singles = singles)
}

translate_index <- function(e, ctx) {
  x <- translate_expr(e$arg, ctx)
  rank <- type_rank(e$arg$type)
  if (length(e$indices) > rank) stop("Stan2R: internal error, too many indices", call. = FALSE)
  ic <- index_codes(e$indices, x$code, rank, ctx)
  base <- wrap(x, PREC_ATOM)
  result_rank <- type_rank(e$type)
  sub <- paste(ic$codes, collapse = ", ")
  if (rank == 1 || result_rank <= 1) return(rx(paste0(base, "[", sub, "]")))
  if (length(ic$singles) == 0) return(rx(paste0(base, "[", sub, ", drop = FALSE]")))
  ctx_use_helper(ctx, "drop_dims")
  rx(paste0("drop_dims(", base, "[", sub, ", drop = FALSE], ", r_int_vector(ic$singles), ")"))
}

r_int_vector <- function(x) if (length(x) == 1) as.character(x) else paste0("c(", paste(x, collapse = ", "), ")")

translate_array_expr <- function(e, ctx) {
  elems <- lapply(e$elems, translate_expr, ctx = ctx)
  codes <- vapply(elems, `[[`, "", "code")
  inner_rank <- type_rank(e$elems[[1]]$type)
  if (inner_rank == 0) return(rx(paste0("c(", paste(codes, collapse = ", "), ")")))
  if (inner_rank == 1) return(rx(paste0("rbind(", paste(codes, collapse = ", "), ")")))
  ctx_use_helper(ctx, "stack_first_dim")
  rx(paste0("stack_first_dim(list(", paste(codes, collapse = ", "), "))"))
}

#### function calls ####

translate_call <- function(e, ctx) {
  if (e$origin == "UserDefined") {
    if (e$suffix == "FnTarget") stan2r_unsupported(sprintf("user-defined _lp function %s()", e$name), ctx$line, "_lp functions")
    args <- lapply(e$args, translate_expr, ctx = ctx)
    ctx$user_calls <- unique(c(ctx$user_calls, e$name))
    return(rx(paste0(r_name(e$name), "(", paste(vapply(args, `[[`, "", "code"), collapse = ", "), ")")))
  }
  parts <- split_dist_name(e$name)
  if (!is.null(parts)) {
    glm <- stan_glm_distributions[parts$dist]
    if (!is.na(glm) || !is.null(stan_distributions()[[parts$dist]])) {
      return(translate_distribution_call(e, parts, ctx))
    }
  }
  fn <- stan_function_map[[e$name]]
  if (is.null(fn)) {
    stan2r_unsupported(sprintf("Stan function %s() has no R translation", e$name), ctx$line,
                       paste0("function ", e$name))
  }
  args <- lapply(e$args, translate_expr, ctx = ctx)
  out <- fn(args, e, ctx)
  if (is.character(out)) out <- rx(out)
  out
}

# R code for the arguments of a distribution, named by the Stan parameter names;
# GLM forms are rewritten to their base distribution with a linear predictor.
distribution_args <- function(dist_name, arg_exprs, ctx) {
  glm <- stan_glm_distributions[dist_name]
  if (!is.na(glm)) {
    x <- arg_exprs[[1]]
    alpha <- translate_expr(arg_exprs[[2]], ctx)
    beta <- translate_expr(arg_exprs[[3]], ctx)
    xr <- translate_expr(x, ctx)
    lin <- if (la_kind(x$type) == "matrix") paste0("as.vector(", wrap(xr, PREC_SPECIAL), " %*% ", wrap(beta, PREC_SPECIAL + 1), ")")
           else paste0("sum(", wrap(xr, PREC_MUL), " * ", wrap(beta, PREC_MUL + 1), ")")
    eta <- rx(paste0(wrap(alpha, PREC_ADD), " + ", lin), PREC_ADD)
    base <- stan_distributions()[[unname(glm)]]
    rest <- lapply(arg_exprs[-(1:3)], translate_expr, ctx = ctx)
    vals <- c(list(eta), rest)
    names(vals) <- base$args
    return(list(def = base, values = vals, name = unname(glm)))
  }
  def <- stan_distributions()[[dist_name]]
  if (length(arg_exprs) != length(def$args)) {
    stan2r_unsupported(sprintf("%s with %d arguments", dist_name, length(arg_exprs)), ctx$line,
                       paste0("distribution ", dist_name))
  }
  vals <- lapply(arg_exprs, translate_expr, ctx = ctx)
  names(vals) <- def$args
  list(def = def, values = vals, name = dist_name)
}

dist_template <- function(da, which, ctx) {
  tpl <- da$def[[which]]
  if (is.null(tpl)) {
    what <- c(r = "random number generation", d = "density", p = "CDF", q = "quantile function")[which]
    stan2r_unsupported(sprintf("%s for distribution %s", what, da$name), ctx$line,
                       paste0(what, " ", da$name))
  }
  ctx_use_helper(ctx, template_helpers(tpl, da$def$helpers))
  tpl
}

# helpers a template actually calls
template_helpers <- function(tpl, helpers) {
  helpers[vapply(helpers, function(h) grepl(paste0(h, "("), tpl, fixed = TRUE), logical(1))]
}

# number of draws for an _rng call: 1 for scalar results, else the length of a container argument
rng_count <- function(e, args) {
  if (type_rank(e$type) == 0) return("1")
  for (i in seq_along(e$args)) if (type_rank(e$args[[i]]$type) > 0) return(paste0("length(", args[[i]]$code, ")"))
  "1"
}

translate_distribution_call <- function(e, parts, ctx) {
  is_lp <- parts$suffix %in% c("lpdf", "lupdf", "lpmf", "lupmf")
  if (parts$suffix == "rng") {
    da <- distribution_args(parts$dist, e$args, ctx)
    if (isTRUE(da$def$multivariate)) {
      like <- multivariate_rng_like(da, e, ctx)
      return(rx(fill_template(dist_template(da, "r", ctx), c(da$values, list(like = like)))))
    }
    n <- rng_count(e, da$values)
    return(rx(fill_template(dist_template(da, "r", ctx), c(da$values, list(n = n)))))
  }
  y <- translate_expr(e$args[[1]], ctx)
  da <- distribution_args(parts$dist, e$args[-1], ctx)
  if (is_lp) {
    d <- fill_template(dist_template(da, "d", ctx), c(da$values, list(y = y$code)))
    return(rx(paste0("sum(", d, ")")))
  }
  if (isTRUE(da$def$discrete) || !isTRUE(da$def$p_native)) {
    p <- fill_template(dist_template(da, "p", ctx), c(da$values, list(y = y$code, opts = "")))
    return(rx(switch(parts$suffix, cdf = paste0("prod(", p, ")"), lcdf = paste0("sum(log(", p, "))"),
                     lccdf = paste0("sum(log1p(-", p, "))"))))
  }
  opts <- switch(parts$suffix, cdf = "", lcdf = ", log.p = TRUE", lccdf = ", lower.tail = FALSE, log.p = TRUE")
  p <- fill_template(dist_template(da, "p", ctx), c(da$values, list(y = y$code, opts = opts)))
  rx(if (parts$suffix == "cdf") paste0("prod(", p, ")") else paste0("sum(", p, ")"))
}

# shape template for multivariate _rng results (an object of the result's shape)
multivariate_rng_like <- function(da, e, ctx) {
  vals <- da$values
  t <- e$type
  if (da$name %in% c("multi_normal", "multi_normal_cholesky", "multi_normal_prec",
                     "multi_student_t", "multi_student_t_cholesky")) {
    mu <- vals$mu
    return(if (t$adims > 0) paste0("matrix(NA_real_, nrow = NROW(", mu$code, "), ncol = NCOL(", mu$code, "))")
           else paste0("rep(NA_real_, length(", mu$code, "))"))
  }
  if (da$name == "dirichlet") return(paste0("rep(NA_real_, length(", vals$alpha$code, "))"))
  if (da$name %in% c("lkj_corr", "lkj_corr_cholesky")) {
    stan2r_unsupported(sprintf("%s_rng (use a declared matrix with a sampling statement)", da$name), ctx$line,
                       paste0(da$name, "_rng"))
  }
  if (da$name %in% c("multinomial", "multinomial_logit")) {
    stan2r_unsupported(sprintf("%s_rng", da$name), ctx$line, paste0(da$name, "_rng"))
  }
  "NULL"
}
