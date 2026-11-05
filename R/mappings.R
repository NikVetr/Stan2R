# Stan library functions -> R code. Each entry is function(args, e, ctx) where
# `args` are translated arguments (list(code, prec, type)) and `e` is the call's
# IR node (its $type is the Stan return type). Distribution functions
# (_lpdf, _rng, _cdf, ...) are handled by distributions.R instead.

fn_args <- function(args) paste(vapply(args, `[[`, "", "code"), collapse = ", ")

# rename: Stan f(a, b) -> R g(a, b)
fn_rename <- function(r_fun) function(args, e, ctx) rx(paste0(r_fun, "(", fn_args(args), ")"))

# template with positional placeholders {1}, {2}, ... and {(1)} for parenthesised
fn_tpl <- function(template, prec = PREC_ATOM, helpers = character(0)) {
  force(template)
  function(args, e, ctx) {
    vals <- args
    names(vals) <- as.character(seq_along(args))
    ctx_use_helper(ctx, helpers)
    rx(fill_template(template, vals), prec)
  }
}

fn_helper <- function(name) {
  force(name)
  function(args, e, ctx) {
    ctx_use_helper(ctx, name)
    rx(paste0(name, "(", fn_args(args), ")"))
  }
}

fn_const <- function(code) function(args, e, ctx) rx(code)

# rows()/cols()/size() depend on the argument's Stan type
fn_rows <- function(args, e, ctx) {
  a <- args[[1]]
  switch(la_kind(a$type), vector = rx(paste0("length(", a$code, ")")), row_vector = rx("1"),
         matrix = rx(paste0("nrow(", a$code, ")")), stop("Stan2R: rows() of ", a$type$base))
}
fn_cols <- function(args, e, ctx) {
  a <- args[[1]]
  switch(la_kind(a$type), vector = rx("1"), row_vector = rx(paste0("length(", a$code, ")")),
         matrix = rx(paste0("ncol(", a$code, ")")), stop("Stan2R: cols() of ", a$type$base))
}
fn_size <- function(args, e, ctx) {
  a <- args[[1]]
  if (a$type$adims > 0 && type_rank(a$type) > 1) rx(paste0("dim(", a$code, ")[1]")) else rx(paste0("length(", a$code, ")"))
}

fn_rep_matrix <- function(args, e, ctx) {
  x <- args[[1]]
  if (length(args) == 3) return(rx(paste0("matrix(", x$code, ", ", args[[2]]$code, ", ", args[[3]]$code, ")")))
  if (la_kind(x$type) == "vector") return(rx(paste0("matrix(", x$code, ", nrow = length(", x$code, "), ncol = ", args[[2]]$code, ")")))
  rx(paste0("matrix(", x$code, ", nrow = ", args[[2]]$code, ", ncol = length(", x$code, "), byrow = TRUE)"))
}

fn_rep_array <- function(args, e, ctx) {
  x <- args[[1]]
  dims <- vapply(args[-1], `[[`, "", "code")
  inner <- type_rank(x$type)
  if (inner == 0 && length(dims) == 1) return(rx(paste0("rep(", x$code, ", ", dims, ")")))
  if (inner == 0) return(rx(paste0("array(", x$code, ", dim = c(", paste(dims, collapse = ", "), "))")))
  if (inner == 1 && length(dims) == 1) return(rx(paste0("matrix(", x$code, ", nrow = ", dims, ", ncol = length(", x$code, "), byrow = TRUE)")))
  stan2r_unsupported("rep_array of containers with more than one array dimension", ctx$line, "rep_array of containers")
}

# Stan arrays are row-major in to_array_1d / to_vector of arrays; R arrays are column-major
fn_flatten <- function(args, e, ctx) {
  a <- args[[1]]
  if (a$type$adims >= 2 || (a$type$adims >= 1 && type_rank(a$type) >= 2)) {
    return(rx(paste0("as.vector(aperm(", a$code, "))")))
  }
  rx(paste0("as.vector(", a$code, ")"))
}

fn_to_matrix <- function(args, e, ctx) {
  a <- args[[1]]
  if (length(args) == 1) {
    if (a$type$adims == 2) return(rx(a$code))
    if (a$type$adims == 1) return(rx(paste0("matrix(", a$code, ", nrow = ", NROW_code(a), ")")))
    return(rx(paste0("as.matrix(", a$code, ")")))
  }
  byrow <- if (length(args) == 4) paste0(", byrow = ", args[[4]]$code, " == 0") else ""
  if (a$type$adims >= 1 && type_rank(a$type) >= 2) {
    return(rx(paste0("matrix(as.vector(aperm(", a$code, ")), ", args[[2]]$code, ", ", args[[3]]$code, byrow, ")")))
  }
  rx(paste0("matrix(", a$code, ", ", args[[2]]$code, ", ", args[[3]]$code, byrow, ")"))
}
NROW_code <- function(a) paste0("NROW(", a$code, ")")

fn_append_row <- function(args, e, ctx) {
  if (la_kind(e$type) == "matrix") rx(paste0("rbind(", fn_args(args), ")")) else rx(paste0("c(", fn_args(args), ")"))
}
fn_append_col <- function(args, e, ctx) {
  if (la_kind(e$type) == "matrix") rx(paste0("cbind(", fn_args(args), ")")) else rx(paste0("c(", fn_args(args), ")"))
}
fn_append_array <- function(args, e, ctx) {
  if (type_rank(e$type) == 1) rx(paste0("c(", fn_args(args), ")")) else if (type_rank(e$type) == 2) rx(paste0("rbind(", fn_args(args), ")"))
  else stan2r_unsupported("append_array of arrays with more than two dimensions", ctx$line, "append_array")
}

fn_reverse <- function(args, e, ctx) {
  a <- args[[1]]
  if (type_rank(a$type) == 1) return(rx(paste0("rev(", a$code, ")")))
  ctx_use_helper(ctx, character(0))
  rx(paste0(wrap(a, PREC_ATOM), "[rev(seq_len(dim(", a$code, ")[1]))", strrep(", ", type_rank(a$type) - 1), ", drop = FALSE]"))
}

fn_min_max <- function(r_fun) function(args, e, ctx) rx(paste0(r_fun, "(", fn_args(args), ")"))

fn_quad_form <- function(args, e, ctx) {
  A <- args[[1]]
  B <- args[[2]]
  if (la_kind(B$type) == "vector") return(rx(paste0("sum(", wrap(B, PREC_MUL), " * (", A$code, " %*% ", wrap(B, PREC_SPECIAL + 1), "))")))
  rx(paste0("t(", B$code, ") %*% ", wrap(A, PREC_SPECIAL + 1), " %*% ", wrap(B, PREC_SPECIAL + 1)), PREC_SPECIAL)
}

fn_eigen <- function(which) function(args, e, ctx) {
  a <- args[[1]]$code
  if (which == "values") rx(paste0("rev(eigen(", a, ", symmetric = TRUE)$values)"))
  else rx(paste0("eigen(", a, ", symmetric = TRUE)$vectors[, rev(seq_len(nrow(", a, "))), drop = FALSE]"))
}

fn_dims <- function(args, e, ctx) {
  a <- args[[1]]
  if (type_rank(a$type) <= 1) rx(paste0("length(", a$code, ")")) else rx(paste0("dim(", a$code, ")"))
}

fn_log0 <- function(base_code) function(args, e, ctx) {
  if (length(args) == 0) rx(paste0("log(", base_code, ")")) else rx(paste0("log", base_code, "(", args[[1]]$code, ")"))
}

fn_segment <- function(args, e, ctx) {
  x <- args[[1]]
  rx(paste0(wrap(x, PREC_ATOM), "[", wrap(args[[2]], PREC_ADD), " - 1 + seq_len(", args[[3]]$code, ")]"))
}

fn_head_tail <- function(which) function(args, e, ctx) {
  x <- args[[1]]
  n <- args[[2]]$code
  if (type_rank(x$type) == 1) return(rx(paste0(which, "(", x$code, ", ", n, ")")))
  rx(paste0(which, "(", x$code, ", ", n, ")"))
}

fn_block <- function(args, e, ctx) {
  m <- wrap(args[[1]], PREC_ATOM)
  rx(paste0(m, "[", wrap(args[[2]], PREC_ADD), " - 1 + seq_len(", args[[4]]$code, "), ",
            wrap(args[[3]], PREC_ADD), " - 1 + seq_len(", args[[5]]$code, "), drop = FALSE]"))
}

fn_sub_col <- function(args, e, ctx) {
  m <- wrap(args[[1]], PREC_ATOM)
  rx(paste0(m, "[", wrap(args[[2]], PREC_ADD), " - 1 + seq_len(", args[[4]]$code, "), ", args[[3]]$code, "]"))
}
fn_sub_row <- function(args, e, ctx) {
  m <- wrap(args[[1]], PREC_ATOM)
  rx(paste0(m, "[", args[[2]]$code, ", ", wrap(args[[3]], PREC_ADD), " - 1 + seq_len(", args[[4]]$code, ")]"))
}

fn_multinomial_rng <- function(prob_code) function(args, e, ctx) {
  if (grepl("softmax", prob_code)) ctx_use_helper(ctx, "softmax")
  rx(paste0("as.vector(rmultinom(1, size = ", args[[2]]$code, ", prob = ", sprintf(prob_code, args[[1]]$code), "))"))
}

# b / A for a row vector or matrix b
fn_mdivide_right <- function(args, e, ctx) {
  code <- paste0(wrap(args[[1]], PREC_SPECIAL), " %*% solve(", args[[2]]$code, ")")
  if (la_kind(e$type) == "row_vector") rx(paste0("as.vector(", code, ")")) else rx(code, PREC_SPECIAL)
}

fn_lkj_rng <- function(helper) function(args, e, ctx) {
  ctx_use_helper(ctx, helper)
  rx(paste0(helper, "(diag(", args[[1]]$code, "), ", args[[2]]$code, ")"))
}

stan_function_map <- list(
  # arithmetic and special functions
  abs = fn_rename("abs"), fabs = fn_rename("abs"), sqrt = fn_rename("sqrt"), cbrt = fn_tpl("sign({1}) * abs({1})^(1 / 3)", PREC_MUL),
  square = fn_tpl("{(1)}^2", PREC_POW), exp = fn_rename("exp"), exp2 = fn_tpl("2^{(1)}", PREC_POW),
  log = fn_rename("log"), log2 = fn_log0("2"), log10 = fn_log0("10"), log1p = fn_rename("log1p"),
  expm1 = fn_rename("expm1"), log1m = fn_tpl("log1p(-{(1)})"), log1p_exp = fn_tpl("log1p(exp({1}))"),
  log1m_exp = fn_tpl("log(-expm1({1}))"), log_diff_exp = fn_tpl("{(1)} + log(-expm1({(2)} - {(1)}))", PREC_ADD),
  log_inv_logit = fn_tpl("plogis({1}, log.p = TRUE)"), log1m_inv_logit = fn_tpl("plogis({1}, lower.tail = FALSE, log.p = TRUE)"),
  inv_logit = fn_rename("plogis"), logit = fn_rename("qlogis"), inv_cloglog = fn_tpl("1 - exp(-exp({1}))", PREC_ADD),
  Phi = fn_rename("pnorm"), inv_Phi = fn_rename("qnorm"), Phi_approx = fn_tpl("plogis(0.07056 * {(1)}^3 + 1.5976 * {(1)})"),
  erf = fn_tpl("2 * pnorm({(1)} * sqrt(2)) - 1", PREC_ADD), erfc = fn_tpl("2 * pnorm(-{(1)} * sqrt(2))", PREC_MUL),
  inv = fn_tpl("1 / {(1)}", PREC_MUL), inv_sqrt = fn_tpl("1 / sqrt({1})", PREC_MUL), inv_square = fn_tpl("1 / {(1)}^2", PREC_MUL),
  pow = fn_tpl("{(1)}^{(2)}", PREC_POW), hypot = fn_tpl("sqrt({(1)}^2 + {(2)}^2)"), fma = fn_tpl("{(1)} * {(2)} + {(3)}", PREC_ADD),
  fdim = fn_tpl("pmax({(1)} - {(2)}, 0)"), fmin = fn_rename("pmin"), fmax = fn_rename("pmax"),
  fmod = fn_tpl("{(1)} - {(2)} * trunc({(1)} / {(2)})", PREC_ADD),
  floor = fn_rename("floor"), ceil = fn_rename("ceiling"), trunc = fn_rename("trunc"), round = fn_helper("stan_round"),
  to_int = fn_rename("trunc"), int_step = fn_tpl("as.integer({(1)} > 0)"), step = fn_tpl("as.numeric({(1)} >= 0)"),
  is_inf = fn_rename("is.infinite"), is_nan = fn_rename("is.nan"),
  lgamma = fn_rename("lgamma"), tgamma = fn_rename("gamma"), digamma = fn_rename("digamma"), trigamma = fn_rename("trigamma"),
  lbeta = fn_rename("lbeta"), beta = fn_rename("beta"), lchoose = fn_rename("lchoose"), choose = fn_rename("choose"),
  binomial_coefficient_log = fn_rename("lchoose"), lmgamma = fn_tpl("0.25 * {(1)} * ({(1)} - 1) * log(pi) + sum(lgamma({(2)} + (1 - seq_len({1})) / 2))", PREC_ADD),
  inc_beta = fn_tpl("pbeta({3}, {1}, {2})"), gamma_p = fn_tpl("pgamma({2}, {1})"), gamma_q = fn_tpl("pgamma({2}, {1}, lower.tail = FALSE)"),
  log_falling_factorial = fn_tpl("lgamma({(1)} + 1) - lgamma({(1)} - {(2)} + 1)", PREC_ADD),
  log_rising_factorial = fn_tpl("lgamma({(1)} + {(2)}) - lgamma({1})", PREC_ADD),
  falling_factorial = fn_tpl("exp(lgamma({(1)} + 1) - lgamma({(1)} - {(2)} + 1))"),
  rising_factorial = fn_tpl("exp(lgamma({(1)} + {(2)}) - lgamma({1}))"),
  logical_negation = fn_tpl("!{(1)}", PREC_NOT),
  cos = fn_rename("cos"), sin = fn_rename("sin"), tan = fn_rename("tan"), acos = fn_rename("acos"),
  asin = fn_rename("asin"), atan = fn_rename("atan"), atan2 = fn_rename("atan2"), cosh = fn_rename("cosh"),
  sinh = fn_rename("sinh"), tanh = fn_rename("tanh"), acosh = fn_rename("acosh"), asinh = fn_rename("asinh"),
  atanh = fn_rename("atanh"), bessel_first_kind = fn_tpl("besselJ({2}, {1})"),
  bessel_second_kind = fn_tpl("besselY({2}, {1})"), modified_bessel_first_kind = fn_tpl("besselI({2}, {1})"),
  modified_bessel_second_kind = fn_tpl("besselK({2}, {1})"),
  log_modified_bessel_first_kind = fn_tpl("log(besselI({2}, {1}, expon.scaled = TRUE)) + {(2)}", PREC_ADD),
  log_inv_logit_diff = fn_tpl("log(plogis({1}) - plogis({2}))"),
  # constants
  pi = fn_const("pi"), e = fn_const("exp(1)"), sqrt2 = fn_const("sqrt(2)"), not_a_number = fn_const("NaN"),
  positive_infinity = fn_const("Inf"), negative_infinity = fn_const("-Inf"), machine_precision = fn_const(".Machine$double.eps"),
  # reductions
  sum = fn_rename("sum"), prod = fn_rename("prod"), mean = fn_rename("mean"), sd = fn_rename("sd"),
  variance = fn_rename("var"), min = fn_min_max("min"), max = fn_min_max("max"),
  log_sum_exp = fn_helper("log_sum_exp"), log_mix = fn_helper("log_mix"),
  softmax = fn_helper("softmax"), log_softmax = fn_helper("log_softmax"), cumulative_sum = fn_rename("cumsum"),
  # sizes
  size = fn_size, num_elements = fn_rename("length"), rows = fn_rows, cols = fn_cols, dims = fn_dims,
  # constructors
  rep_vector = fn_rename("rep"), rep_row_vector = fn_rename("rep"), rep_array = fn_rep_array, rep_matrix = fn_rep_matrix,
  zeros_vector = fn_tpl("rep(0, {1})"), zeros_row_vector = fn_tpl("rep(0, {1})"), zeros_array = fn_tpl("rep(0, {1})"),
  zeros_int_array = fn_tpl("rep(0L, {1})"), ones_vector = fn_tpl("rep(1, {1})"), ones_row_vector = fn_tpl("rep(1, {1})"),
  ones_array = fn_tpl("rep(1, {1})"), ones_int_array = fn_tpl("rep(1L, {1})"),
  uniform_simplex = fn_tpl("rep(1 / {(1)}, {1})"),
  linspaced_vector = fn_tpl("seq({2}, {3}, length.out = {1})"), linspaced_row_vector = fn_tpl("seq({2}, {3}, length.out = {1})"),
  linspaced_array = fn_tpl("seq({2}, {3}, length.out = {1})"), linspaced_int_array = fn_tpl("round(seq({2}, {3}, length.out = {1}))"),
  one_hot_vector = fn_tpl("replace(rep(0, {1}), {2}, 1)"), one_hot_row_vector = fn_tpl("replace(rep(0, {1}), {2}, 1)"),
  one_hot_array = fn_tpl("replace(rep(0, {1}), {2}, 1)"), one_hot_int_array = fn_tpl("replace(rep(0L, {1}), {2}, 1L)"),
  identity_matrix = fn_tpl("diag({1})"), diag_matrix = fn_tpl("diag({1}, nrow = length({1}))"),
  # conversions
  to_vector = fn_flatten, to_row_vector = fn_flatten, to_array_1d = fn_flatten, to_matrix = fn_to_matrix,
  to_array_2d = fn_tpl("{1}"), to_int = fn_rename("trunc"),
  # slicing and joining
  head = fn_head_tail("head"), tail = fn_head_tail("tail"), segment = fn_segment, block = fn_block,
  col = fn_tpl("{(1)}[, {2}]"), row = fn_tpl("{(1)}[{2}, ]"), sub_col = fn_sub_col, sub_row = fn_sub_row,
  append_row = fn_append_row, append_col = fn_append_col, append_array = fn_append_array, reverse = fn_reverse,
  sort_asc = fn_rename("sort"), sort_desc = fn_tpl("sort({1}, decreasing = TRUE)"),
  sort_indices_asc = fn_tpl("order({1})"), sort_indices_desc = fn_tpl("order({1}, decreasing = TRUE)"),
  rank = fn_tpl("sum({(1)} < {(1)}[{2}])"),
  # linear algebra
  dot_product = fn_tpl("sum({(1)} * {(2)})"), dot_self = fn_tpl("sum({(1)}^2)"),
  columns_dot_product = fn_tpl("colSums({(1)} * {(2)})"), rows_dot_product = fn_tpl("rowSums({(1)} * {(2)})"),
  columns_dot_self = fn_tpl("colSums({(1)}^2)"), rows_dot_self = fn_tpl("rowSums({(1)}^2)"),
  crossprod = fn_rename("crossprod"), tcrossprod = fn_rename("tcrossprod"),
  quad_form = fn_quad_form, quad_form_sym = fn_quad_form, quad_form_diag = fn_tpl("{(1)} * outer({2}, {2})", PREC_MUL),
  diag_pre_multiply = fn_tpl("{(1)} * {(2)}", PREC_MUL), diag_post_multiply = fn_tpl("t(t({1}) * {(2)})"),
  multiply_lower_tri_self_transpose = fn_tpl("tcrossprod({(1)} * lower.tri({1}, diag = TRUE))"),
  diagonal = fn_rename("diag"), add_diag = fn_tpl("{(1)} + diag({2}, nrow = nrow({1}))", PREC_ADD),
  trace = fn_tpl("sum(diag({1}))"), determinant = fn_rename("det"), log_determinant = fn_tpl("as.numeric(determinant({1})$modulus)"),
  log_determinant_spd = fn_tpl("as.numeric(determinant({1})$modulus)"),
  inverse = fn_rename("solve"), inverse_spd = fn_rename("solve"), chol2inv = fn_tpl("chol2inv(t({1}))"),
  cholesky_decompose = fn_tpl("t(chol({1}))"), mdivide_left = fn_rename("solve"), mdivide_left_spd = fn_rename("solve"),
  mdivide_left_tri_low = fn_tpl("as.vector(forwardsolve({1}, {2}))"), mdivide_right_tri_low = fn_mdivide_right,
  mdivide_right = fn_mdivide_right, mdivide_right_spd = fn_mdivide_right,
  eigenvalues_sym = fn_eigen("values"), eigenvectors_sym = fn_eigen("vectors"),
  singular_values = fn_tpl("svd({1})$d"), qr_thin_Q = fn_helper("qr_thin_Q"), qr_thin_R = fn_helper("qr_thin_R"),
  qr_Q = fn_helper("qr_Q"), qr_R = fn_helper("qr_R"), gp_exp_quad_cov = fn_tpl("{(2)}^2 * exp(-as.matrix(dist({1}))^2 / (2 * {(3)}^2))", PREC_MUL),
  cov_exp_quad = fn_tpl("{(2)}^2 * exp(-as.matrix(dist({1}))^2 / (2 * {(3)}^2))", PREC_MUL), matrix_exp = fn_tpl("expm::expm({1})"), matrix_power = fn_helper("matrix_power"),
  distance = fn_tpl("sqrt(sum(({(1)} - {(2)})^2))"), squared_distance = fn_tpl("sum(({(1)} - {(2)})^2)"),
  norm1 = fn_tpl("sum(abs({1}))"), norm2 = fn_tpl("sqrt(sum({(1)}^2))"),
  # random number generators whose signature differs from the distribution's
  multinomial_rng = fn_multinomial_rng("%s"), multinomial_logit_rng = fn_multinomial_rng("softmax(%s)"),
  lkj_corr_rng = fn_lkj_rng("rlkj_corr"), lkj_corr_cholesky_rng = fn_lkj_rng("rlkj_corr_cholesky")
)
