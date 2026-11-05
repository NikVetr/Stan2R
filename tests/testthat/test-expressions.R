# Each case: Stan transformed-data statements -> expected R lines and values.

expr_case <- function(data_decls, td_stmts, dat) {
  code <- tr(td_program(data_decls, td_stmts), dat = dat)
  list(code = code, lines = code_lines(code), out = run_code(code, dat))
}

test_that("matrix products use %*% and keep Stan's shapes", {
  skip_if_no_stanc()
  X <- matrix(1:6, 3, 2)
  b <- c(0.5, -1)
  r <- expr_case(c("matrix[3, 2] X;", "vector[2] b;", "row_vector[3] w;"),
                 c("vector[3] Xb = X * b;", "row_vector[2] wX = w * X;", "real wb = w * (X * b);",
                   "matrix[3, 3] outer_b = (X * b) * w;", "matrix[2, 2] XtX = X' * X;",
                   "vector[3] elem = (X * b) .* (X * b);", "real scaled = 2 * w * (X * b);"),
                 list(X = X, b = b, w = c(1, 2, 3)))
  expect_true("Xb <- as.vector(X %*% b)" %in% r$lines)
  expect_true("XtX <- t(X) %*% X" %in% r$lines)
  expect_equal(r$out$Xb, as.vector(X %*% b))
  expect_equal(r$out$wX, as.vector(c(1, 2, 3) %*% X))
  expect_equal(r$out$wb, sum(c(1, 2, 3) * (X %*% b)))
  expect_equal(r$out$outer_b, outer(as.vector(X %*% b), c(1, 2, 3)))
  expect_equal(r$out$XtX, crossprod(X))
  expect_equal(r$out$elem, as.vector(X %*% b)^2)
  expect_equal(r$out$scaled, 2 * sum(c(1, 2, 3) * (X %*% b)))
})

test_that("integer division, modulus and precedence follow Stan", {
  skip_if_no_stanc()
  r <- expr_case(c("int a;", "int b;", "real x;"),
                 c("int q = a / b;", "int q2 = a %/% b;", "int m = a % b;", "real neg = -x^2;",
                   "real grouped = (x + 1) * (x - 1);", "real ratio = x / (a + b);",
                   "real t = a > b ? x : -x;", "int both = (a > 0) && (b > 0);"),
                 list(a = -7, b = 2, x = 3))
  expect_equal(r$out$q, -3)
  expect_equal(r$out$q2, -3)
  expect_equal(r$out$m, -1)
  expect_equal(r$out$neg, -9)
  expect_equal(r$out$grouped, 8)
  expect_equal(r$out$ratio, 3 / -5)
  expect_equal(r$out$t, -3)
  expect_true(r$out$both == 0)
})

test_that("indexing, slicing and ranges follow Stan", {
  skip_if_no_stanc()
  M <- matrix(1:12, 3, 4)
  A <- array(1:24, c(2, 3, 4))
  r <- expr_case(c("matrix[3, 4] M;", "array[2] matrix[3, 4] A;", "array[2] int idx;", "int n0;"),
                 c("row_vector[4] row2 = M[2];", "vector[3] col3 = M[:, 3];", "matrix[2, 4] rows = M[idx];",
                   "matrix[3, 2] tail_cols = M[, 3:];", "real a = A[2, 3, 4];", "matrix[3, 4] A1 = A[1];",
                   "matrix[1, 2] sub = A[2, 1:1, 2:3];", "vector[0] empty = M[1:n0, 1];",
                   "real total = 0;", "for (i in 1:n0) total += 1;"),
                 list(M = M, A = A, idx = c(3, 1), n0 = 0))
  expect_equal(r$out$row2, M[2, ])
  expect_equal(r$out$col3, M[, 3])
  expect_equal(r$out$rows, M[c(3, 1), , drop = FALSE])
  expect_equal(r$out$tail_cols, M[, 3:4])
  expect_equal(r$out$a, A[2, 3, 4])
  expect_equal(r$out$A1, A[1, , ])
  expect_equal(r$out$sub, matrix(A[2, 1, 2:3], 1, 2))
  expect_length(r$out$empty, 0)
  expect_equal(r$out$total, 0)
})

test_that("library functions with type-dependent meaning", {
  skip_if_no_stanc()
  v <- c(3, 1, 2)
  r <- expr_case(c("vector[3] v;", "array[2, 3] real a;"),
                 c("matrix[3, 2] rv = rep_matrix(v, 2);", "matrix[2, 3] rr = rep_matrix(v', 2);",
                   "array[6] real flat = to_array_1d(a);", "int n_rows = rows(v);", "int n_cols = cols(v');",
                   "vector[3] sm = softmax(v);", "real lse = log_sum_exp(v);", "matrix[3, 3] D = diag_matrix(v);",
                   "matrix[3, 3] L = cholesky_decompose(add_diag(rep_matrix(0.5, 3, 3), 1));"),
                 list(v = v, a = matrix(1:6, 2, 3, byrow = TRUE)))
  expect_equal(r$out$rv, matrix(v, 3, 2))
  expect_equal(r$out$rr, matrix(v, 2, 3, byrow = TRUE))
  expect_equal(r$out$flat, 1:6)
  expect_equal(c(r$out$n_rows, r$out$n_cols), c(3, 3))
  expect_equal(r$out$sm, exp(v) / sum(exp(v)))
  expect_equal(r$out$lse, log(sum(exp(v))))
  expect_equal(r$out$D, diag(v))
  S <- matrix(0.5, 3, 3) + diag(3)
  expect_equal(r$out$L %*% t(r$out$L), S)
  expect_true(all(r$out$L[upper.tri(r$out$L)] == 0))
})

test_that("unknown library functions are explicit errors naming the line", {
  skip_if_no_stanc()
  code <- td_program("vector[2] v;", c("real a = 1;", "real b = owens_t(a, 2);"))
  err <- tryCatch(tr(code), error = function(e) e)
  expect_s3_class(err, "stan2r_unsupported")
  expect_equal(err$line, 6L)
  expect_match(conditionMessage(err), "owens_t")
})
