test_that("s-expressions nest atoms and lists", {
  x <- Stan2R:::read_sexp("(a (b c) () (d (e)))")
  expect_equal(x, list("a", list("b", "c"), list(), list("d", list("e"))))
})

test_that("quoted atoms keep spaces, parentheses and escapes", {
  x <- Stan2R:::read_sexp('(PString "\\"a (b) \\\\n\\"")')
  expect_equal(x[[2]], '"a (b) \\n"')
})

test_that("unbalanced input is an error", {
  expect_error(Stan2R:::read_sexp("(a (b)"), "unbalanced")
  expect_error(Stan2R:::read_sexp("(a))"), "unbalanced")
})
