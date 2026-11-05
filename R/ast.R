# Convert the s-expression read from stanc's typed AST into Stan2R's
# intermediate representation (IR): plain R lists with a `node` field.
# This is the only file that knows the stanc dump format; an unrecognised
# node is an error that names the node, so format drift fails loudly.

ast_fail <- function(what, x) {
  shown <- if (is.character(x)) x else if (is.list(x) && length(x) && is.character(x[[1]])) x[[1]] else "<list>"
  stop(sprintf("Stan2R: unrecognised stanc AST %s '%s' (stanc %s). The stanc dump format may have changed.",
               what, shown, stanc_version_string()), call. = FALSE)
}

# a record prints as ((field value) (field value) ...); return a named list
sexp_record <- function(x) {
  if (!is.list(x)) ast_fail("record", x)
  out <- list()
  for (el in x) {
    if (!is.list(el) || length(el) < 1 || !is.character(el[[1]])) ast_fail("record field", el)
    out[el[[1]]] <- list(if (length(el) == 2) el[[2]] else el[-1])
  }
  out
}

# inline-record constructor: (Name (field value) ...) -> named list of fields
sexp_inline <- function(x) sexp_record(x[-1])

sexp_head <- function(x) if (is.character(x)) x else if (is.list(x) && length(x) && is.character(x[[1]])) x[[1]] else ""

sexp_option <- function(x) if (length(x) == 0) NULL else x[[1]]

#### types ####

# unsized type: list(base = "int"|"real"|"complex"|"vector"|"row_vector"|"matrix"|
#                       "complex_vector"|..., adims = number of array dimensions)
ast_utype <- function(x) {
  adims <- 0L
  while (is.list(x) && identical(x[[1]], "UArray")) {
    adims <- adims + 1L
    x <- x[[2]]
  }
  if (is.list(x)) {
    base <- switch(x[[1]], UTuple = "tuple", UFun = "function", ast_fail("type", x))
  } else {
    base <- switch(x,
                   UInt = "int", UReal = "real", UComplex = "complex",
                   UVector = "vector", URowVector = "row_vector", UMatrix = "matrix",
                   UComplexVector = "complex_vector", UComplexRowVector = "complex_row_vector",
                   UComplexMatrix = "complex_matrix", UMathLibraryFunction = "function",
                   ast_fail("type", x))
  }
  list(base = base, adims = adims)
}

# sized type: base, container dims (size expressions), array dims (outermost first)
ast_stype <- function(x) {
  array_dims <- list()
  while (is.list(x) && identical(x[[1]], "SArray")) {
    array_dims <- c(array_dims, list(ast_expr(x[[3]])))
    x <- x[[2]]
  }
  if (is.character(x)) {
    base <- switch(x, SInt = "int", SReal = "real", SComplex = "complex", ast_fail("declared type", x))
    dims <- list()
  } else {
    head <- x[[1]]
    base <- switch(head,
                   SVector = "vector", SRowVector = "row_vector", SMatrix = "matrix",
                   SComplexVector = "complex_vector", SComplexRowVector = "complex_row_vector",
                   SComplexMatrix = "complex_matrix", STuple = "tuple",
                   ast_fail("declared type", x))
    if (base == "tuple") {
      dims <- list()
    } else {
      # (SVector AoS e) / (SMatrix AoS e1 e2); the memory-pattern atom is skipped
      size_args <- x[-1]
      size_args <- size_args[!vapply(size_args, is.character, logical(1))]
      dims <- lapply(size_args, ast_expr)
    }
  }
  list(base = base, dims = dims, array_dims = array_dims)
}

ast_transformation <- function(x) {
  head <- sexp_head(x)
  with_exprs <- c("Lower", "Upper", "LowerUpper", "Offset", "Multiplier", "OffsetMultiplier")
  args <- if (head %in% with_exprs) lapply(x[-1], ast_expr) else list()
  switch(head,
         Identity = list(kind = "none"),
         Lower = list(kind = "bounds", lower = args[[1]], upper = NULL),
         Upper = list(kind = "bounds", lower = NULL, upper = args[[1]]),
         LowerUpper = list(kind = "bounds", lower = args[[1]], upper = args[[2]]),
         Offset = list(kind = "affine", offset = args[[1]], multiplier = NULL),
         Multiplier = list(kind = "affine", offset = NULL, multiplier = args[[1]]),
         OffsetMultiplier = list(kind = "affine", offset = args[[1]], multiplier = args[[2]]),
         Ordered = list(kind = "ordered"),
         PositiveOrdered = list(kind = "positive_ordered"),
         Simplex = list(kind = "simplex"),
         UnitVector = list(kind = "unit_vector"),
         SumToZero = list(kind = "sum_to_zero_vector"),
         SumToZeroMatrix = list(kind = "sum_to_zero_matrix"),
         CholeskyCorr = list(kind = "cholesky_factor_corr"),
         CholeskyCov = list(kind = "cholesky_factor_cov"),
         Correlation = list(kind = "corr_matrix"),
         Covariance = list(kind = "cov_matrix"),
         StochasticColumn = list(kind = "column_stochastic_matrix"),
         StochasticRow = list(kind = "row_stochastic_matrix"),
         TupleTransformation = list(kind = "tuple"),
         ast_fail("transformation", x))
}

#### expressions ####

ast_operator <- function(x) {
  ops <- c(Plus = "+", PPlus = "+", Minus = "-", PMinus = "-", Times = "*", Divide = "/",
           IntDivide = "%/%", Modulo = "%", LDivide = "\\", EltTimes = ".*", EltDivide = "./",
           Pow = "^", EltPow = ".^", Or = "||", And = "&&", Equals = "==", NEquals = "!=",
           Less = "<", Leq = "<=", Greater = ">", Geq = ">=", PNot = "!", Transpose = "'")
  if (!is.character(x) || !(x %in% names(ops))) ast_fail("operator", x)
  unname(ops[x])
}

ast_fun_kind <- function(x) {
  # (StanLib FnPlain) / (UserDefined (FnLpdf true))
  origin <- x[[1]]
  if (!(origin %in% c("StanLib", "UserDefined"))) ast_fail("function kind", x)
  suffix <- x[[2]]
  suffix_name <- sexp_head(suffix)
  if (!(suffix_name %in% c("FnPlain", "FnRng", "FnLpdf", "FnLpmf", "FnTarget", "FnJacobian")))
    ast_fail("function suffix", suffix)
  list(origin = origin, suffix = suffix_name)
}

ast_index <- function(x) {
  head <- sexp_head(x)
  switch(head,
         All = list(kind = "all"),
         Single = {
           e <- ast_expr(x[[2]])
           list(kind = if (e$type$adims > 0) "multi" else "single", expr = e)
         },
         Upfrom = list(kind = "upfrom", lower = ast_expr(x[[2]])),
         Downfrom = list(kind = "downfrom", upper = ast_expr(x[[2]])),
         Between = list(kind = "between", lower = ast_expr(x[[2]]), upper = ast_expr(x[[3]])),
         ast_fail("index", x))
}

ast_identifier_name <- function(x) sexp_record(x)$name

ast_expr <- function(x) {
  rec <- sexp_record(x)
  meta <- sexp_record(rec$emeta)
  type <- ast_utype(meta$type_)
  ad <- meta$ad_level
  e <- rec$expr
  head <- sexp_head(e)
  node <- switch(head,
    Variable = list(node = "Var", name = ast_identifier_name(e[[2]])),
    IntNumeral = list(node = "Int", value = e[[2]]),
    RealNumeral = list(node = "Real", value = e[[2]]),
    ImagNumeral = list(node = "Imag", value = e[[2]]),
    BinOp = list(node = "BinOp", lhs = ast_expr(e[[2]]), op = ast_operator(e[[3]]), rhs = ast_expr(e[[4]])),
    PrefixOp = list(node = "PrefixOp", op = ast_operator(e[[2]]), arg = ast_expr(e[[3]])),
    PostfixOp = list(node = "PostfixOp", arg = ast_expr(e[[2]]), op = ast_operator(e[[3]])),
    TernaryIf = list(node = "Ternary", cond = ast_expr(e[[2]]), yes = ast_expr(e[[3]]), no = ast_expr(e[[4]])),
    FunApp = ,
    CondDistApp = {
      kind <- ast_fun_kind(e[[2]])
      list(node = "Call", name = ast_identifier_name(e[[3]]), origin = kind$origin,
           suffix = kind$suffix, conditional = head == "CondDistApp",
           args = lapply(e[[4]], ast_expr))
    },
    Promotion = {
      inner <- ast_expr(e[[2]])
      to <- ast_utype(e[[3]][[1]])
      if (grepl("complex", to$base)) inner$promote_complex <- TRUE
      return(within_type(inner, type, ad))
    },
    Paren = return(within_type(ast_expr(e[[2]]), type, ad, paren = TRUE)),
    Indexed = list(node = "Index", arg = ast_expr(e[[2]]), indices = lapply(e[[3]], ast_index)),
    ArrayExpr = list(node = "ArrayExpr", elems = lapply(e[[2]], ast_expr)),
    RowVectorExpr = list(node = "RowVectorExpr", elems = lapply(e[[2]], ast_expr)),
    GetTarget = list(node = "GetTarget"),
    TupleProjection = list(node = "TupleProjection", arg = ast_expr(e[[2]]), index = e[[3]]),
    TupleExpr = list(node = "TupleExpr", elems = lapply(e[[2]], ast_expr)),
    ast_fail("expression", e))
  node$type <- type
  node$ad <- ad
  node
}

# a promoted or parenthesised expression keeps its inner node but takes the outer type
within_type <- function(inner, type, ad, paren = FALSE) {
  if (!identical(inner$type, type)) inner$promoted_from <- inner$type
  inner$type <- type
  inner$ad <- ad
  if (paren) inner$paren <- TRUE
  inner
}

#### statements ####

ast_lvalue <- function(x) {
  # (LValue ((lval L) (lmeta M))) where L is LVariable / LIndexed / LTupleProjection
  head <- sexp_head(x)
  if (head == "LTuplePack") ast_fail("lvalue (tuple unpacking)", x)
  if (head != "LValue") ast_fail("lvalue", x)
  ast_lval(x[[2]])
}

ast_lval <- function(x) {
  rec <- sexp_record(x)
  meta <- sexp_record(rec$lmeta)
  l <- rec$lval
  head <- sexp_head(l)
  node <- switch(head,
    LVariable = list(node = "Var", name = ast_identifier_name(l[[2]])),
    LIndexed = list(node = "Index", arg = ast_lval(l[[2]]), indices = lapply(l[[3]], ast_index)),
    LTupleProjection = list(node = "TupleProjection", arg = ast_lval(l[[2]]), index = l[[3]]),
    ast_fail("lvalue", l))
  node$type <- ast_utype(meta$type_)
  node$ad <- meta$ad_level
  node
}

ast_printables <- function(x) {
  lapply(x, function(p) {
    switch(sexp_head(p),
           PString = list(node = "Str", value = sexp_unquote_stan(p[[2]]), type = list(base = "string", adims = 0L), ad = "DataOnly"),
           PExpr = ast_expr(p[[2]]),
           ast_fail("printable", p))
  })
}

# PString atoms keep Stan's own quotes, e.g. "\"n=\"" -> n=
sexp_unquote_stan <- function(s) sub('^"(.*)"$', "\\1", s)

ast_truncation <- function(x) {
  head <- sexp_head(x)
  switch(head,
         NoTruncate = NULL,
         TruncateUpFrom = list(lower = ast_expr(x[[2]]), upper = NULL),
         TruncateDownFrom = list(lower = NULL, upper = ast_expr(x[[2]])),
         TruncateBetween = list(lower = ast_expr(x[[2]]), upper = ast_expr(x[[3]])),
         ast_fail("truncation", x))
}

ast_stmt <- function(x) {
  rec <- sexp_record(x)
  s <- rec$stmt
  head <- sexp_head(s)
  switch(head,
    Skip = list(node = "Skip"),
    Break = list(node = "Break"),
    Continue = list(node = "Continue"),
    ReturnVoid = list(node = "Return", value = NULL),
    Return = list(node = "Return", value = ast_expr(s[[2]])),
    Print = list(node = "Print", items = ast_printables(s[[2]])),
    Reject = list(node = "Reject", items = ast_printables(s[[2]])),
    FatalError = list(node = "FatalError", items = ast_printables(s[[2]])),
    TargetPE = list(node = "TargetPE", expr = ast_expr(s[[2]])),
    JacobianPE = list(node = "JacobianPE", expr = ast_expr(s[[2]])),
    NRFunApp = {
      kind <- ast_fun_kind(s[[2]])
      list(node = "NRFunApp", name = ast_identifier_name(s[[3]]), origin = kind$origin,
           suffix = kind$suffix, args = lapply(s[[4]], ast_expr))
    },
    Assignment = {
      f <- sexp_inline(s)
      op <- f$assign_op
      op <- if (identical(op, "Assign")) "=" else paste0(ast_operator(op[[2]]), "=")
      list(node = "Assign", lhs = ast_lvalue(f$assign_lhs), op = op, rhs = ast_expr(f$assign_rhs))
    },
    Tilde = {
      f <- sexp_inline(s)
      kind <- ast_fun_kind(f$kind)
      list(node = "Tilde", arg = ast_expr(f$arg), dist = ast_identifier_name(f$distribution),
           origin = kind$origin, suffix = kind$suffix, args = lapply(f$args, ast_expr),
           truncation = ast_truncation(f$truncation))
    },
    IfThenElse = list(node = "If", cond = ast_expr(s[[2]]), then = ast_stmt(s[[3]]),
                      otherwise = if (length(s[[4]]) == 0) NULL else ast_stmt(s[[4]][[1]])),
    While = list(node = "While", cond = ast_expr(s[[2]]), body = ast_stmt(s[[3]])),
    For = {
      f <- sexp_inline(s)
      list(node = "For", var = ast_identifier_name(f$loop_variable), lower = ast_expr(f$lower_bound),
           upper = ast_expr(f$upper_bound), body = ast_stmt(f$loop_body))
    },
    ForEach = list(node = "ForEach", var = ast_identifier_name(s[[2]]), iter = ast_expr(s[[3]]),
                   body = ast_stmt(s[[4]])),
    Profile = list(node = "Profile", name = sexp_unquote_stan(s[[2]]), stmts = lapply(s[[3]], ast_stmt)),
    Block = list(node = "Block", stmts = lapply(s[[2]], ast_stmt)),
    VarDecl = {
      f <- sexp_inline(s)
      vars <- lapply(f$variables, function(v) {
        vr <- sexp_record(v)
        list(name = ast_identifier_name(vr$identifier),
             init = if (length(vr$initial_value)) ast_expr(vr$initial_value[[1]]) else NULL)
      })
      list(node = "Decl", stype = ast_stype(f$decl_type),
           transform = ast_transformation(f$transformation), vars = vars)
    },
    FunDef = {
      f <- sexp_inline(s)
      rt <- f$returntype
      returns <- if (identical(rt, "Void")) NULL else ast_utype(rt[[2]])
      args <- lapply(f$arguments, function(a) {
        # (AutoDiffable UReal ((name x) (id_loc <opaque>))) ; ad level may be DataOnly
        list(ad = a[[1]], type = ast_utype(a[[2]]), name = ast_identifier_name(a[[3]]))
      })
      body <- ast_stmt(f$body)
      list(node = "FunDef", name = ast_identifier_name(f$funname), returns = returns,
           args = args, body = body, forward = identical(body$node, "Skip"))
    },
    ast_fail("statement", s))
}

ast_block_names <- c(functionblock = "functions", datablock = "data",
                     transformeddatablock = "transformed data", parametersblock = "parameters",
                     transformedparametersblock = "transformed parameters", modelblock = "model",
                     generatedquantitiesblock = "generated quantities")

# program -> list of blocks, each NULL (absent) or list of IR statements
ast_program <- function(sexp) {
  rec <- sexp_record(sexp)
  missing_fields <- setdiff(names(ast_block_names), names(rec))
  if (length(missing_fields)) ast_fail("program field", missing_fields[1])
  blocks <- lapply(names(ast_block_names), function(field) {
    b <- rec[[field]]
    if (length(b) == 0) return(NULL)
    stmts <- sexp_record(b[[1]])$stmts
    lapply(stmts, ast_stmt)
  })
  names(blocks) <- unname(ast_block_names)
  blocks
}
