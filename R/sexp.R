# Reader for the s-expressions printed by `stanc --debug-decorated-ast`
# (stanc.js option `debug-typed-ast`). Atoms become character scalars and lists
# become unnamed R lists, so `(a (b c))` reads as list("a", list("b", "c")).

sexp_tokenize <- function(text) {
  pattern <- '\\(|\\)|"(?:[^"\\\\]|\\\\.)*"|[^\\s()"]+'
  regmatches(text, gregexpr(pattern, text, perl = TRUE))[[1]]
}

sexp_unquote <- function(atom) {
  if (!startsWith(atom, "\"")) return(atom)
  inner <- substr(atom, 2, nchar(atom) - 1)
  # sexplib escapes: \" \\ \n \t; decode in one pass so "\\\\n" stays literal
  m <- gregexpr("\\\\(.)", inner, perl = TRUE)
  regmatches(inner, m) <- lapply(regmatches(inner, m), function(esc) {
    ch <- substr(esc, 2, 2)
    out <- ch
    out[ch == "n"] <- "\n"
    out[ch == "t"] <- "\t"
    out
  })
  inner
}

read_sexp <- function(text) {
  tokens <- sexp_tokenize(text)
  if (length(tokens) == 0) stop("Stan2R: empty s-expression from stanc")
  stack <- vector("list", 64)
  depth <- 0L
  current <- list()
  n_items <- 0L
  result <- NULL
  for (tok in tokens) {
    if (tok == "(") {
      depth <- depth + 1L
      if (depth > length(stack)) length(stack) <- 2L * length(stack)
      stack[[depth]] <- list(items = current, n = n_items)
      current <- list()
      n_items <- 0L
    } else if (tok == ")") {
      if (depth == 0L) stop("Stan2R: unbalanced ')' in stanc s-expression")
      finished <- current
      parent <- stack[[depth]]
      depth <- depth - 1L
      current <- parent$items
      n_items <- parent$n + 1L
      current[n_items] <- list(finished)
      if (depth == 0L) result <- c(result, list(finished))
    } else {
      n_items <- n_items + 1L
      current[[n_items]] <- sexp_unquote(tok)
    }
  }
  if (depth != 0L) stop("Stan2R: unbalanced '(' in stanc s-expression")
  if (length(result) != 1L) stop("Stan2R: expected one top-level s-expression from stanc")
  result[[1]]
}
