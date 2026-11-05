# Source map: line spans and comments for every statement.
# stanc's AST dump hides locations (`<opaque>`), so the original text is lexed
# here and scanned for statement boundaries only (no expression parsing). The
# scanned statements are then matched one-to-one with the AST statements; any
# disagreement is an internal error rather than a silently shifted line number.

stan_lex <- function(text) {
  pattern <- paste0(
    "//[^\\n]*", "|/\\*[\\s\\S]*?\\*/", "|#include[^\\n]*",
    '|"(?:[^"\\\\\\n]|\\\\.)*"',
    "|[A-Za-z_][A-Za-z0-9_]*",
    "|(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:[eE][+-]?[0-9]+)?i?",
    "|\\.\\*=|\\./=|%/%|\\.\\*|\\./|\\.\\^|\\+=|-=|\\*=|/=|==|!=|<=|>=|&&|\\|\\||<-",
    "|\\S")
  m <- gregexpr(pattern, text, perl = TRUE)[[1]]
  if (m[1] == -1) {
    return(data.frame(text = character(0), start = integer(0), end = integer(0),
                      line = integer(0), end_line = integer(0), col = integer(0),
                      kind = character(0), stringsAsFactors = FALSE))
  }
  starts <- as.integer(m)
  ends <- starts + attr(m, "match.length") - 1L
  toks <- regmatches(text, list(m))[[1]]
  newlines <- gregexpr("\n", text, fixed = TRUE)[[1]]
  newlines <- if (newlines[1] == -1) integer(0) else as.integer(newlines)
  line_of <- function(pos) findInterval(pos - 1L, newlines) + 1L
  line_start <- c(1L, newlines + 1L)
  line <- line_of(starts)
  kind <- ifelse(startsWith(toks, "//") | startsWith(toks, "/*"), "comment",
                 ifelse(startsWith(toks, "#include"), "include",
                        ifelse(startsWith(toks, "\""), "string",
                               ifelse(grepl("^[A-Za-z_]", toks), "ident", "other"))))
  data.frame(text = toks, start = starts, end = ends, line = line,
             end_line = line_of(ends), col = starts - line_start[line] + 1L,
             kind = kind, stringsAsFactors = FALSE)
}

stan_block_headers <- list(
  "functions" = "functions", "data" = "data", "transformed data" = c("transformed", "data"),
  "parameters" = "parameters", "transformed parameters" = c("transformed", "parameters"),
  "model" = "model", "generated quantities" = c("generated", "quantities"))

# Scan tokens (comments removed) into a statement tree with token spans.
stan_scan <- function(tokens) {
  code <- tokens[tokens$kind != "comment", , drop = FALSE]
  if (any(code$kind == "include")) {
    stop("Stan2R: internal error, #include must be expanded before scanning", call. = FALSE)
  }
  tx <- code$text
  n <- length(tx)
  pos <- 1L

  scan_error <- function(msg) {
    where <- if (pos <= n) sprintf(" at line %d", code$line[pos]) else " at end of file"
    stop("Stan2R: could not scan statement structure", where, ": ", msg, call. = FALSE)
  }
  expect <- function(tok) {
    if (pos > n || tx[pos] != tok) scan_error(sprintf("expected '%s'", tok))
    pos <<- pos + 1L
  }
  skip_balanced <- function(open, close) {
    expect(open)
    depth <- 1L
    while (depth > 0L) {
      if (pos > n) scan_error(sprintf("unbalanced '%s'", open))
      if (tx[pos] == open) depth <- depth + 1L
      if (tx[pos] == close) depth <- depth - 1L
      pos <<- pos + 1L
    }
  }
  # advance to the first token in `stops` at bracket depth 0
  scan_to <- function(stops) {
    depth <- 0L
    while (pos <= n) {
      t <- tx[pos]
      if (depth == 0L && t %in% stops) return(t)
      if (t %in% c("(", "[", "{")) depth <- depth + 1L
      if (t %in% c(")", "]", "}")) depth <- depth - 1L
      pos <<- pos + 1L
    }
    scan_error(sprintf("expected one of %s", paste(stops, collapse = " ")))
  }
  node <- function(kind, first, ...) list(kind = kind, first = first, last = pos - 1L, ...)

  stmt <- function(in_functions = FALSE) {
    first <- pos
    t <- tx[pos]
    if (in_functions) {
      stop_tok <- scan_to(c("{", ";"))
      if (stop_tok == ";") {
        pos <<- pos + 1L
        return(node("fundef", first, body = NULL))
      }
      body <- stmt()
      return(node("fundef", first, body = body))
    }
    if (t == "{") {
      pos <<- pos + 1L
      kids <- list()
      while (pos <= n && tx[pos] != "}") kids[[length(kids) + 1L]] <- stmt()
      close <- pos
      expect("}")
      return(node("block", first, children = kids, close = close))
    }
    if (t %in% c("for", "while")) {
      pos <<- pos + 1L
      skip_balanced("(", ")")
      body <- stmt()
      return(node(t, first, body = body))
    }
    if (t == "if") {
      pos <<- pos + 1L
      skip_balanced("(", ")")
      then <- stmt()
      otherwise <- NULL
      if (pos <= n && tx[pos] == "else") {
        pos <<- pos + 1L
        otherwise <- stmt()
      }
      return(node("if", first, then = then, otherwise = otherwise))
    }
    if (t == "profile") {
      pos <<- pos + 1L
      skip_balanced("(", ")")
      expect("{")
      kids <- list()
      while (pos <= n && tx[pos] != "}") kids[[length(kids) + 1L]] <- stmt()
      close <- pos
      expect("}")
      return(node("profile", first, children = kids, close = close))
    }
    if (t == ";") {
      pos <<- pos + 1L
      return(node("skip", first))
    }
    scan_to(";")
    pos <<- pos + 1L
    node("simple", first)
  }

  blocks <- list()
  while (pos <= n) {
    first <- pos
    name <- NULL
    for (b in names(stan_block_headers)) {
      words <- stan_block_headers[[b]]
      k <- length(words)
      if (pos + k - 1L <= n && all(tx[pos:(pos + k - 1L)] == words) &&
          pos + k <= n && tx[pos + k] == "{") {
        name <- b
        pos <- pos + k
        break
      }
    }
    if (is.null(name)) scan_error(sprintf("unexpected '%s' outside a program block", tx[pos]))
    expect("{")
    kids <- list()
    while (pos <= n && tx[pos] != "}") kids[[length(kids) + 1L]] <- stmt(in_functions = name == "functions")
    close <- pos
    expect("}")
    blocks[[name]] <- list(kind = "program_block", first = first, last = pos - 1L,
                           children = kids, close = close)
  }
  list(blocks = blocks, code = code)
}

# Attach line spans, comments, and blank-line flags to IR statements.
attach_source_map <- function(blocks, stan_lines) {
  text <- paste(stan_lines, collapse = "\n")
  tokens <- stan_lex(text)
  scanned <- stan_scan(tokens)
  code <- scanned$code
  comments <- tokens[tokens$kind == "comment", , drop = FALSE]
  comments$claimed <- rep(FALSE, nrow(comments))

  # blank line before each token/comment, in source order
  all_items <- tokens[order(tokens$start), ]
  prev_end <- c(0L, all_items$end_line[-nrow(all_items)])
  all_items$blank_before <- all_items$line - prev_end > 1L & seq_len(nrow(all_items)) > 1L
  blank_at <- setNames(all_items$blank_before, all_items$start)

  comment_text <- function(i) {
    raw <- comments$text[i]
    if (startsWith(raw, "//")) return(sub("^//\\s?", "", raw))
    body <- sub("^/\\*+", "", sub("\\*+/$", "", raw))
    body <- strsplit(body, "\n", fixed = TRUE)[[1]]
    body <- trimws(sub("^\\s*\\*", "", body))
    body <- body[nzchar(body) | seq_along(body) %in% range(seq_along(body))]
    body[nzchar(body) | length(body) == 1L]
  }
  take_before <- function(offset) {
    idx <- which(!comments$claimed & comments$start < offset)
    comments$claimed[idx] <<- TRUE
    lapply(idx, function(i) list(text = comment_text(i), blank_before = unname(blank_at[as.character(comments$start[i])])))
  }
  take_trailing <- function(last_tok) {
    end_line <- code$end_line[last_tok]
    next_start <- if (last_tok < nrow(code)) code$start[last_tok + 1L] else Inf
    idx <- which(!comments$claimed & comments$start > code$end[last_tok] &
                   comments$line == end_line & comments$start < next_start)
    if (length(idx) == 0) return(NULL)
    comments$claimed[idx] <<- TRUE
    unlist(lapply(idx, comment_text))
  }

  mismatch <- function(where, ir_kind, scan_kind) {
    stop(sprintf("Stan2R: internal error aligning source with stanc AST in %s (AST %s vs source %s).",
                 where, ir_kind, scan_kind), call. = FALSE)
  }
  scan_kind_of <- function(ir) {
    switch(ir$node, Block = "block", For = , ForEach = "for", While = "while", If = "if",
           Profile = "profile", Skip = "skip", FunDef = "fundef", "simple")
  }

  align <- function(ir, sc, where) {
    if (scan_kind_of(ir) != sc$kind) mismatch(where, ir$node, sc$kind)
    start_off <- code$start[sc$first]
    ir$comments_before <- take_before(start_off)
    ir$blank_before <- unname(blank_at[as.character(start_off)])
    ir$line <- code$line[sc$first]
    ir$end_line <- code$end_line[sc$last]
    ir$col <- code$col[sc$first]
    ir$span <- c(start = start_off, end = code$end[sc$last])
    if (ir$node %in% c("Block", "Profile")) {
      if (length(ir$stmts) != length(sc$children)) mismatch(where, paste(length(ir$stmts), "statements"), paste(length(sc$children), "statements"))
      ir$stmts <- Map(align, ir$stmts, sc$children, where)
      ir$closing_comments <- take_before(code$start[sc$close])
    } else if (ir$node %in% c("For", "ForEach", "While")) {
      ir$body <- align(ir$body, sc$body, where)
    } else if (ir$node == "If") {
      ir$then <- align(ir$then, sc$then, where)
      if (xor(is.null(ir$otherwise), is.null(sc$otherwise))) mismatch(where, "if/else", "if/else")
      if (!is.null(ir$otherwise)) ir$otherwise <- align(ir$otherwise, sc$otherwise, where)
    } else if (ir$node == "FunDef") {
      if (!ir$forward) {
        if (is.null(sc$body)) mismatch(where, "function body", "forward declaration")
        ir$body <- align(ir$body, sc$body, where)
      }
    }
    if (!(ir$node %in% c("Block", "Profile"))) {
      ir$comment_after <- take_trailing(sc$last)
    } else {
      ir$comment_after <- take_trailing(sc$close)
    }
    ir
  }

  present <- names(blocks)[!vapply(blocks, is.null, logical(1))]
  if (!setequal(present, names(scanned$blocks))) {
    stop("Stan2R: internal error, program blocks found in the source (",
         paste(names(scanned$blocks), collapse = ", "), ") differ from stanc's (",
         paste(present, collapse = ", "), ").", call. = FALSE)
  }
  block_info <- list()
  for (b in present) {
    sc <- scanned$blocks[[b]]
    if (length(blocks[[b]]) != length(sc$children)) {
      mismatch(paste(b, "block"), paste(length(blocks[[b]]), "statements"), paste(length(sc$children), "statements"))
    }
    leading <- take_before(code$start[sc$first])
    blocks[[b]] <- Map(align, blocks[[b]], sc$children, paste(b, "block"))
    block_info[[b]] <- list(comments_before = leading,
                            closing_comments = take_before(code$start[sc$close]),
                            line = code$line[sc$first], end_line = code$line[sc$close],
                            span = c(start = code$start[sc$first], end = code$end[sc$close]),
                            open_brace = code$end[sc$first + length(stan_block_headers[[b]])],
                            close_brace = code$start[sc$close])
  }
  list(blocks = blocks, block_info = block_info,
       trailing_comments = take_before(Inf), text = text)
}
