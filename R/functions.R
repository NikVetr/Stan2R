#### general functions ####

# subset_samps <- function(var_name, samps) {
#   # Regex pattern to match exact variable name and indexed names
#   pattern <- paste0("^", var_name, "(\\.\\d+)*(\\.?)$")
#   matched_cols <- grep(pattern, colnames(samps), value = TRUE)
#   return(samps[, matched_cols, drop = FALSE])
# }

subset_samps <- function(var_name, samps) {
  pattern <- paste0("^", var_name, "(\\.\\d+)*(\\.?)$")
  matched_cols <- grep(pattern, names(samps), value = TRUE)
  return(samps[, ..matched_cols])
}

munge_samps <- function(var_name, df) {
  # Check if the data frame has only one row
  if (nrow(df) == 1) {
    # Process a single row to construct the appropriate data structure
    # Determine the structure of the data (scalar, vector, matrix, or array)
    col_names <- colnames(df)
    col_names_without_var <- gsub(paste0("^", var_name, "\\."), "", col_names)
    if (all(grepl("\\.\\d+\\.", col_names))) {
      # Handle vectors, matrices, and arrays
      # Extract indices from column names
      indices <- lapply(strsplit(col_names_without_var, "\\."), function(x) as.numeric(x))
      
      # Determine if it's a vector, matrix, or array
      max_dims <- sapply(indices, length)
      if (all(max_dims == 1)) {
        # Vector
        return(unlist(df))
      } else if (all(max_dims == 2)) {
        # Matrix
        dim_order <- do.call(rbind, indices)
        mat <- matrix(NA, nrow = max(dim_order[,1]), ncol = max(dim_order[,2]))
        mat[cbind(dim_order[,1], dim_order[,2])] <- unlist(df)
        return(mat)
      } else {
        # Array (higher dimensions)
        array_dims <- apply(do.call(rbind, indices), 2, max)
        arr <- array(NA, dim = array_dims)
        array_indices <- do.call(expand.grid, lapply(array_dims, seq_len))
        arr[as.matrix(array_indices)] <- unlist(df)
        return(arr)
      }
    } else {
      # Scalar or simple vector
      return(unlist(df))
    }
  } else {
    # Recursively apply to each row
    return(lapply(seq_len(nrow(df)), function(i) munge_samps(var_name = var_name, df = df[i, , drop = FALSE])))
  }
}

clean_stan <- function(stan_code, return_annot = F, R_style = T) {
  
  # Regular expressions for detecting types of lines
  declaration_regex <- "^\\s*(int|real|vector|row_vector|matrix|array)"
  lpmf_lpdf_regex <- "_lp(mf|df)\\s*\\("
  target_assignment_regex <- "^\\s*target\\s*(\\+=|=)"
  
  # create a data object to modify
  annot <- matrix("", nrow = length(stan_code), ncol = 5)
  colnames(annot) <- c("leading_ws", "trailing_ws", "comments", "preceding_content", "n_empty_preceding")
  annot <- data.frame(annot)
  annot$n_empty_preceding <- 0
  lines <- stan_code
  
  # Remove comments (text after //)
  has_comments <- grepl("//", lines)
  annot$comments[has_comments] <- paste0(ifelse(R_style, "#", "//"), sub(".*?// ?", "", lines[has_comments]))
  lines <- gsub("//.*", "", lines)
  
  #trim whitespace
  lws <- regmatches(lines, gregexpr("^\\s+", lines))
  lws[sapply(lws, length) == 0] <- ""
  tws <- regmatches(lines, gregexpr("\\s+$", lines))
  tws[sapply(tws, length) == 0] <- ""
  annot$leading_ws <- unlist(lws)
  annot$trailing_ws <- unlist(tws)
  lines <- trimws(lines)
  
  #remove one chunk of the preceding whitespace to reflect different style conventions
  if(R_style){
    annot$leading_ws <- gsub("\t", "  ", annot$leading_ws)
    smallest_lws <- nchar(annot$leading_ws)
    smallest_lws <- paste0(rep(x = " ", min(smallest_lws[smallest_lws > 0])), collapse = "")
    annot$leading_ws <- sub(smallest_lws, "", annot$leading_ws)  
  }
  
  # Remove empty lines
  empty_lines <- nchar(lines) == 0
  rleel <- rle(empty_lines)
  el_blocks <- cumsum(rleel$lengths)
  el_nextline_inds <- el_blocks[rleel$values] + 1
  el_preceding_len <- rleel$lengths[rleel$values]
  annot$n_empty_preceding[el_nextline_inds] <- el_preceding_len
  annot$preceding_content[el_nextline_inds] <- sapply(el_nextline_inds, function(i){
    #if previous line is just whitespace sans comments, we do not need to have ws, just \n?
    prev_lines_content <- annot[(i-annot$n_empty_preceding[i]):(i-1), 
                                c("leading_ws", "comments")]
    prev_lines_content <- apply(prev_lines_content, 1, function(x) 
      paste0(unlist(x), collapse = ""))
    prev_lines_content <- paste0(paste0(prev_lines_content, "\n"), collapse = "")
    prev_lines_content
  })
  annot <- annot[!empty_lines, c("leading_ws", "trailing_ws", 
                                 "comments", "preceding_content")]
  lines <- lines[!empty_lines]
  
  # Ensure closing brackets are on their own lines
  straggling_close_brackets <- grepl("}", lines) & nchar(lines) != 1
  if(any(straggling_close_brackets)){
    close_bracket_inds <- which(straggling_close_brackets)
    n_lines_inserted <- 0
    for(i in 1:sum(straggling_close_brackets)){
      current_line <- (close_bracket_inds[i] + n_lines_inserted)
      insert_lines <- strsplit(x = lines[current_line], 
                               split = paste0("(?=", "}", ")"), 
                               perl = TRUE)[[1]]
      
      #insert the lines
      lines <- append(lines, insert_lines, current_line)
      insert_annot <- rbind(annot[current_line,], 
                            matrix("", length(insert_lines)-1, ncol(annot), 
                                   dimnames = list(NULL, colnames(annot)))
      )
      insert_annot$leading_ws <- insert_annot$leading_ws[1]
      insert_annot$trailing_ws <- insert_annot$trailing_ws[1]
      annot <- rbind(annot[1:current_line,], 
                     insert_annot,
                     annot[(current_line+1):nrow(annot),])
      
      #and now delete the old lines
      lines <- lines[-current_line]
      annot <- annot[-current_line,]
      n_lines_inserted <- n_lines_inserted + length(insert_lines) - 1
    }
  }
  
  # Ensure opening brackets are not followed by anything
  # ie only one opening bracket per line
  open_bracket_lines <- grep("\\{", lines)
  straggling_open_bracket_inds <- open_bracket_lines[
    sapply(open_bracket_lines, function(i){
      sum(gregexpr("\\{", lines[i])[[1]] > 0) > 1
    })]
  
  if(length(straggling_open_bracket_inds) > 0){
    n_lines_inserted <- 0
    for(i in 1:length(straggling_open_bracket_inds)){
      current_line <- (straggling_open_bracket_inds[i] + n_lines_inserted)
      insert_lines <- trimws(strsplit(x = lines[current_line], 
                                      split = paste0("(?=", "{", ")"), 
                                      perl = TRUE)[[1]])
      
      if(insert_lines[1] == "{"){
        insert_lines <- c("", insert_lines)
      }
      
      last_bit_exists <- F
      if(insert_lines[length(insert_lines)] != "{"){
        last_bit <- insert_lines[length(insert_lines)]
        last_bit_exists <- T
      }
      
      insert_lines <- trimws(paste(insert_lines[which(insert_lines == "{")-1],
                                   insert_lines[which(insert_lines == "{")]
      ))
      
      if(last_bit_exists){
        insert_lines <- c(insert_lines, last_bit)
      }
      
      #insert the lines
      lines <- append(lines, insert_lines, current_line)
      
      #modify the annotation
      insert_annot <- rbind(annot[current_line,], 
                            matrix("", length(insert_lines)-1, ncol(annot), 
                                   dimnames = list(NULL, colnames(annot)))
                            )
      insert_annot$leading_ws <- insert_annot$leading_ws[1]
      insert_annot$trailing_ws <- insert_annot$trailing_ws[1]
      annot <- rbind(annot[1:current_line,], 
                     insert_annot,
                     annot[(current_line+1):nrow(annot),])
      
      #and now delete the old lines
      lines <- lines[-current_line]
      annot <- annot[-current_line,]
      n_lines_inserted <- n_lines_inserted + length(insert_lines) - 1
     
    }
  }
  
  # get multiline statements on the same line
  # if it's not a "} or ending in a "{" or ";",
  # concatenate lines up to the most recent line after a ";"
  continuation_lines <- !(substr(lines, nchar(lines), nchar(lines)) %in% c("{", ";", "}"))
  if(any(continuation_lines)){
    rlel <- rle(continuation_lines)
    rle_df <- data.frame(val = rlel$values, 
                         start = cumsum(c(0, rlel$lengths)[-length(rlel$lengths)]) + 1,
                         stop = cumsum(rlel$lengths))
    rle_df <- rle_df[rle_df$val,]
    n_lines_removed <- 0
    for(i in 1:nrow(rle_df)){
      line_bounds <- unlist(rle_df[i,2:3] + c(0,1)) - n_lines_removed
      remove_lines <- lines[line_bounds[1]:line_bounds[2]]
      insert_line <- paste0(remove_lines, collapse = " ")
      
      #insert the line
      lines <- append(lines, insert_line, line_bounds[2])
      
      #propagate later comments and preceding content
      insert_annot <- annot[line_bounds[1],]
      insert_annot_comments <- annot$comments[line_bounds[1]:line_bounds[2]]
      insert_annot_comments[sapply(trimws(insert_annot_comments), nchar) > 0] <- paste0(
        insert_annot_comments[sapply(trimws(insert_annot_comments), nchar) > 0], " "
      )
      insert_annot$comments <- paste0(insert_annot_comments, collapse = "")
      insert_annot_preceding_content <- annot$preceding_content[line_bounds[1]:line_bounds[2]]
      # insert_annot_preceding_content[sapply(trimws(insert_annot_preceding_content), nchar) > 0]
      insert_annot$preceding_content <- paste0(insert_annot_preceding_content, collapse = "")
      annot <- rbind(annot[1:line_bounds[2],], 
                     insert_annot,
                     annot[(line_bounds[2]+1):nrow(annot),])
      
      #and now delete the old lines
      lines <- lines[-(line_bounds[1]:line_bounds[2])]
      annot <- annot[-(line_bounds[1]:line_bounds[2]),]
      n_lines_removed <- n_lines_removed + length(remove_lines) - 1
    }
  }
  
  #split variable declaration and assignment
  assignment_and_declaration_inds <- which(grepl(declaration_regex, lines) & 
                                             grepl("=", gsub("<[^>]*>", "BOUND_PLACEHOLDER", lines)))
  if(length(assignment_and_declaration_inds) > 0){
    n_lines_inserted <- 0
    for(i in 1:length(assignment_and_declaration_inds)){
      current_line <- assignment_and_declaration_inds[i] + n_lines_inserted
      line <- lines[current_line]
      
      rhs <- paste0(strsplit(gsub("<[^>]*>", "BOUND_PLACEHOLDER", line), "(?==)", perl = T)[[1]][-1], collapse = "")
      lhs <- trimws(substr(line, 1, nchar(line) - nchar(rhs)))
      insert_lines <- c(paste0(lhs, ";"), 
                        paste(tail(strsplit(lhs, " ")[[1]], 1), rhs))
      
      #insert the lines
      lines <- append(lines, insert_lines, current_line)
      insert_annot <- rbind(annot[current_line,], 
                            matrix("", length(insert_lines)-1, ncol(annot), 
                                   dimnames = list(NULL, colnames(annot)))
      )
      insert_annot$leading_ws <- insert_annot$leading_ws[1]
      insert_annot$trailing_ws <- insert_annot$trailing_ws[1]
      annot <- rbind(annot[1:current_line,], 
                     insert_annot,
                     annot[(current_line+1):nrow(annot),])
      
      #and now delete the old lines
      lines <- lines[-current_line]
      annot <- annot[-current_line,]
      
      n_lines_inserted <- n_lines_inserted + length(insert_lines) - 1
    }
    annot$curr_line <- lines
  }
  
  if(return_annot){
    return(list(code = lines, annot = annot[,c("leading_ws", "trailing_ws", "comments", "preceding_content")]))  
  } else {
    return(lines)
  }
  
}

#### parsing functions ####

parse_stan_blocks <- function(stan_code, return_annot_mod = F) {
  
  # Split the code into lines
  lines <- stan_code
  
  # Initialize variables
  block_mod <- integer(0)
  current_block <- NULL
  blocks <- list()
  bracket_count <- 0  # Counter for open curly brackets
  block_regex <- "^\\s*(data|transformed data|parameters|transformed parameters|model|functions|generated quantities)\\s*\\{"
  
  # Iterate over the lines
  for (line_i in seq_along(lines)) {
    line <- lines[line_i]
    
    # Check for block start
    if (grepl(block_regex, line)) {
      current_block <- sub("^\\s*(\\w+(\\s+\\w+)?)\\s*\\{.*", "\\1", line)
      blocks[[current_block]] <- character(0)
      bracket_count <- 1
      block_mod <- c(block_mod, line_i)
      
    } else if (!is.null(current_block)) {
      
      # Check for open and close curly brackets
      bracket_count <- bracket_count + 
        sum(gregexpr("\\{", line)[[1]] > 0) - 
        sum(gregexpr("\\}", line)[[1]] > 0)
      
      if (bracket_count == 0) {
        current_block <- NULL
        block_mod <- c(block_mod, line_i)
      } else {
        # If we are inside a block, add the line to the current block
        blocks[[current_block]] <-  c(blocks[[current_block]], line)
      }
    }
  }
  
  #return objects per request
  if(return_annot_mod){
    return(list(block_contents = blocks, block_mod = block_mod))  
  } else {
    return(blocks)
  }
  
}

parse_declaration <- function(line) {
  # Initialize variables to store extracted information
  declaration_type <- ""
  var_name <- ""
  bounds <- ""
  dimensions <- ""
  distribution <- NA  # Not applicable for declaration lines
  parameters <- NA    # Not applicable for declaration lines
  
  # Regular expressions for different variable types
  # Updated regex to ensure bounds are captured correctly without angle brackets
  scalar_regex <- "^(int|real)\\s*(<([^>]*)>)?\\s+(\\w+);$"
  vector_regex <- "^(vector|row_vector)(<([^>]*)>)?\\[(.+)\\]\\s+(\\w+);$"
  matrix_regex <- "^matrix(<([^>]*)>)?\\[(.+)\\]\\s+(\\w+);$"
  array_regex <- "^array\\[(.+)\\]\\s+(int|real)\\s*(<([^>]*)>)?\\s+(\\w+);$"
  
  # Check and extract information for scalar types
  if (grepl(scalar_regex, line)) {
    declaration_type <- sub(scalar_regex, "\\1", line)
    bounds <- sub(scalar_regex, "\\3", line)  # Captures contents within angle brackets
    var_name <- sub(scalar_regex, "\\4", line)
    dimensions <- "1"  # Scalars have a dimension of 1
  } else if (grepl(vector_regex, line)) {  # Vector types
    declaration_type <- sub(vector_regex, "\\1", line) 
    bounds <- sub(vector_regex, "\\3", line)  # Updated to capture bounds correctly          
    dimensions <- sub(vector_regex, "\\4", line)      
    var_name <- sub(vector_regex, "\\5", line)        
  } else if (grepl(matrix_regex, line)) { 
    declaration_type <- "matrix"
    bounds <- sub(matrix_regex, "\\2", line)  
    dimensions <- sub(matrix_regex, "\\3", line)
    var_name <- sub(matrix_regex, "\\4", line) 
  } else if (grepl(array_regex, line)) {
    declaration_type <- sub(array_regex, "\\2", line)  
    dimensions <- sub(array_regex, "\\1", line) 
    bounds <- sub(array_regex, "\\4", line) 
    var_name <- sub(array_regex, "\\5", line)
  }
  
  #add in explicit bounds if necessary
  if (bounds == "") {
    bounds <- "lower=-Inf, upper=Inf"
  } else if (!grepl("upper", bounds)) {
    bounds <- paste0(bounds, ", upper=Inf")
  } else if (!grepl("lower", bounds)) {
    bounds <- paste0("lower=-Inf, ", bounds)
  }
  
  return(list(declaration_type = declaration_type, var_name = var_name, bounds = bounds, dimensions = dimensions, distribution = distribution, parameters = parameters))
}

parse_assignment <- function(line) {
  # Initialize variables to store extracted information
  lhs_name <- ""
  lhs_index <- NA
  rhs_expression <- ""
  
  # Regular expression for assignment statements
  # This regex captures the LHS (including index), and the RHS expression
  assignment_regex <- "^(\\w+)(\\[.*\\])?\\s*=\\s*(.+);$"
  
  # Check if the line matches the assignment statement pattern
  if (grepl(assignment_regex, line)) {
    lhs_name <- sub(assignment_regex, "\\1", line)
    lhs_index <- sub(assignment_regex, "\\2", line)
    rhs_expression <- sub(assignment_regex, "\\3", line)
    
    # Removing any additional brackets around the index
    lhs_index <- gsub("\\[|\\]", "", lhs_index)
    
    return(list(lhs_name = lhs_name, lhs_index = lhs_index, 
                rhs_expression = rhs_expression))
    
  }
  
  return(list(lhs_name = lhs_name, lhs_index = lhs_index, rhs_expression = rhs_expression))
}

parse_lse <- function(line) {
  # Regular expression to capture LHS and RHS in 'target += log_sum_exp(...);' statements
  lse_regex <- "^(\\w+)(\\[.*\\])?\\s*\\+=\\s*log_sum_exp\\((.+)\\);$"
  
  if (grepl(lse_regex, line)) {
    lhs_name <- sub(lse_regex, "\\1", line)
    lse_argument <- sub(lse_regex, "\\3", line)
    
    return(list(lhs_name = lhs_name, lse_argument = lse_argument))
  }
  
  return(list(lhs_name = NA, lse_argument = NA))
}

parse_operation <- function(expression) {
  # Replace Stan-specific operations with R equivalents
  expression <- gsub("\\.\\*", "*", expression)  # Element-wise multiplication
  expression <- gsub("\\./", "/", expression)    # Element-wise division
  
  
  # Handle some of Stan's built-in functions, translating them to R equivalents
  # Example: Stan's normal() becomes R's dnorm()
  # Add similar replacements for other Stan functions as needed
  expression <- gsub("normal\\(", "dnorm(", expression)
  expression <- gsub("bernoulli_logit\\(", "dbinom(", expression)
  
  # ... [Additional replacements for other Stan-specific functions and operations]
  
  return(expression)
}

parse_loop <- function(line) {
  loop_type <- NA
  loop_variable <- NA
  loop_range <- NA
  loop_action <- if (grepl("^\\s*\\}", line)) "close" else "open"
  
  # Handle opening of loops
  if (loop_action == "open") {
    for_loop_regex <- "^for\\s*\\((\\w+)\\s+in\\s+(.+)\\)\\s*\\{"
    while_loop_regex <- "^while\\s*\\((.+)\\)\\s*\\{"
    
    if (grepl(for_loop_regex, line)) {
      loop_type <- "for"
      loop_variable <- sub(for_loop_regex, "\\1", line)
      loop_range <- sub(for_loop_regex, "\\2", line)
    } else if (grepl(while_loop_regex, line)) {
      loop_type <- "while"
      loop_range <- sub(while_loop_regex, "\\1", line)
    }
  }
  
  return(list(loop_type = loop_type, loop_variable = loop_variable, loop_range = loop_range, loop_action = loop_action))
}

parse_sampling <- function(line) {
  # Initialize variables to store extracted information
  variable_name <- ""
  distribution <- ""
  parameters <- ""
  declaration_type <- NA  # Not applicable for sampling lines
  bounds <- NA           # Not applicable for sampling lines
  dimensions <- NA       # Not applicable for sampling lines
  
  # Regular expression for sampling statements
  # sampling_regex <- "^(\\w+)\\s*~\\s*(\\w+)\\(([^)]*)\\);"
  sampling_regex <- "^(.+?)\\s*~\\s*(\\w+)\\(([^)]*)\\);"
  
  # Check if the line matches the sampling statement pattern
  if (grepl(sampling_regex, line)) {
    variable_name <- sub(sampling_regex, "\\1", line)
    distribution <- sub(sampling_regex, "\\2", line)
    parameters <- sub(sampling_regex, "\\3", line)
  }
  
  return(list(declaration_type = declaration_type, var_name = variable_name, bounds = bounds, dimensions = dimensions, distribution = distribution, parameters = parameters))
}

parse_lpmf_lpdf <- function(line) {
  # Extract the RHS of the assignment
  rhs <- sub(".*=\\s*(.+);", "\\1", line)
  lhs <- sub("\\s*=.*", "", line)
  
  # Isolate the _lpmf or _lpdf function call
  lpmf_lpdf_call <- regmatches(rhs, regexpr("\\w+_lp(mf|df)\\s*\\([^\\)]+\\)", rhs))
  
  # Initialize variables for extracted information
  distribution <- ""
  outcome_variable <- ""
  parameters <- ""
  mixture_probability <- NA  # Default to 0
  
  # Check if there's a mixture probability component
  if (grepl("\\+", rhs)) {
    mixture_probability <- sub("\\s*\\+\\s*\\w+_lp(mf|df).*", "", rhs)
    mixture_probability <- sub("log\\s*\\((.+)\\)", "\\1", mixture_probability)
  }
  
  # Process the _lpmf or _lpdf call
  if (length(lpmf_lpdf_call) > 0) {
    # Extract the distribution
    distribution <- sub("(.+)_lp(mf|df).*", "\\1", lpmf_lpdf_call)
    
    # Extract the outcome variable and parameters
    args <- sub(".*\\(([^|]+)\\s*\\|\\s*(.+)\\)", "\\1,\\2", lpmf_lpdf_call)
    args_split <- strsplit(args, ",")[[1]]
    if (length(args_split) >= 2) {
      outcome_variable <- trimws(args_split[1])
      parameters <- paste(args_split[-1], collapse = ", ")
    }
  }
  
  return(list(distribution = distribution, outcome_variable = outcome_variable, 
              parameters = parameters, mixture_probability = mixture_probability,
              lhs = lhs, rhs = rhs))
}



parse_target_assignment <- function(line) {
  
  # Initialize variables for extracted information
  lhs <- "target"
  rhs <- ""
  
  # Regular expression patterns
  target_plus_equal_regex <- "^\\s*target\\s*\\+=\\s*(.+);"
  target_equal_regex <- "^\\s*target\\s*=\\s*target\\s*\\+\\s*(.+);"
  
  # Check for 'target += ...' pattern
  if (grepl(target_plus_equal_regex, line)) {
    rhs <- sub(target_plus_equal_regex, "\\1", line)
  }
  # Check for 'target = target + ...' pattern
  else if (grepl(target_equal_regex, line)) {
    rhs <- sub(target_equal_regex, "\\1", line)
  }
  
  # Further processing if RHS contains lpmf or lpdf
  parsed_info <- if (grepl("_lp(mf|df)\\s*\\(", rhs)) {
    parse_lpmf_lpdf(rhs)
  } else {
    list(rhs_expression = rhs)
  }
  parsed_info["lhs"] <- lhs
  
  return(parsed_info)
}

parse_stan_lines <- function(block_content) {
  
  # create a modifiable lines object
  lines <- block_content
  
  if(length(lines) == 0) return(lines)
  
  # Initialize a data frame to hold line information
  line_info <- data.frame(line = character(), type = character(), loop_depth = integer(), loop_start_line = integer(), stringsAsFactors = FALSE)
  
  # Initialize loop depth counter and stack
  loop_depth <- 0
  open_brackets_stack <- list()
  
  # Regular expressions for detecting types of lines
  declaration_regex <- "^\\s*(int|real|vector|row_vector|matrix|array)"
  lpmf_lpdf_regex <- "_lp(mf|df)\\s*\\("
  target_assignment_regex <- "^\\s*target\\s*(\\+=|=)"
  
  # Iterate over the lines
  for (i in 1:length(lines)) {
    line <- lines[i]
    
    # Initialize variables to store parsed information
    parsed_type <- ""
    parsed_info <- list()
    
    # Determine the type of line and parse accordingly
    if (grepl(declaration_regex, line)) {
      parsed_type <- "declaration"
      parsed_info <- parse_declaration(line)
    } else if (grepl(lpmf_lpdf_regex, line)) {
      parsed_type <- "lpmf_lpdf"
      parsed_info <- parse_lpmf_lpdf(line)
    } else if (grepl("log_sum_exp", line)) {
      parsed_type <- "log_sum_exp"
      parsed_info <- parse_lse(line)
    } else if (grepl(target_assignment_regex, line)) {
      parsed_type <- "sampling_target"
      parsed_info <- parse_target_assignment(line)
    } else if (grepl("=", line)) {
      parsed_type <- "assignment"
      parsed_info <- parse_assignment(line)
    } else if (grepl("~", line)) {
      parsed_type <- "sampling"
      parsed_info <- parse_sampling(line)
    } else if (grepl("^for\\s*\\(|^while\\s*\\(", line)) {
      parsed_type <- "loop"
      parsed_info <- parse_loop(line)
      loop_depth <- loop_depth + 1
      open_brackets_stack[[length(open_brackets_stack) + 1]] <- i
    } else if (grepl("^\\s*\\}", line)) {
      if (length(open_brackets_stack) > 0) {
        start_line <- open_brackets_stack[[length(open_brackets_stack)]]
        parsed_type <- "loop"
        parsed_info <- list(loop_action = "close", 
                            start_line = start_line,
                            loop_type = line_info$loop_type[start_line])
        open_brackets_stack <- open_brackets_stack[-length(open_brackets_stack)]
        loop_depth <- loop_depth - 1
      } else {
        parsed_type <- "operation"
      }
    } else {
      parsed_type <- "other"
    }
    
    # Add parsed information to the data frame
    line_info <- bind_rows(line_info, 
                           c(list(line = line, 
                                  type = parsed_type, 
                                  loop_depth = loop_depth - ifelse(parsed_type == "loop" && parsed_info$loop_action == "open", 1, 0)
                           ), 
                           parsed_info))
    
    line_info
  }
  
  if(any(line_info$type == "loop")){
    loop_info <- data.frame(open_ind = which(line_info$loop_action == "open"),
                            close_ind = which(line_info$loop_action == "close"),
                            loop_var = line_info$loop_variable[which(line_info$loop_action == "open")],
                            loop_range = line_info$loop_range[which(line_info$loop_action == "open")])
    
    for(i in 1:nrow(loop_info)){
      loop_inds <- (loop_info$open_ind[i]+1):(loop_info$close_ind[i]-1)
      if(all(is.na(line_info$loop_variable[loop_inds]))){
        line_info$loop_variable[loop_inds] <- loop_info$loop_var[i]
        line_info$loop_range[loop_inds] <- loop_info$loop_range[i]
      } else {
        line_info$loop_variable[loop_inds] <- paste0(line_info$loop_variable[loop_inds], ", ", loop_info$loop_var[i])
        line_info$loop_range[loop_inds] <- paste0(line_info$loop_range[loop_inds], ", ", loop_info$loop_range[i])
      }
    }
  }
  
  return(line_info)
}

#### generative functions ####

interpret_assignment <- function(line) {
  # Construct the LHS of the assignment
  lhs <- line$lhs_name
  if (nchar(line$lhs_index) > 0) {
    lhs <- paste0(lhs, "[", line$lhs_index, "]")
  }
  
  # Process the RHS of the assignment
  rhs <- interpret_operation(line$rhs_expression)
  
  # Combine LHS and RHS to form the assignment statement
  return(paste0(lhs, " <- ", rhs))
}


interpret_operation <- function(stan_expression) {
  
  #extract function calls if any exist in operation
  functions_used <- regmatches(stan_expression, 
                               gregexpr("\\b\\w+(?=\\()", 
                                        stan_expression, 
                                        perl = TRUE))[[1]]
  function_mapping <- list(rep_vector = "rep",
                           to_vector = "",
                           inv_logit = "plogis")
  
  #check if functions are in map list or in the current environment and warn user accordingly
  if(length(functions_used) >= 1){
    f_in_R <- sapply(functions_used, exists, mode = "function")
    f_in_map <- functions_used %in% names(function_mapping)
    f_not_anywhere <- !f_in_R & !f_in_map
    f_in_R_not_map <- f_in_R & !f_in_map
    f_in_R_not_map_envirs <- sapply(functions_used[f_in_R_not_map], function(fname)
      environmentName(environment(get(fname))))
    f_in_R_not_map_envirs[f_in_R_not_map_envirs == ""] <- "base"
    
    #warnings for using R funcs in place of Stan funcs
    if(any(f_in_R_not_map)){
      nfunc <- sum(f_in_R_not_map)
      if(nfunc == 1){
        function_string_concat <- paste0("<<", functions_used[f_in_R_not_map], "()>>")
        function_string_concat_R <- paste0("<<", f_in_R_not_map_envirs, "::", 
                                           functions_used[f_in_R_not_map], "()>>")
      } else if(nfunc == 2){
        function_string_concat <- paste0("<<", functions_used[f_in_R_not_map][1], "()>> and <<",
                                         functions_used[f_in_R_not_map][2], "()>>")
        function_string_concat_R <- paste0("<<", f_in_R_not_map_envirs[1], "::", 
                                           functions_used[f_in_R_not_map][1], "()>> and <<", 
                                           f_in_R_not_map_envirs[2], "::",
                                           functions_used[f_in_R_not_map][2], "()>>")
      } else {
        function_string_concat <- paste0(paste0(sapply(functions_used[f_in_R_not_map][-nfunc], 
                                                       function(fi){paste0("<<", fi, "()>>")}), collapse = ", "), 
                                         ", and <<", functions_used[f_in_R_not_map][nfunc], "()>>")
        function_string_concat_R <- paste0(paste0(sapply(functions_used[f_in_R_not_map][-nfunc], 
                                                         function(fi){
                                                           paste0("<<", f_in_R_not_map_envirs[fi], 
                                                                  "::", 
                                                                  fi, 
                                                                  "()>>")}), collapse = ", "), 
                                           ", and <<", f_in_R_not_map_envirs[nfunc], 
                                           "::", functions_used[f_in_R_not_map][nfunc], "()>>")
      }
      warning(paste0("Stan function", ifelse(nfunc>1, "s ", " "), function_string_concat, 
                     " not found in map but found in R_GlobalEnv. Using R-function", 
                     ifelse(nfunc>1, "s ", " "), 
                     function_string_concat_R, " in their place."))
    }
    
    #warnings for not having any matching functions at all
    if(any(f_not_anywhere)){
      
      nfunc <- sum(f_not_anywhere)
      if(nfunc == 1){
        function_string_concat <- paste0("<<", functions_used[f_not_anywhere], "()>>")
      } else if(nfunc == 2){
        function_string_concat <- paste0("<<", functions_used[f_not_anywhere][1], "()>> and <<",
                                         functions_used[f_not_anywhere][2], "()>>")
      } else {
        function_string_concat <- paste0(paste0(sapply(functions_used[f_not_anywhere][-nfunc], 
                                                       function(fi){paste0("<<", fi, "()>>")}), collapse = ", "), 
                                         ", and <<", functions_used[f_not_anywhere][nfunc], "()>>")
      }
      warning(paste0("Stan function", ifelse(nfunc>1, " ", "s "), function_string_concat, 
                     " not found in map or in R_GlobalEnv. Implement before running code."))
    }
    
    # Replace specific Stan functions with equivalent R functions or code
    if(any(f_in_map)){
      for(fi in which(f_in_map)){
        stan_expression <- gsub(paste0(functions_used[fi], "\\("), 
                                paste0(function_mapping[[functions_used[fi]]], "\\("), stan_expression)
      }
    }
  }
  
  # Handle element-wise operations (*, /, +, -) - In R, these are by default element-wise
  # Stan's element-wise operations (e.g., .*) are the same as R's default, so we can remove the '.'
  stan_expression <- gsub("\\.\\*", "*", stan_expression)
  stan_expression <- gsub("\\./", "/", stan_expression)
  stan_expression <- gsub("\\.\\+", "+", stan_expression)
  stan_expression <- gsub("\\.-", "-", stan_expression)
  
  # Return the translated stan_expression
  return(stan_expression)
}


interpret_loop <- function(line) {
  
  if (line$loop_action == "open") {
    # Construct the loop header
    loop_variable <- line$loop_variable
    loop_range <- line$loop_range
    if(grepl(",", loop_variable)){
      loop_variable <- trimws(strsplit(loop_variable, ",")[[1]])[1]
      loop_range <- trimws(strsplit(loop_range, ",")[[1]])[1]
    }
    
    loop_header <- paste0("for (", loop_variable, " in ", loop_range, ") {")
    return(loop_header)
  } else if (line$loop_action == "close") {
    # Close the loop
    return("}")
  } else {
    # Handle unexpected cases
    return("# Unexpected loop structure")
  }
}


interpret_sampling <- function(line, dat, samps, sample_index = 1, post_pred_sim = TRUE, sim = TRUE, bound_declarations = NA) {
  
  # LHS of the sampling statement
  lhs <- line$var_name
  if(is.na(line$var_name)){
    lhs <- line$outcome_variable
  }
  lhs_name <- gsub("\\[.*\\]", "", lhs)
  lhs_index <- gsub(lhs_name, "", lhs)
  
  # Check if the LHS is in the data block and handle accordingly
  if (lhs_name %in% names(dat)) {
    if (sim) {
      # Sample from the specified distribution (prior predictive simulation)
      return(paste0(lhs, " <- ", generate_sampling_code(line)))
    } else {
      # Find the density or mass for the parameter in the specified distribution
      return(paste0(lhs, " <- ", generate_density_code(line)))
    }
    
  } else {
    if (post_pred_sim) {
      
      # Posterior predictive simulation
      
      #retrieve sampled param value
      sampled_param <- munge_samps(lhs_name, subset_samps(lhs_name, samps[sample_index, , drop = FALSE]))
      
      # generate appropriate R code
      if (is.vector(sampled_param)) {
        # Vector: Use c() to create a vector in R code
        r_code <- paste0(lhs, " <- c(", paste(sampled_param, collapse = ", "), ")", lhs_index)
      } else if (is.matrix(sampled_param)) {
        # Matrix: Use matrix() to create a matrix in R code
        matrix_values <- paste(sampled_param, collapse = ", ")
        if(nchar(matrix_values) > 1000){
          param_breaks <- unique(c(seq(1, length(sampled_param), by = 200), length(sampled_param) + 1))
          matrix_values_chunked <- sapply(2:length(param_breaks), function(i){
            param_breaks[i-1]:(param_breaks[i]-1)
            paste(sampled_param[param_breaks[i-1]:(param_breaks[i]-1)], collapse = ", ")
          })
          matrix_values <- paste0(c(paste0(matrix_values_chunked[1:(length(matrix_values_chunked)-1)], ","), 
                                    matrix_values_chunked[length(matrix_values_chunked)]), collapse = "\n")
        }
        r_code <- paste0(lhs, " <- matrix(c(", matrix_values, "), nrow = ", nrow(sampled_param), ", ncol = ", ncol(sampled_param), ")", lhs_index)
      } else if (is.array(sampled_param)) {
        # Array: More complex, handle as needed
        r_code <- "# Array assignment not yet implemented"
      } else {
        # Scalar or other types
        r_code <- paste0(lhs, " <- ", sampled_param, lhs_index)
      }
      return(r_code)
      
    } else {
      # Prior predictive simulation
      return(paste0(lhs, " <- ", generate_sampling_code(line))) 
    }
  }
}

distribution_maps <- function(type = c("r", "q", "p", "d")[1]){
  if(type == "p"){
    return(list(
      std_normal = "pnorm(bound, mean = 0, sd = 1)",
      normal = "pnorm(bound, mean = param1, sd = param2)",  # Assuming 'mean, sd' parameterization
      beta = "pbeta(bound, shape1 = param1, shape2 = param2)",
      binomial = "pbinom(bound, size = param1, prob = param2)",
      gamma = "pgamma(bound, shape = param1, rate = param2)",  # Assuming 'shape, rate' parameterization
      student_t = "pt(bound, df = param1, ncp = param2)",
      double_exponential = "plaplace(bound, location = param1, scale = param2)",  # Note: plaplace may not be directly available in base R
      exponential = "pexp(bound, rate = param1)",
      lognormal = "plnorm(bound, meanlog = param1, sdlog = param2)",
      chi_square = "pchisq(bound, df = param1)",
      uniform = "punif(bound, min = param1, max = param2)",
      poisson = "ppois(bound, lambda = param1)",
      beta_binomial = "extraDistr::pbetabinom(bound, size = param1, prob = param2, prob2 = param3)"  # Note: pbetabinom may not be directly available as named; placeholder for conceptual use
      # Add other distributions as needed
    ))
  } else if(type == "q"){
    return(list(
      std_normal = "qnorm(target_quantile, mean = 0, sd = 1)",
      normal = "qnorm(target_quantile, mean = param1, sd = param2)",
      beta = "qbeta(target_quantile, shape1 = param1, shape2 = param2)",
      binomial = "qbinom(target_quantile, size = param1, prob = param2)",
      gamma = "qgamma(target_quantile, shape = param1, rate = param2)",
      student_t = "qt(target_quantile, df = param1, ncp = param2)",
      exponential = "qexp(target_quantile, rate = param1)",
      lognormal = "qlnorm(target_quantile, meanlog = param1, sdlog = param2)",
      chi_square = "qchisq(target_quantile, df = param1)",
      uniform = "qunif(target_quantile, min = param1, max = param2)",
      poisson = "qpois(target_quantile, lambda = param1)",
      double_exponential = "qaplace(target_quantile, location = param1, scale = param2)",
      beta_binomial = "extraDistr::qbbinom(target_quantile, size = param1, alpha = param2, beta = param3)"
      # Add other distributions as needed
    ))
  } else if(type == "r"){
    return(list(
      std_normal = "rnorm(n = lhs_dim, mean = 0, sd = 1)",
      normal = "rnorm(n = lhs_dim, mean = param1, sd = param2)",  # Assuming 'mean, sd' parameterization
      beta = "rbeta(n = lhs_dim, shape1 = param1, shape2 = param2)",
      binomial = "rbinom(n = lhs_dim, size = param1, prob = param2)",
      gamma = "rgamma(n = lhs_dim, shape = param1, rate = param2)",  # Assuming 'shape, rate' parameterization
      student_t = "rt(n = lhs_dim, df = param1, ncp = param2)",
      double_exponential = "rlaplace(n = lhs_dim, location = param1, scale = param2)",
      exponential = "rexp(n = lhs_dim, rate = param1)",
      lognormal = "rlnorm(n = lhs_dim, meanlog = param1, sdlog = param2)",
      chi_square = "rchisq(n = lhs_dim, df = param1)",
      uniform = "runif(n = lhs_dim, min = param1, max = param2)",
      poisson = "rpois(n = lhs_dim, lambda = param1)",
      beta_binomial = "extraDistr::rbbinom(n = lhs_dim, size = param1, alpha = param2, beta = param3)"
      # Add other distributions as needed
    ))
  } else if(type == "d"){
    return(list(
      normal = "dnorm(x = dat$VAR_NAME, mean = param1, sd = param2)",
      beta = "dbeta(x = dat$VAR_NAME, shape1 = param1, shape2 = param2)",
      binomial = "dbinom(x = dat$VAR_NAME, size = param1, prob = param2)",
      gamma = "dgamma(x = dat$VAR_NAME, shape = param1, rate = param2)",
      student_t = "dt(x = dat$VAR_NAME, df = param1, ncp = param2)",
      double_exponential = "dexp(x = dat$VAR_NAME, rate = param1)",  # Assuming 'rate' parameterization
      exponential = "dexp(x = dat$VAR_NAME, rate = param1)",
      lognormal = "dlnorm(x = dat$VAR_NAME, meanlog = param1, sdlog = param2)",
      chi_square = "dchisq(x = dat$VAR_NAME, df = param1)",
      uniform = "dunif(x = dat$VAR_NAME, min = param1, max = param2)",
      poisson = "dpois(x = dat$VAR_NAME, lambda = param1)",
      beta_binomial = "extraDistr::dbbinom(x = dat$VAR_NAME, size = param1, alpha = param2, beta = param3)"
      # Add other distributions as needed
    ))
  }
}

generate_sampling_code <- function(line) {
  distribution <- line$distribution
  parameters <- line$parameters
  param_names <- paste0("param", seq_along(strsplit(parameters, ",")[[1]]))
  dimensions <- line$dimensions
  if(is.na(dimensions) | dimensions == ""){
    dimensions <- 1
  } else if(!(is.na(line$lhs_index) | line$lhs_index == "")){
    dimensions <- paste0("prod(dim(", line$var_name, ")")
  }
  
  #check if bounds exist
  bounded <- F
  if((line$lower_bound != "-Inf" || line$upper_bound != "Inf") & 
     distribution %in% c("normal", "beta", "gamma", "exponential", "lognormal", "chi_square", "uniform", "poisson", "std_normal")){
    bounded <- T
    lb <- line$lower_bound
    ub <- line$upper_bound
  }
  
  if(bounded){
    # Mapping Stan distributions to R functions for evaluating upper and lower bounds
    p_distribution_map <- distribution_maps("p")
    if(!(distribution %in% names(p_distribution_map))){
      warning(paste0("Distribution named <<", distribution, 
                     ">> not found in conversion table, requiring implementation."))
    }
    
    #evaluate quantiles of bounds
    ub_code <- p_distribution_map[[distribution]]
    lb_code <- p_distribution_map[[distribution]]
    if(parameters != ""){
      for (i in seq_along(param_names)) {
        ub_code <- gsub(param_names[i], strsplit(parameters, ",")[[1]][i], ub_code)
        lb_code <- gsub(param_names[i], strsplit(parameters, ",")[[1]][i], lb_code)
      }  
    }
    ub_code <- gsub("bound", ub, ub_code)
    lb_code <- gsub("bound", lb, lb_code)
    bounds_diff_code <- paste0("(", ub_code, " - ", lb_code, ")")
    
    #execute inverse transform sampling -- first find the target quantile
    uniform_sample_code <- paste0("runif(n = ", dimensions, ", min = 0, max = 1)")
    target_quantile_code <- paste0(lb_code, " + ", bounds_diff_code, " * ", uniform_sample_code)
    
    #now find where this quantile falls in the target distribution
    q_distribution_map <- distribution_maps("q")
    if(!(distribution %in% names(q_distribution_map))){
      warning(paste0("Distribution named <<", distribution, 
                     ">> not found in conversion table, requiring implementation."))
    }
    
    inverse_transform_code <- q_distribution_map[[distribution]]
    if(parameters != ""){
      for (i in seq_along(param_names)) {
        inverse_transform_code <- gsub(param_names[i], strsplit(parameters, ",")[[1]][i], inverse_transform_code)
        inverse_transform_code <- gsub(param_names[i], strsplit(parameters, ",")[[1]][i], inverse_transform_code)
      }
    }
    inverse_transform_code <- gsub("target_quantile", target_quantile_code, inverse_transform_code)
    
    
    # Combine the codes to sample within bounds and transform back
    r_function <- inverse_transform_code
    return(r_function)
    
  } else {
    # Mapping Stan distributions to R functions for random number gen
    r_distribution_map <- distribution_maps("r")
    if(!(distribution %in% names(r_distribution_map))){
      warning(paste0("Distribution named <<", distribution, 
                     ">> not found in conversion table, requiring implementation."))
    }
    
    # Generate R code for sampling
    r_function <- r_distribution_map[[distribution]]
    if (is.null(r_function)) {
      return(paste("# No R equivalent for", r_distribution_map, "distribution"))
    } else {
      # Replace placeholders with actual parameters
      for (i in seq_along(param_names)) {
        r_function <- gsub(param_names[i], strsplit(parameters, ",")[[1]][i], r_function)
      }
      r_function <- gsub("lhs_dim", dimensions, r_function)
      return(r_function)
    }
  }
  
}


generate_density_code <- function(line, dat) {
  distribution <- line$distribution
  parameters <- line$parameters
  var_name <- line$var_name
  bounded <- F
  if((line$lower_bound != "-Inf" || line$upper_bound != "Inf") & 
     distribution %in% c("normal", "beta", "gamma", "exponential", "lognormal", "chi_square", "uniform", "poisson", "std_normal")){
    bounded <- T
    lb <- line$lower_bound
    ub <- line$upper_bound
  }
  
  # Mapping Stan distributions to R density/mass functions
  d_distribution_map <- distribution_maps("d")
  if(!(distribution %in% names(d_distribution_map))){
    warning(paste0("Distribution named <<", distribution, 
                   ">> not found in conversion table, requiring implementation."))
  }
  
  # Generate R code for density/mass calculation
  r_function <- d_distribution_map[[distribution]]
  if (is.null(r_function)) {
    return(paste("# No R equivalent for", distribution, "distribution"))
  } else {
    # Replace placeholder parameters with actual parameters
    param_names <- paste0("param", seq_along(strsplit(parameters, ",")[[1]]))
    for (i in seq_along(param_names)) {
      r_function <- gsub(param_names[i], strsplit(parameters, ",")[[1]][i], r_function)
    }
    
    # Replace placeholder random variable name with actual random variable name
    r_function <- gsub("VAR_NAME", var_name, r_function)
    
    #transform density or mass to respect bounds (ie, multiply by the total slice of probability in the distribution)
    if(bounded){
      p_distribution_map <- distribution_maps("p")
      if(!(distribution %in% names(p_distribution_map))){
        warning(paste0("Distribution named <<", distribution, 
                       ">> not found in conversion table, requiring implementation."))
      }
      
      #evaluate quantiles of bounds
      ub_code <- p_distribution_map[[distribution]]
      lb_code <- p_distribution_map[[distribution]]
      if(parameters != ""){
        for (i in seq_along(param_names)) {
          ub_code <- gsub(param_names[i], strsplit(parameters, ",")[[1]][i], ub_code)
          lb_code <- gsub(param_names[i], strsplit(parameters, ",")[[1]][i], lb_code)
        }  
      }
      ub_code <- gsub("bound", ub, ub_code)
      lb_code <- gsub("bound", lb, lb_code)
      bounds_diff_code <- paste0("(", ub_code, " - ", lb_code, ")")
      
      r_function <- paste0("(", r_function, ") / ", bounds_diff_code)
    }
    
    return(r_function)
  }
}



interpret_declaration <- function(line, dat, samps) {
  if (line$block_name == "data") {
    
    # Read from input data
    var_name <- line$var_name
    r_code <- paste0(var_name, " <- dat$", var_name)
    
  } else {
    
    # Initialize the variable in R
    var_name <- line$var_name
    declaration_type <- line$declaration_type
    dimensions <- line$dimensions
    # Default case for unknown types
    r_code <- paste0("# Unknown type for ", var_name)
    
    # Handle different declaration types
    if (declaration_type == "int" || declaration_type == "real") {
      if (dimensions == "1") {
        # Scalar
        r_code <- paste0(var_name, " <- NA")
      } else {
        # Vector
        if(grepl(",", dimensions)){
          declaration_type <- "array"
        } else {
          declaration_type <- "vector"
        }
      }
    }
    
    if (declaration_type == "vector") {
      # Vector
      r_code <- paste0(var_name, " <- rep(NA, ", dimensions, ")")
    } else if (declaration_type == "matrix") {
      # Matrix
      matrix_dims <- strsplit(dimensions, ",")[[1]]
      r_code <- paste0(var_name, " <- matrix(NA, nrow = ", matrix_dims[1], ", ncol = ", matrix_dims[2], ")")
    } else if (declaration_type == "array") {
      # Array with more than 2 dimensions
      array_dims <- paste("dim = c(", dimensions, ")", sep = "")
      r_code <- paste0(var_name, " <- array(NA, ", array_dims, ")")
    }
  }
  
  return(r_code)
}

interpret_lpmf_lpdf <- function(line, dat, samps, sample_index = 1, post_pred_sim = TRUE, sim = TRUE, bound_declarations = NA) {
  sampling_line <- line
  sampling_line$var_name <- line$outcome_variable
  
  
  if (!is.na(line$mixture_probability)) {
    
    # Extract indices from lhs and outcome_variable
    lhs_var <- gsub("\\[.*\\]", "", line$lhs)
    lhs_index <- gsub(".*\\[|\\].*", "", line$lhs)
    outcome_var_index <- gsub(".*\\[|\\].*", "", line$outcome_variable)
    outcome_base <- gsub("\\[.*\\]", "", line$outcome_variable)
    
    if(lhs_index == 1){
      #change lp class to numeric_mixture
      class_def_code <- "if(!isClass('numeric_mixture')){\n\tsetClass(\n\t\t\"numeric_mixture\",\n\t\tcontains = \"numeric\",\n\t\tslots = c(probabilities = \"numeric\", outcome_var = \"character\")\n\t)\n}\n"
      
      # class_set_code <- paste0("class(", lhs_var, ") <- \"numeric_mixture\"\n") #hm this is giving me errors when accessing attributes
      class_set_code <- paste0(lhs_var, " <- new(\"numeric_mixture\", .Data = ", lhs_var, ", probabilities = numeric(length(", lhs_var, ")), outcome_var = \"", outcome_base, "\")\n\n")
      
      # outcome_var_code <- paste0(lhs_var, "@outcome_var <- \"", outcome_base, "\"\n\n")
      outcome_var_code <- ""
      class_code <- paste0(class_def_code, class_set_code, outcome_var_code)
    } else {
      class_code <- ""
    }
    
    # Generate R code for mixture probability and component
    mixture_prob_code <- interpret_operation(line$mixture_probability)
    component_code <- interpret_sampling(line = sampling_line, dat, samps, sample_index, post_pred_sim = F, sim = sim, bound_declarations)
    
    
    if(!is.na(line$loop_variable)){
      
      numeric_code <- paste0(lhs_var, "[", lhs_index, "] <- ", 
                             gsub(paste0(gsub("\\]", "\\\\]", gsub("\\[", "\\\\[", line$outcome_variable)), 
                                         " <- "), "", component_code), "\n")
      prob_code <- paste0(lhs_var, "@probabilities[", lhs_index, "] <- ", mixture_prob_code, "\n")
      r_code <- paste0(class_code, numeric_code, prob_code)
      
      return(r_code)
    }
    
  }
  
  r_code <- interpret_sampling(sampling_line, dat, samps, sample_index, post_pred_sim, sim, bound_declarations)
  
  return(r_code)
}


interpret_lse <- function(line, dat, samps, sim = TRUE) {
  
  if(!is.na(line$loop_variable)){
    
    if (sim) {
      # Sample from the mixture components
      
      #this records the mixture component index sampled
      r_code_1 <- paste0("if(is.null(", "out[[paste0(", line$lse_argument, "@outcome_var, \"_mixcomp\")]])){\n\tout[[paste0(", 
                       line$lse_argument, "@outcome_var, \"_mixcomp\")]] <- numeric(length(", line$lse_argument, "@outcome_var))\n}\n")
      r_code_2 <- paste0("out[[paste0(", line$lse_argument, "@outcome_var, \"_mixcomp\")]][", line$loop_variable,"] <- ", 
                       "sample(1:", "length(", line$lse_argument, ")", ", 1, prob = ", line$lse_argument, "@probabilities)\n")
      r_code_3 <- paste0("out[[", line$lse_argument, "@outcome_var]][", line$loop_variable,"] <- ", line$lse_argument, 
                       "[", "out[[paste0(", line$lse_argument, "@outcome_var, \"_mixcomp\")]][", line$loop_variable,"]","]\n")
      r_code <- paste0(r_code_1, r_code_2, r_code_3)
      
      # r_code <- paste0("out[[", line$lse_argument, "@outcome_var]][", line$loop_variable,"] <- ", line$lse_argument,
      #                  "[", "sample(1:", "length(", line$lse_argument, ")", ", 1, prob = ", line$lse_argument, "@probabilities)]\n") #this does not record component
      
    } else {
      # Compute weighted average of mixture components
      r_code <- paste0("out[[", line$lse_argument, "@outcome_var]][", line$loop_variable,"] <- sum(", line$lse_argument, " * ", line$lse_argument, "@probabilities)\n")
      
    }
    
    return(r_code)
  }
  
  return("#error\n\n")
}



interpret_lines <- function(line, dat, samps, sample_index, post_pred_sim, sim) {
  switch(line$type,
         "declaration" = interpret_declaration(line, dat, samps),
         "assignment" = interpret_assignment(line),
         # "operation" = interpret_operation(line, dat, samps),
         "loop" = interpret_loop(line),
         "sampling" = interpret_sampling(line, dat, samps, sample_index, post_pred_sim, sim),
         "log_sum_exp" = interpret_lse(line, dat, samps, sim = sim),
         "lpmf_lpdf" = interpret_lpmf_lpdf(line, dat, samps, sample_index, post_pred_sim, sim),
         "Unknown type")
}




# code for implementing and troubleshooting line types

# #declarations
# dcls <- parsed_lines[parsed_lines$type == "declaration",]
# sapply(1:nrow(dcls), function(i) cbind(dcls$line[i], dcls$block_name[i], interpret_declaration(dcls[i,], dat, samps)))
# 
# #assignments
# asms <- parsed_lines[parsed_lines$type == "assignment",]
# sapply(1:nrow(asms), function(i) cbind(asms$line[i], asms$block_name[i], interpret_assignment(asms[i,])))
# 
# #loops
# loops <- parsed_lines[parsed_lines$type == "loop",]
# sapply(1:nrow(loops), function(i) cbind(loops$line[i], loops$block_name[i], interpret_loop(loops[i,])))
# 
# #sampling statements
# sss <- parsed_lines[parsed_lines$type == "sampling",]
# head(sss[,!apply(apply(sss,2,is.na), 2, all)])
# sapply(1:nrow(sss), function(i){print(i); (#cbind(sss$line[i], sss$block_name[i], 
#   interpret_sampling(line = sss[i,], dat, samps, 
#                      sample_index = 1, 
#                      post_pred_sim = TRUE, sim = TRUE, bound_declarations = NA))
# })
# 
# interpret_sampling(line = sss[i,], dat, samps, 
#                    sample_index = 1, 
#                    post_pred_sim = F, sim = TRUE, bound_declarations = NA)
# 
#log prob statements
# lps <- parsed_lines[parsed_lines$type == "lpmf_lpdf",]
# head(lps[,!apply(apply(lps,2,is.na), 2, all)])
# sapply(1:nrow(lps), function(i){print(i); (#cbind(lps$line[i], lps$block_name[i],
#   interpret_lpmf_lpdf(line = lps[i,], dat, samps,
#                       sample_index = 1,
#                       post_pred_sim = F, sim = TRUE, bound_declarations = NA))
# })
# i = 2
# line <- lps[i,]
# line
# interpret_lpmf_lpdf(line = line, dat, samps,
#                     sample_index = 1,
#                     post_pred_sim = F, sim = F, bound_declarations = NA)


# #log-sum-exp statements
# lses <- parsed_lines[parsed_lines$type == "log_sum_exp",]
# line <- lses[1,]
# 
# cat(interpret_lse(line))


#### composite functions ####

retrieve_block_objects <- function(stan_code, block = c("data", 
                                                             "transformed data", 
                                                             "parameters", 
                                                             "transformed parameters", 
                                                             "model")[1],
                                      type = "declaration"){
  block_contents <- parse_stan_blocks(clean_stan(stan_code))
  parsed_data_lines <- lapply(block_contents[block], parse_stan_lines)
  parsed_data_lines[[block]][parsed_data_lines[[block]]$type == type, 
                             c("var_name", "declaration_type", "bounds")]
}

loop_over_param <- function(var_name, dims, set_as = 1) {
  # Splitting the dimension string into components and trimming any whitespace
  dim_components <- trimws(strsplit(dims, ",")[[1]])
  
  # Initializing variables to build the for-loop string
  loop_start <- ""
  loop_end <- ""
  indent_level <- ""
  
  # Building the loop structure dynamically based on the number of dimensions
  for (i in seq_along(dim_components)) {
    # Creating the variable index name
    index_var <- paste0(var_name, "_i_", i)
    # Adding to the start of the loop string (opening for-loops)
    loop_start <- paste0(loop_start, indent_level, "for(", index_var, " in 1:", dim_components[i], "){\n")
    # Preparing the end of the loop string (closing braces)
    loop_end <- paste0(indent_level, "}\n", loop_end)
    # Increasing the indent for the next level of nesting
    indent_level <- paste0(indent_level, " ")
  }
  
  # Constructing the assignment line with the appropriate indentation
  index_vars <- paste0(var_name, "_i_", seq_along(dim_components), collapse = ", ")
  assignment_line <- paste0(indent_level, var_name, "[", index_vars, "] = ", set_as, ";\n")
  
  # Combining all parts to form the complete loop structure
  loop_code <- paste0(loop_start, assignment_line, loop_end)
  
  return(loop_code)
}


flatten_model <- function(stan_code, scale_identifier = "sd_", set_as = 1, put_params_in_data = F){
  
  #this function parses a Stan model, identifies scale parameters,
  #and sets them all to 1 (or whatever is in <set_as>), instead of sampling
  #alternatively, if put_params_in_data = T, it includes those parameters in the data block
  
  #some basic pre-processing
  clean_stan_info <- clean_stan(stan_code, return_annot = T, R_style = F)
  stan_annot <- clean_stan_info$annot
  
  #retrieve blocks and modify annotation
  block_contents_info <- parse_stan_blocks(clean_stan_info$code, return_annot_mod = T)
  block_contents <- block_contents_info$block_contents
  stan_annot <- stan_annot[-block_contents_info$block_mod,]
  block_contents <- block_contents[sapply(block_contents, length) > 0]
  
  #parse block content and retrieve all declared params
  parsed_lines <- lapply(block_contents, parse_stan_lines)
  all_params <- data.frame(name = parsed_lines$parameters$var_name, 
                           dim = parsed_lines$parameters$dimensions)
  
  #combine parsed_lines with annotation 
  #(should have done this in parse_Stan whoops, will refactor later)
  nlines_per_block <- unlist(lapply(parsed_lines, nrow))
  block_starts <- cumsum(c(1, nlines_per_block[-length(nlines_per_block)]))
  names(block_starts) <- names(nlines_per_block)
  parsed_lines <- lapply(setNames(seq_along(parsed_lines), names(parsed_lines)), function(i){
    cbind(parsed_lines[[i]], stan_annot[block_starts[i]:(block_starts[i] + nlines_per_block[i] - 1),])
  })
  
  #ID lines to modify
  scale_declaration_inds <- grepl(pattern = scale_identifier, all_params$name) &
    parsed_lines$parameters$type == "declaration"  
  scale_params <- all_params[scale_declaration_inds,]
  scale_sampling_inds <- parsed_lines$model$type == "sampling" & 
    parsed_lines$model$var_name %in% scale_params$name
  
  #create modified lines
  if(put_params_in_data){
    new_lines <- parsed_lines
    new_lines$data <- rbind(new_lines$data, new_lines$parameters[scale_declaration_inds,])
    new_lines$parameters <- new_lines$parameters[!scale_declaration_inds,]
    new_lines$model <- new_lines$model[!scale_sampling_inds,]
    
  } else {
    if(length(set_as) == 1){
      new_declarations <- character(length(scale_params$name))
      new_declarations[scale_params$dim == "1"] <- paste0(scale_params$name, " = ", set_as, ";")
      if(any(scale_params$dim != "1")){
        new_declarations[scale_params$dim != "1"] <- sapply(which(scale_params$dim != "1"), function(flpi){
          loop_over_param(var_name = scale_params$name[flpi], dims = scale_params$dim[flpi], set_as = set_as)
        })  
      }
      constrain_free_declarations <- gsub("<[^>]*>", "", parsed_lines$parameters$line[scale_declaration_inds])
      new_declarations <- paste0(constrain_free_declarations, "\n", 
                                 parsed_lines$parameters$leading_ws[scale_declaration_inds],
                                 new_declarations)  
    } else {
      new_declarations <- setNames(character(length(scale_params$name)), scale_params$name)
      new_declarations[scale_params$dim == "1"] <- paste0(scale_params$name, " = ", set_as[scale_params$name], ";")
      if(any(scale_params$dim != "1")){
        new_declarations[scale_params$dim != "1"] <- sapply(which(scale_params$dim != "1"), function(flpi){
          loop_over_param(var_name = scale_params$name[flpi], dims = scale_params$dim[flpi], set_as = set_as[scale_params$name[flpi]])
        })  
      }
      constrain_free_declarations <- gsub("<[^>]*>", "", parsed_lines$parameters$line[scale_declaration_inds])
      
      new_declarations <- paste0(constrain_free_declarations, "\n", 
                                 parsed_lines$parameters$leading_ws[scale_declaration_inds],
                                 new_declarations)
    }
    
    #integrate into the original program
    new_lines <- parsed_lines
    new_lines$parameters$line[scale_declaration_inds] <- paste0("// ", 
              new_lines$parameters$line[scale_declaration_inds])
    new_lines$model$line[scale_sampling_inds] <- new_declarations
  }
  
  #create new lines inside blocks
  for(block_name in names(new_lines)){
    new_lines[[block_name]]$new_line <- paste0(new_lines[[block_name]]$preceding_content, 
                                                 new_lines[[block_name]]$leading_ws, 
                                                 new_lines[[block_name]]$line,
                                                 new_lines[[block_name]]$trailing_ws, 
                                                 new_lines[[block_name]]$comments)
  }
  
  #construct new block text
  new_stan_code <- lapply(names(new_lines), function(block_name){
    paste0(block_name, " {\n", 
           paste0(new_lines[[block_name]]$new_line, collapse = "\n"), 
           "\n}")
  })
  new_stan_code <- paste0(unlist(new_stan_code), collapse = "\n\n")
  new_stan_code
}

parse_Stan <- function(stan_code, dat = NA, samps = NA, output_file = NA, sample_index = 1, post_pred_sim = TRUE, sim = TRUE){
  
  #preprocess the code for easier interpretation
  clean_stan_info <- clean_stan(stan_code, T)
  clean_stan_annot <- clean_stan_info$annot
  clean_stan_code <- clean_stan_info$code
  
  #extract block information
  block_contents_info <- parse_stan_blocks(clean_stan_code, T)
  block_contents <- block_contents_info$block_contents
  empty_blocks <- sapply(block_contents, length) == 0
  
  #modify annot for block removal
  clean_stan_annot <- clean_stan_annot[-block_contents_info$block_mod,]
  
  #remove empty blocks if there are any
  if(any(empty_blocks)){
    block_contents <- block_contents[!empty_blocks]  
  }
  
  #parse lines in each block and concatenate
  parsed_lines <- lapply(block_contents, parse_stan_lines)
  parsed_lines <- cbind(bind_rows(parsed_lines), 
                        block_name = rep(names(parsed_lines), sapply(parsed_lines, nrow)))
  
  #get bounds in separately
  ul_bounds <- trimws(do.call(rbind, strsplit(parsed_lines$bounds, ",")))
  parsed_lines$lower_bound <- gsub(".*lower\\s*=\\s*([^, ]+).*", "\\1", parsed_lines$bounds)
  parsed_lines$upper_bound <- gsub(".*upper\\s*=\\s*([^, ]+).*", "\\1", parsed_lines$bounds)
  
  #if sampling line lhs does not have an index, make sure you're using the appropriate size 
  #(ie, assigning the whole param), as well as relaying bound info from the declaration
  all_decl <- parsed_lines[parsed_lines$type == "declaration", c("var_name", "dimensions", "bounds",
                                                                 "lower_bound", "upper_bound")]
  sampling_inds <- which(parsed_lines$type == "sampling")
  sampling_inds_no_index <- sampling_inds[is.na(parsed_lines$lhs_index[sampling_inds]) | 
                                            parsed_lines$lhs_index[sampling_inds] == ""]
  if(length(sampling_inds_no_index) > 0){
    parsed_lines[sampling_inds_no_index,
                 c("dimensions", "bounds", "lower_bound", "upper_bound")] <-
      all_decl[match(parsed_lines[sampling_inds_no_index,]$var_name, all_decl$var_name),
               c("dimensions", "bounds", "lower_bound", "upper_bound")]
  }
  
  processed_code <- character(nrow(parsed_lines))
  for(i in 1:nrow(parsed_lines)){
    processed_code[i] <- interpret_lines(line = parsed_lines[i,], dat = dat, samps = samps, 
                                         sample_index = sample_index, post_pred_sim = post_pred_sim, sim = sim)
  }
  
  #if there are transformed parameters, need to sample them before transforming parameters
  if("transformed parameters" %in% names(block_contents)){
    
    model_param_inds <- which(parsed_lines$block_name == "parameters" & parsed_lines$type == "declaration")
    param_names <- parsed_lines$var_name[model_param_inds]
    sampling_inds <- which(parsed_lines$block_name == "model" &
      parsed_lines$type == "sampling")
    sampling_inds <- sampling_inds[parsed_lines$var_name[sampling_inds] %in% param_names]
    insert_after_line <- max(which(parsed_lines$block_name == "parameters"))
    
    #construct new indices and reorder all relevant data objects
    reorder_vec <- append(x = (1:length(processed_code))[-sampling_inds], sampling_inds, insert_after_line)
    reordered_processed_code <- processed_code[reorder_vec]
    reordered_parsed_lines <- parsed_lines
    reordered_parsed_lines$block_name[sampling_inds] <- "parameters"
    reordered_parsed_lines <- reordered_parsed_lines[reorder_vec,]
    reordered_annot <- clean_stan_annot[reorder_vec,]
    
    #add in new line after insertion
    reordered_processed_code[insert_after_line] <- 
      paste0(reordered_processed_code[insert_after_line], "\n\n")
    
    #assign to standard data objects
    processed_code <- reordered_processed_code  
    parsed_lines <- reordered_parsed_lines
    clean_stan_annot <- reordered_annot
  }
  
  #add in tabs for whitespace... need to also get newline characters handled
  processed_code <- sapply(1:length(parsed_lines$loop_depth), function(i){
    tabs <- paste0(rep(x = "\t", times = parsed_lines$loop_depth[i]), collapse = "")
    paste0(tabs, gsub(pattern = "\n", replacement = paste0("\n", tabs), x = processed_code[i]), collapse = "")
  })

  #add in list object "out" to record dynamic variables
  newline_processed_code <- strsplit(processed_code, "\n")
  newline_processed_code <- sapply(newline_processed_code, function(lines){
    new_lines <- sapply(lines, function(line){
      if(grepl("<-", line) & !grepl("out\\[\\[", line)){
        lb <- strsplit(line, " <- ")[[1]]
        if(length(lb) == 2){
          lhs <- lb[1]
          rhs <- lb[2]
          
          if(grepl("class\\(", lhs)){
            class_prefix <- gsub("class\\(.*", "", lhs)
            lhs <- sub("^.*class\\(([^)]+)\\).*", "\\1", lhs)
          }
          
          var_name <- sub("([^\\[$@]*)[\\[$@].*", "\\1", lhs)
          var_subsetting <- gsub(var_name, "", lhs)
          out_var_name <- paste0("out[[\"", var_name, "\"]]")
          out_var_name <- gsub('\t|\n|\r', "", out_var_name)
          out_var_name <- paste0(out_var_name, var_subsetting)
          
          if(grepl("class\\(", lb[1])){
            out_var_name <- paste0("class(", out_var_name, ")")
            lhs <- paste0(class_prefix, "class(", lhs, ")")
          }
          
          new_line <- paste0(lhs, " <- ", out_var_name, " <- ", rhs)
          return(new_line)
          
        } else {
          print(line)
          stop("need two components to linebreak")
        }
      } else {
        return(line)
      }
    })
    return(paste0(new_lines, collapse = "\n"))
    
  })
  
  #add in original comments and whitespace
  processed_code <- paste0(unlist(clean_stan_annot$preceding_content),
                           unlist(clean_stan_annot$leading_ws),
                           newline_processed_code,
                           unlist(clean_stan_annot$trailing_ws),
                           unlist(clean_stan_annot$comments))
  
  #add in block comments
  rleb <- rle(parsed_lines$block_name)
  inds_to_insert_at <- cumsum(c(0, rleb$lengths[-length(rleb$lengths)])) + 1
  inds_to_insert_at <- inds_to_insert_at + 1:length(inds_to_insert_at) - 1
  block_names_to_insert <- parsed_lines$block_name[cumsum(rleb$lengths)]
  block_names_to_insert <- paste0("#### ", block_names_to_insert, " ####")
  block_names_to_insert <- sapply(block_names_to_insert, function(x)
    paste0(
      "\n#", paste0(rep("~", nchar(x)-2), collapse = ""), "#\n",
      x, "\n",
      "#", paste0(rep("~", nchar(x)-2), collapse = ""), "#\n",
      collapse = ""
    )
  )
  
  for(i in 1:length(inds_to_insert_at)){
    processed_code <- append(processed_code, block_names_to_insert[i], inds_to_insert_at[i]-1)  
  }
  
  #make sure we call data.table and initialize out at the start
  processed_code <- c("library(data.table)\nout <- list()\n", processed_code)
  
  #collapse into one line string
  collapsed_code <- paste0(processed_code, collapse = "\n")
  collapsed_code <- gsub("\n\n\n", "\n\n", collapsed_code)
  
  if(!is.na(output_file)){
    sink(output_file)
    cat(collapsed_code)
    sink()  
  }
  
  return(collapsed_code)
}

#### plotting functions ####

prettyScientificExpr <- function(numbers, digits = 2, pretty = T, scientific = T, exp_scientific = 4) {
  if(length(numbers) > 1){
    
    if(pretty) numbers <- pretty(numbers, length(numbers))
    
    if(scientific){
      out <- lapply(numbers, prettyScientificExpr, digits = digits)
      return(list(values = sapply(out, function(x) x[["values"]]), 
                  labels = lapply(out, function(x) x[["labels"]])))
      return()
    } else {
      return(list(values = numbers, 
                  labels = as.list(numbers)))
    }
    
  } else {
    
    split_number <- strsplit(formatC(numbers, format = "e", digits = digits), "e")[[1]]
    coefficient <- split_number[1]
    exponent <- as.numeric(split_number[2])
    
    
    
    if (exponent == 0) {
      result <- parse(text = coefficient)
    } else if(exponent >= exp_scientific || exponent < 0){
      result <- parse(text = paste0("list(", coefficient, " %*% 10^", exponent, ")"))
    } else {
      result <- parse(text = as.character(numbers))
    }
    
    return(list(values = numbers, 
                labels = result))
    
  }
  
}

# prettyScientificExpr(7:10)


add_continuous_legend <- function(colors, labels, positions, n_rect = 100, x = NA, y = NA, h = NA, w = NA, 
                                  vertical = TRUE, xpd = NA, n_labels = 5, 
                                  left_below = TRUE, log_scale = FALSE,
                                  n_digits = 3, pretty = T, scientific = T, main = NA) {
  
  usr <- par("usr")
  xyr <- diff(range(usr[1:2])) / par("pin")[1] / diff(range(usr[3:4])) * par("pin")[2]
  # Calculate plot dimensions if not provided
  if (is.na(h) || is.na(w)) {
    plot_width <- usr[2] - usr[1]
    plot_height <- usr[4] - usr[3]
  }
  
  # Adjust size for orientation
  if (vertical) {
    h <- ifelse(is.na(h), plot_height * 0.75, h)
    w <- ifelse(is.na(w), plot_width * 0.05, w)
  } else {
    h <- ifelse(is.na(h), plot_height * 0.05, h)
    w <- ifelse(is.na(w), plot_width * 0.75, w)
  }
  
  # Calculate legend position if not provided
  if (is.na(x) || is.na(y)) {
    if (vertical) {
      x <- usr[2] - w * 1.5
      y <- usr[4] - w * 0.5 / xyr
    } else {
      y <- usr[4] - h * 1.5
      x <- usr[2] - w - h * 1.5 * xyr
    }
  }
  
  # Adjust the color spectrum
  if (length(colors) != n_rect) {
    colors <- colorRampPalette(colors)(n_rect)
  }
  
  # Add rectangles
  rect_width <- w / n_rect
  rect_height <- h / n_rect
  
  for (i in 1:n_rect) {
    if (vertical) {
      rect_x_left <- x
      rect_x_right <- x + w
      rect_y_bottom <- y - (i-1) * rect_height
      rect_y_top <- y - i * rect_height
    } else {
      rect_x_left <- x + (i-1) * rect_width
      rect_x_right <- x + i * rect_width
      rect_y_bottom <- y
      rect_y_top <- y - h
    }
    
    rect(rect_x_left, rect_y_bottom, rect_x_right, rect_y_top, col = colors[i], border = NA, xpd = xpd)
  }
  rect(x, y-h, x+w, y, xpd = xpd)
  
  # Handling labels
  if (!all(is.na(labels)) && !all(is.na(positions))) {
    if (!all(is.na(n_labels)) && n_labels > 1) {
      
      if (log_scale) {
        min_val <- log10(min(labels))
        max_val <- log10(max(labels))
      } else {
        min_val <- min(labels)
        max_val <- max(labels)
      }
      max_pos <- positions[which.max(labels)]
      min_pos <- positions[which.min(labels)]
      diff_minmax <- max_val - min_val
      diff_scale <- max_pos - min_pos
      slope <- diff_minmax / diff_scale
      
      label_positions <- seq(0, 1, length.out = n_labels)
      label_relative_positions <- label_positions - min_pos
      extrapolated_values <- min_val + label_relative_positions * slope
      if (log_scale) {
        extrapolated_values <- 10^(extrapolated_values)
      }
      
      #now make pretty if we want that
      if(log_scale){
        order_mag <- (10^floor(log10(extrapolated_values)))
        lab_logmod10 <- extrapolated_values / order_mag
        if(pretty){
          prettylab_logmod10 <- ceiling(lab_logmod10 - 1E-6)
          prettylab_logmod10[length(prettylab_logmod10)] <- floor(lab_logmod10[length(lab_logmod10)] + 1E-6)
          labels_values <- prettyScientificExpr(prettylab_logmod10 * order_mag, n_digits, 
                                                pretty = F, scientific = scientific)
        } else {
          labels_values <- prettyScientificExpr(lab_logmod10 * order_mag, n_digits, 
                                                pretty = F, scientific = scientific)
        }
        
      } else {
        labels_values <- prettyScientificExpr(extrapolated_values, n_digits, 
                                              pretty = pretty, scientific = scientific)  
      }
      
      
      label_positions <- sapply(labels_values$values, function(xi){
        if(log_scale){
          (log10(xi) - min_val) / slope + min_pos
        } else {
          (xi - min_val) / slope + min_pos
        }
      })
      
      #check that none of the labels_vals are outside the range of the scale
      in_range <- label_positions >= 0 & label_positions <= 1
      labels_values$values <- labels_values$values[in_range]
      labels_values$labels <- labels_values$labels[in_range]
      label_positions <- label_positions[in_range]
      label_vals <- labels_values[["labels"]]
      
      # #if we got rid of boundary labels, add them back in?
      # interlab_dist <- 1 / n_labels
      # boundary_lab_relative_dists <- c(label_positions[1], 1 - tail(label_positions, 1)) / interlab_dist
      # if(any(boundary_lab_relative_dists > 0.5)){
      #   labpos_to_add <- c(0,1)[boundary_lab_relative_dists > 0.5]
      #   label_positions <- c(label_positions, labpos_to_add)
      #   labels_values
      # }
      
      #plot the labels
      lb <- left_below
      for (i in seq_along(label_vals)) {
        label_pos <- label_positions[i]
        label_val <- label_vals[[i]]
        if (vertical) {
          text_x <- x - (strwidth(label_val) / 2 + w / 2) * ifelse(lb, 1, -1) + ifelse(lb, 0, w)
          text_y <- y - h + h * label_pos
          segments(x0 = x + w * ifelse(lb, 0, 1), 
                   x1 = x + ifelse(lb, -w / 4, 5 * w / 4), 
                   y0 = text_y, 
                   y1 = text_y, 
                   xpd = xpd)
        } else {
          text_x <- x + label_pos * w
          text_y <- y - (h * ifelse(lb, 1, 0) + strheight(lb) / 2 + h / 2) * ifelse(lb, 1, -1)
          segments(x0 = text_x, 
                   x1 = text_x, 
                   y0 = y - h * ifelse(lb, 1, 0), 
                   y1 = y + ifelse(lb, - 5 * h / 4, h/ 4),
                   xpd = xpd)
        }
        text(text_x, text_y, label_val, xpd = xpd)
      }
    }
  }
  
  #write in the title
  if(!is.na(main)){
    text(x = x + w / 2, y = y + strheight(main) / 4 * 3, labels = main, xpd = xpd, font = 2)  
  }
  
}
