make_latex_entry <- function(a, b) {
  if (is.na(a) || is.na(b)) {
    return(NA_character_)
  } else {
    sprintf("$%.2f\\pm%.3f$", a, b)
  }
}



MODEL_SELECTION <- function(files_path, n_list, tau_list, cop_TYPES, all_criteria, output_filename){
   # Get all file names
  all_files <- list.files(files_path, pattern = "\\.csv$", full.names = TRUE)
  
  # Open main output file
  sink(output_filename)
  
  for (tau in tau_list){
    for (n in n_list){
      # Regular expression for given tau and n
      pattern <- paste0("tau_", tau, "_sample_size_", n, "\\.csv$")
      
      # Selecting files that match a pattern
      matching_files <- grep(pattern, all_files, value = TRUE)
      
      # Check if the file, which follows the pattern, exists
      if (length(matching_files) > 0) {
        # Load and merge all matching files into one data.frame
        df_n_tau <- lapply(matching_files, read.csv, check.names = FALSE) %>% bind_rows()
        
        # create empty dataframe of results
        df_result <- skeleton_MODEL_SELECTION(cop_TYPES, all_criteria)
        
        for (i in 1:nrow(df_n_tau)) {
          true_cop <- df_n_tau$cop_name[i]
          for (crit in all_criteria) {
            # prefixes like "AIC_", "xv1_", "xvK2_"
            prefix <- paste0(crit, "_")
            
            # columns matching the given criterion
            crit_cols <- grep(paste0("^", prefix), names(df_n_tau), value = TRUE)
            
            # name of selected copula
            selected_cop <- sub(prefix, "", crit_cols[which.max(df_n_tau[i, crit_cols])])
            
            # find the matching row in the df_result and increment the count
            row_index <- which(df_result$d.cop == true_cop & df_result$IC == crit)
            df_result[row_index, selected_cop] <- df_result[row_index, selected_cop] + 1
          }
        }
        
        caption_text <- sprintf("Copula selection using different information criteria ($n = %d$, $\\tau = %.2f$)", n, tau)
        
        # Adding dashlines to your table
        line_positions <- seq(10, nrow(df_result) - 1, by = 10)
        add_to_row <- list()
        add_to_row$pos <- as.list(line_positions)
        add_to_row$command <- rep("\\hdashline\n", length(line_positions)) 
        
        display_names <- c(
          "AIC" = "AIC",
          "xv1" = "$\\text{xv}_{1}$",
          "xv_CIC" = "$\\text{xv}\\_\\text{CIC}$",
          "xvK_2" = "$K_{n}=2$",
          "xvK_5" = "$K_{n}=5$",
          "xvK_10" = "$K_{n}=10$",
          "xvK_n:2" = "$K_{n}=\\tfrac{n}{2}$",
          "xvK_n:logn" = "$K_{n}=\\tfrac{n}{\\text{log}(n)}$",
          "xvK_n05" = "$K_{n}=\\sqrt{n}$",
          "xvK_n05:logn" = "$K_{n}=\\tfrac{\\sqrt{n}}{\\text{log}(n)}$"
        )
        
        df_result$IC <- display_names[as.character(df_result$IC)]
        
        print(
          xtable(df_result, caption = caption_text, digits = c(0, 0, 0, 0, 0, 0, 0, 0)),
          include.rownames = FALSE,
          caption.placement = "bottom",
          sanitize.text.function = identity,
          comment = FALSE,
          add.to.row = add_to_row,
          table.placement = "htbp"
        )
        cat("\n\n")  # space between tables
      }
    }
  }
  sink()
}



skeleton_MODEL_SELECTION <- function(cop_TYPES, all_criteria){
  # Create an empty data.frame with all combinations of copula names and information criteria
  df <- expand.grid(d.cop = cop_TYPES, IC = all_criteria, stringsAsFactors = FALSE)
  
  # Add 0 for each copula
  for (cop in cop_TYPES) {
    df[[cop]] <- 0
  }
  df$d.cop <- factor(df$d.cop, levels = cop_TYPES)
  df$IC <- factor(df$IC, levels = all_criteria)
  
  # Sort by both columns
  df <- df[order(df$d.cop, df$IC), ]
  
  # Reset row names
  rownames(df) <- NULL
  return(df)
}



COINC_PERC <- function(files_path, n_list, tau_list, cop_TYPES, all_criteria, criterion1, criterion2){
  # Get all file names
  all_files <- list.files(files_path, pattern = "\\.csv$", full.names = TRUE)
  
  agreement_matrix <- matrix(NA, nrow = length(n_list), ncol = length(tau_list),
                             dimnames = list(n_list, tau_list))
  
  for (n in n_list){
    
    for (tau in tau_list){
      
      # Regular expression for given tau and n
      pattern <- paste0("tau_", tau, "_sample_size_", n, "\\.csv$")
      
      # Selecting files that match a pattern
      matching_files <- grep(pattern, all_files, value = TRUE)
      
      # Check if the file, which follows the pattern, exists
      if (length(matching_files) > 0) {
        
        # Load and merge all matching files into one data.frame
        df_n_tau <- lapply(matching_files, read.csv) %>% bind_rows()
        N <- nrow(df_n_tau)
        
        # number of times when two criteria selected the same copula
        agree_count <- 0
        
        for (i in 1:N) {
          
          # which copula is selected by criterion1
          criterion1_cols <- grep(paste0("^", criterion1, "_"), names(df_n_tau), value = TRUE)
          criterion1_selected <- sub(paste0(criterion1, "_"), "", criterion1_cols[which.max(df_n_tau[i, criterion1_cols])])
          
          # which copula is selected by criterion2
          criterion2_cols <- grep(paste0("^", criterion2, "_"), names(df_n_tau), value = TRUE)
          criterion2_selected <- sub(paste0(criterion2, "_"), "", criterion2_cols[which.max(df_n_tau[i, criterion2_cols])])
          
          # add 1 if two criteria selected the same copula
          if (criterion1_selected == criterion2_selected) {
            agree_count <- agree_count + 1
          }
        }
        
        # result write to matrix
        agreement_matrix[as.character(n), as.character(tau)] <- agree_count
      }
    }
    
  }
  
  options(digits = 15) 
  
  agreement_df <- as.data.frame(agreement_matrix)
  All <-  100*(rowSums(agreement_df, na.rm = TRUE) / (rowSums(!is.na(agreement_df)) * N))
  
  agreement_df <- 100*(agreement_df/N)
  agreement_df$All <- All
  
  is_num <- sapply(agreement_df, is.numeric)
  
  agreement_df[is_num] <- lapply(
    agreement_df[is_num],
    function(x) ifelse(is.na(x), NA, sprintf("%.3f", x))
  )
  return(agreement_df)
}



SE_COINC_PERC <- function(files_path, n_list, tau_list, cop_TYPES, all_criteria, criterion1, criterion2, num_repl){
  
  # Get all file names
  all_files <- list.files(files_path, pattern = "\\.csv$", full.names = TRUE)
  
  se_matrix <- matrix(NA, nrow = length(n_list), ncol = length(tau_list),
                      dimnames = list(n_list, tau_list))
  se_matrix <- cbind(se_matrix, All = rep(NA, length(n_list)))
  
  N <- num_repl*length(cop_TYPES)
  
  for (n in n_list){
    
    all_vector <- rep(NA, N*length(tau_list))
    
    for (tau in tau_list){
      
      # Regular expression for given tau and n
      pattern <- paste0("tau_", tau, "_sample_size_", n, "\\.csv$")
      
      # Selecting files that match a pattern
      matching_files <- grep(pattern, all_files, value = TRUE)
      
      if (length(matching_files) > 0) {
        # Load and merge into one data.frame
        df_n_tau <- lapply(matching_files, read.csv) %>% bind_rows()
        
        cell_vector <- rep(0, N)
        
        for (i in 1:N) {
          
          idx_all_vector <- (which(tau_list == tau)-1)*N + i
          all_vector[idx_all_vector] <- 0
          
          # which copula is selected by criterion1
          criterion1_cols <- grep(paste0("^", criterion1, "_"), names(df_n_tau), value = TRUE)
          criterion1_selected <- sub(paste0(criterion1, "_"), "", criterion1_cols[which.max(df_n_tau[i, criterion1_cols])])
          
          # which copula is selected by criterion2
          criterion2_cols <- grep(paste0("^", criterion2, "_"), names(df_n_tau), value = TRUE)
          criterion2_selected <- sub(paste0(criterion2, "_"), "", criterion2_cols[which.max(df_n_tau[i, criterion2_cols])])
          
          if (criterion1_selected == criterion2_selected) {
            cell_vector[i] <- 1
            all_vector[idx_all_vector] <- 1
          }
        }
        
        # write results to the matrix
        se_matrix[as.character(n), as.character(tau)] <- sd(cell_vector)/sqrt(length(cell_vector))
      }
    }
    
    all_vector <- na.omit(all_vector)
    se_matrix[as.character(n), "All"] <- sd(all_vector)/sqrt(length(all_vector))
  }
  
  options(digits = 15) 
  
  se_df <- 100*(qnorm(0.975)*as.data.frame(se_matrix))
  
  is_num <- sapply(se_df, is.numeric)
  
  se_df[is_num] <- lapply(
    se_df[is_num],
    function(x) ifelse(is.na(x), NA, sprintf("%.4f", x))
  )
  
  return(se_df)
}



skeleton_HIT_RATES <- function(all_criteria, cop_TYPES){
  df <- data.frame(
    IC = all_criteria,
    matrix(0, nrow = length(all_criteria), ncol = length(cop_TYPES))
  )
  colnames(df)[-1] <- cop_TYPES
  return(df)
}



HIT_RATES <- function(files_path, n_list, tau_list, cop_TYPES, all_criteria, output_filename){
  
  # Get all file names
  all_files <- list.files(files_path, pattern = "\\.csv$", full.names = TRUE)
  
  sink(output_filename)
  
  for (tau in tau_list){
    for (n in n_list){
      
      # Regular expression for given tau and n
      pattern <- paste0("tau_", tau, "_sample_size_", n, "\\.csv$")
      
      # Selecting files that match a pattern
      matching_files <- grep(pattern, all_files, value = TRUE)
      
      if (length(matching_files) > 0) {
        
        # Load and merge into one data.frame
        df_n_tau <- lapply(matching_files, read.csv) %>% bind_rows()
        
        # dataframe of hit rates
        df_hit_rates <- skeleton_HIT_RATES(all_criteria, cop_TYPES)
        
        for (i in 1:nrow(df_n_tau)) {
          true_cop <- df_n_tau$cop_name[i]
          
          for (crit in all_criteria) {
            # prefixes like "AIC_", "xv1_", "xv_nv_"
            prefix <- paste0(crit, "_")
            
            # columns matching the given criterion
            crit_cols <- grep(paste0("^", prefix), names(df_n_tau), value = TRUE)
            
            # name of selected copula
            selected_cop <- sub(prefix, "", crit_cols[which.max(df_n_tau[i, crit_cols])])
            
            if (selected_cop == true_cop){
              row_index <- which(df_hit_rates$IC == crit)
              df_hit_rates[row_index, selected_cop] <- df_hit_rates[row_index, selected_cop] + 1
            }
            #print(c(crit, i))
            
          }
        }
        
        # dataframe of standard errors
        df_se_hit_rates <- skeleton_HIT_RATES(all_criteria, cop_TYPES)
        
        for (cop in cop_TYPES){
          df_cop_n_tau <- df_n_tau[df_n_tau$cop_name == cop, ]
          
          my_list <- setNames(lapply(all_criteria, function(x) rep(0, nrow(df_cop_n_tau))), all_criteria)
          
          for (j in 1:nrow(df_cop_n_tau)) {
            
            for (crit in all_criteria) {
              # prefixes like "AIC_", "xv1_", "xv_nv_"
              prefix <- paste0(crit, "_")
              
              # columns matching the given criterion
              crit_cols <- grep(paste0("^", prefix), names(df_cop_n_tau), value = TRUE)
              
              # name of selected copula
              selected_cop <- sub(prefix, "", crit_cols[which.max(df_cop_n_tau[j, crit_cols])])
              
              if (selected_cop == cop){
                my_list[[crit]][j] <- 1 
              }
            }
          }
          
          df_se_hit_rates[,cop] <- unname(unlist(lapply(my_list, function(x) sd(x) / sqrt(length(x)))))
        }
    
        options(digits = 15)  
        
        df_hit_rates <- 100*(df_hit_rates[,-1]/(nrow(df_n_tau)/length(cop_TYPES)))
        df_se_hit_rates <- 100*(df_se_hit_rates[,-1]*qnorm(0.975))
        
        # Firstly, save resulting dataframes as strings
        
        df_hit_rates_STRING <- as.data.frame(lapply(
          df_hit_rates,
          function(x) ifelse(is.na(x), NA, sprintf("%.4f", x))
        ))
        
        df_se_hit_rates_STRING <- as.data.frame(lapply(
          df_se_hit_rates,
          function(x) ifelse(is.na(x), NA, sprintf("%.4f", x))
        ))
        
        # Secondly, save resulting dataframes as numeric
        
        df_hit_rates_NUMERIC <- round(
          apply(df_hit_rates_STRING, 2, function(x) as.numeric(as.character(x))),
          digits = 2
        )
        
        df_se_hit_rates_NUMERIC <- round(
          apply(df_se_hit_rates_STRING, 2, function(x) as.numeric(as.character(x))),
          digits = 2
        )
        
        n_tau_numeric <- matrix(
          mapply(make_latex_entry,
                 as.matrix(df_hit_rates_NUMERIC),
                 as.matrix(df_se_hit_rates_NUMERIC)),
          nrow = nrow(df_hit_rates_NUMERIC),
          ncol = ncol(df_hit_rates_NUMERIC),
          byrow = FALSE,
          dimnames = list(rownames(df_hit_rates_NUMERIC), colnames(df_hit_rates_NUMERIC))
        )
        
        # Combine with IC column
        n_tau <- cbind(IC = c("AIC", "$\\text{xv}_{1}$", "$\\text{xv}_\\text{CIC}$", "$K=2$", "$K=5$",
                              "$K=10$", "$K_{n} = \\frac{n}{2}$", "$K_{n} = \\frac{n}{\\log(n)}$", "$K_{n}=\\sqrt{n}}$",
                              "$K_{n} = \\frac{\\sqrt{n}}{\\log(n)}$"), n_tau_numeric)
        
        caption_text <- sprintf("Hit rates ($n = %d$, $\\tau = %.2f$), with 95 \\%% confidence intervals (all values multiplied by 100).", n, tau)
        
        print(xtable(n_tau, align = rep("c", ncol(n_tau) + 1), caption = caption_text), 
              sanitize.text.function = identity, 
              caption.placement = "bottom",
              include.rownames = FALSE, 
              na.string = "")
        
        cat("\n\n")  # empty space between tables
      }
    }
  }
  sink()
}
