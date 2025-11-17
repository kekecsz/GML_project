### Data analysis code for the GML project

#########################################
#          Custom functions             #
#########################################

### Function that pre-processes data that is already segmented per lab

preprocessor <- function(data){
  remainder = nrow(data) %% 50
  last_study_n = 50 + remainder
  study_ns = c(rep(50, (nrow(data) - last_study_n)/50), last_study_n)
  
  # Template sub-list
  template <- list(
    lab_ID = NA,
    study_row_indices = NA,
    results = data.frame(X = NA, O = NA, SX = NA, SO = NA),
    study_n = data.frame(X = NA, O = NA, SX = NA, SO = NA),
    num_hits = data.frame(X = NA, O = NA, SX = NA, SO = NA),
    p_value = data.frame(X = NA, O = NA, SX = NA, SO = NA),
    Z_score = data.frame(X = NA, O = NA, SX = NA, SO = NA)
  )
  
  # Create sub-lists
  results_list <- replicate(length(study_ns), template, simplify = FALSE)
  
  ### extract row indices, hits, and study sample size for each "sub-study" 
  last_number = 0
  for(i in 1:length(study_ns)){
    results_list[[i]][["lab_ID"]] = data$lab_ID[1]
    results_list[[i]][["study_row_indices"]] = (last_number+1):(last_number+study_ns[i])
    last_number = last_number+study_ns[i]
    
    results_list[[i]][["results"]] = data.frame(
      X = data[results_list[[i]][["study_row_indices"]], "hits_X"],
      O = data[results_list[[i]][["study_row_indices"]], "hits_O"],
      SX = data[results_list[[i]][["study_row_indices"]], "hits_SX"],
      SO = data[results_list[[i]][["study_row_indices"]], "hits_SO"]
    )
    
    results_list[[i]][["study_n"]] = data.frame(X = study_ns[i], O = study_ns[i], 
                                                SX = study_ns[i], SO = study_ns[i])
    results_list[[i]][["num_hits"]] = data.frame(
      X = sum(results_list[[i]][["results"]][["X"]]),
      O = sum(results_list[[i]][["results"]][["O"]]),
      SX = sum(results_list[[i]][["results"]][["SX"]]),
      SO = sum(results_list[[i]][["results"]][["SO"]])
    )
  }
  
  return(results_list)
}

### Function calculating Z scores from number of hits, number of trials, and H0 probability

calculate_z_score <- function(x, n, p0) {
  # Expected number of hits under H0
  expected_hits <- n * p0
  
  # Standard error under H0
  se <- sqrt(n * p0 * (1 - p0))
  
  # Z-score (normal approximation)
  z <- (x - expected_hits) / se
  
  return(z)
}

### Function calculating weighted Stouffer's Z from Z-scores

calculate_stouffer_z <- function(z_scores, ns) {
  z_scores = as.numeric(z_scores)
  ns = as.numeric(ns)
  Z_combined <- sum(z_scores * sqrt(ns))/sqrt(sum(ns))
  return(Z_combined)
}

### Function that analyzes every sub-study to get binomial test p-values and Z scores

per_study_analysis_function <- function(results_list){
  conditions <- c("X", "O", "SX", "SO")
  
  for(i in 1:length(results_list)){
    for(cond in conditions){
      results_list[[i]][["p_value"]][cond] <- round(binom.test(
        x = as.numeric(results_list[[i]][["num_hits"]][cond]), 
        n = as.numeric(results_list[[i]][["study_n"]][cond]), 
        p = 0.25, 
        alternative = "greater"
      )$p.value, 2)
      
      results_list[[i]][["Z_score"]][cond] <- round(calculate_z_score(
        x = as.numeric(results_list[[i]][["num_hits"]][cond]), 
        n = as.numeric(results_list[[i]][["study_n"]][cond]), 
        p0 = 0.25), 2)
    }
  }
  
  return(results_list)
}

### Function that analyzes the per lab results to get global confirmatory analysis data

global_analysis_function = function(results_list, prob_of_18_hits){

  global_results_table = data.frame(outcome_name = NA, X = NA, O = NA, SX = NA, SO = NA)
  global_results_table[1:6,"outcome_name"] = c("sig_binomtests", 
                                            "binom_of_num_sig_discoveries_p",
                                            "hits_19plus",
                                            "global_hits",
                                            "global_binom_p",
                                            "stouffer_z")
                                            
  
  perstudy_binom_p_values = t(sapply(results_list, function(x) x$p_value))
  num_sig_binomtests = apply(perstudy_binom_p_values, 2, function(x) sum(x<=0.05))
  global_results_table[global_results_table$outcome_name=="sig_binomtests",c("X", "O", "SX", "SO")] = num_sig_binomtests
  
  binom_of_num_sig_discoveries_p_values = sapply(num_sig_binomtests, function(x) round(binom.test(x, 
                                                                         n = length(results_list), 
                                                                         p = prob_of_18_hits, 
                                                                         alternative = "greater")$p.value, 2))
  global_results_table[global_results_table$outcome_name=="binom_of_num_sig_discoveries_p",c("X", "O", "SX", "SO")] = binom_of_num_sig_discoveries_p_values
  
  num_perstudy_hits = t(sapply(results_list, function(x) x$num_hits))
  num_19plus_hits = apply(num_perstudy_hits, 2, function(x) sum(x>18))
  global_results_table[global_results_table$outcome_name=="hits_19plus",c("X", "O", "SX", "SO")] = num_19plus_hits

  hits_data <- do.call(rbind, lapply(results_list, function(x) x$results))
  
  num_global_hits = apply(hits_data, 2, sum)
  global_results_table[global_results_table$outcome_name=="global_hits",c("X", "O", "SX", "SO")] = num_global_hits

  global_binom_p_values = sapply(num_global_hits, function(x) round(binom.test(x, 
                                                                         n = nrow(hits_data), 
                                                                         p = 0.25, 
                                                                         alternative = "greater")$p.value, 2))
  global_results_table[global_results_table$outcome_name=="global_binom_p",c("X", "O", "SX", "SO")] = global_binom_p_values

  
  ns = data.frame(t(sapply(results_list, function(x) x$study_n)))
  z_scores = data.frame(t(sapply(results_list, function(x) x$Z_score)))
  
  stouffer_z <- round(mapply(calculate_stouffer_z, 
                       z_scores, 
                       ns, 
                       SIMPLIFY = TRUE), 2)

  global_results_table[global_results_table$outcome_name=="stouffer_z",c("X", "O", "SX", "SO")] = stouffer_z

  return(global_results_table)
}

### Wrapper function that does all preprocessing and data analysis at once to get the final confirmatory analysis data

do_GML_confirmatory_analysis = function(data){
  data_segmented_by_lab = split(data, data$lab_ID)
  
  results_list_preprocessed = NULL
  
  for(i in 1:length(data_segmented_by_lab)){
    results_list_preprocessed_lab = preprocessor(data_segmented_by_lab[[i]])
    results_list_preprocessed = c(results_list_preprocessed, results_list_preprocessed_lab)
  }
  
  results_list = per_study_analysis_function(results_list_preprocessed)
  
  
  p_value_50 <- pbinom(18, size = 50, prob = 0.25, lower.tail = FALSE)
  p_value_50
  p_value_51 <- pbinom(18, size = 51, prob = 0.25, lower.tail = FALSE)
  p_value_51
  
  
  prob_of_18_hits = (20*p_value_50+8*p_value_51)/28
  
  global_results_table = global_analysis_function(results_list, prob_of_18_hits)
  
  return(global_results_table)
}

#########################################
#               Load data               #
#########################################


### save data into an object
### in the example below the analysis runs on simulated data, this needs to be replaced with the real data

hits_data = read.csv("https://raw.githubusercontent.com/kekecsz/GML_project/refs/heads/main/simulated_hits_data.csv")


#########################################
#             Run analysis              #
#########################################

### run the analysis
results = do_GML_confirmatory_analysis(hits_data)


results
