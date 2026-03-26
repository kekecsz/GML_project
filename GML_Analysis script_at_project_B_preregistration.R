### Data analysis code for the GML project

#########################################
#             Libraries                 #
#########################################
library(tidyverse)

#########################################
#          Custom functions             #
#########################################

### Function to see if a number is even

is_even <- function(x) x %% 2 == 0

### Function that pre-processes data that is already segmented per lab

preprocessor_1 <- function(data_pre){
  #################################################################
  #      Invalidate shams trials that come before reals trials    #
  #                  and extracting hits data                     #
  #################################################################
  
  session_IDs_with_pairs = NULL
  session_IDs_that_cannot_be_paired = NULL
  trial_pairs = data.frame(real_session_ID = NULL, sham_session_ID = NULL, hits_X = NULL, hits_O = NULL, hits_SX = NULL, hits_SO = NULL, lab_IDs = NULL)
  
  lab_IDs = unique(hits_data_valid$lab_ID)
  
  for(j in 1:length(lab_IDs)){
    nrow_real = hits_data_valid %>% filter(E_session_type == "real", lab_ID == lab_IDs[j]) %>% nrow()
    nrow_sham = hits_data_valid %>% filter(E_session_type == "sham", lab_ID == lab_IDs[j]) %>% nrow()
    
    if(nrow_real == 0){
      print(paste0("no real session in lab ", lab_IDs[j]))
      shams_without_pairs = hits_data_valid %>% filter(E_session_type == "sham", lab_ID == lab_IDs[j]) %>% pull(session_ID)
      session_IDs_that_cannot_be_paired = c(session_IDs_that_cannot_be_paired, shams_without_pairs)
    } else if(nrow_sham == 0){
      print(paste0("no sham session in lab ", lab_IDs[j]))
      reals_without_pairs = hits_data_valid %>% filter(E_session_type == "real", lab_ID == lab_IDs[j]) %>% pull(session_ID)
      session_IDs_that_cannot_be_paired = c(session_IDs_that_cannot_be_paired, reals_without_pairs)
    } else {
      for(i in 1:nrow_real){
        # only search for pairs in the dataset of trials that we don't yet have pairs for in the given lab
        print(paste0("finding pair for trial ", i, " in lab ", lab_IDs[j]))
        
        data_to_search_in = hits_data_valid %>% 
          filter(!session_ID %in% session_IDs_with_pairs,
                 !session_ID %in% session_IDs_that_cannot_be_paired,
                 lab_ID == lab_IDs[j])
        
        
        real_trial_dates = data_to_search_in %>% filter(E_session_type == "real") %>% pull(E_in_lab_experiment_start_time)
        
        first_real_trial_date = min(real_trial_dates)
        which_first_real_trial_date = which.min(real_trial_dates)
        real_session_ID = data_to_search_in %>% filter(E_session_type == "real") %>% slice(which_first_real_trial_date) %>% pull(session_ID)
        hits_X = data_to_search_in %>% filter(E_session_type == "real") %>% slice(which_first_real_trial_date) %>% pull(hits_X)
        hits_O = data_to_search_in %>% filter(E_session_type == "real") %>% slice(which_first_real_trial_date) %>% pull(hits_O)
        
        sham_trial_dates = data_to_search_in %>% filter(E_session_type == "sham") %>% pull(E_in_lab_experiment_start_time)
        
        ### Extract the index of the lowest of the sham_trial_dates that are higher than the first_real_trial_date
        which_first_higher_lowest_date = which_first_lowest(first_real_trial_date, sham_trial_dates)
        
        if(is.na(which_first_higher_lowest_date)){
          session_IDs_that_cannot_be_paired = c(session_IDs_that_cannot_be_paired, real_session_ID)
        }
        
        if(!is.na(which_first_higher_lowest_date)){
          sham_session_ID = data_to_search_in %>% filter(E_session_type == "sham") %>% slice(which_first_higher_lowest_date) %>% pull(session_ID)
          hits_SX = data_to_search_in %>% filter(E_session_type == "sham") %>% slice(which_first_higher_lowest_date) %>% pull(hits_X)
          hits_SO = data_to_search_in %>% filter(E_session_type == "sham") %>% slice(which_first_higher_lowest_date) %>% pull(hits_O)
          session_IDs_with_pairs = c(session_IDs_with_pairs, real_session_ID, sham_session_ID)
          
          trial_pairs = rbind(trial_pairs, data.frame(real_session_ID = real_session_ID, sham_session_ID = sham_session_ID, hits_X = hits_X, hits_O = hits_O, hits_SX = hits_SX, hits_SO = hits_SO, lab_IDs = lab_IDs[j]))
          
          ### check if there are any sham sessions that have happened earlier than the real session of interest
          first_higher_lowest_date = data_to_search_in %>% filter(session_ID == sham_session_ID) %>% pull(E_in_lab_experiment_start_time)
          
          if(sum(sham_trial_dates < first_higher_lowest_date) > 0) {
            which_sham_trial_dates_before_real = which(sham_trial_dates < first_higher_lowest_date)
            shams_without_pairs = data_to_search_in %>% filter(E_session_type == "sham") %>% slice(which_sham_trial_dates_before_real) %>% pull(session_ID)
            session_IDs_that_cannot_be_paired = c(session_IDs_that_cannot_be_paired, shams_without_pairs)
          }
          
          if((sum(sham_trial_dates > first_higher_lowest_date) & (length(real_trial_dates) == 1)) > 0) {
            
            which_are_remainder_sham_trials = which(sham_trial_dates > first_higher_lowest_date)
            shams_without_pairs = data_to_search_in %>% filter(E_session_type == "sham") %>% slice(which_are_remainder_sham_trials) %>% pull(session_ID)
            session_IDs_that_cannot_be_paired = c(session_IDs_that_cannot_be_paired, shams_without_pairs)
          }
        }
      }
    }
  }
  
  trial_pairs_by_lab = split(trial_pairs, trial_pairs$lab_IDs)
  
  return(trial_pairs_by_lab)
}
  
preprocessor_2 <- function(data, study_ns){

 ### extract row indices, hits, and study sample size for each "sub-study" 
  
  study_row_indices = NULL
  last_number = 0
  data_list = NULL
  
  for(i in 1:length(study_ns)){
    study_row_indices[[i]] = (last_number+1):(last_number+study_ns[i])
    last_number = last_number+study_ns[i]
    
    data_list[[i]] = data[study_row_indices[[i]],]

  }
  return(data_list)
}


preprocessor_3 <- function(data, study_ns){
  
  ### extract row indices, hits, and study sample size for each "sub-study" 
  
  study_row_indices = NULL
  last_number = 0
  data_list = NULL
  
  for(i in 1:length(study_ns)){
    study_row_indices[[i]] = (last_number+1):(last_number+study_ns[i])
    last_number = last_number+study_ns[i]
    
    data_study = data[study_row_indices[[i]],]
    data_study$lab_IDs = paste(unique(data_study$lab_IDs), collapse = "_")

    data_list[[i]] = data_study
    
  }
  return(data_list)
}


preprocessor_4 <- function(data){
  
  # Template sub-list
  template <- list(
    lab_ID = NA,
    results = data.frame(X = NA, O = NA, SX = NA, SO = NA),
    study_n = data.frame(X = NA, O = NA, SX = NA, SO = NA),
    num_hits = data.frame(X = NA, O = NA, SX = NA, SO = NA),
    p_value = data.frame(X = NA, O = NA, SX = NA, SO = NA),
    Z_score = data.frame(X = NA, O = NA, SX = NA, SO = NA)
  )
  
  # Create sub-lists
  results_list <- replicate(length(data), template, simplify = FALSE)
  
  ### extract hits, and study sample size for each "sub-study" 

  for(i in 1:length(data)){
    results_list[[i]][["lab_ID"]] = data[[i]]$lab_ID[1]
    
    results_list[[i]][["results"]] = data.frame(
      X = data[[i]][, "hits_X"],
      O = data[[i]][, "hits_O"],
      SX = data[[i]][, "hits_SX"],
      SO = data[[i]][, "hits_SO"]
    )
    
    results_list[[i]][["study_n"]] = data.frame(X = nrow(data[[i]]), O = nrow(data[[i]]), 
                                                SX = nrow(data[[i]]), SO = nrow(data[[i]]))
    results_list[[i]][["num_hits"]] = data.frame(
      X = sum(results_list[[i]][["results"]][["X"]]),
      O = sum(results_list[[i]][["results"]][["O"]]),
      SX = sum(results_list[[i]][["results"]][["SX"]]),
      SO = sum(results_list[[i]][["results"]][["SO"]])
    )
  }
  
  return(results_list)
}

### Extract the index of the lowest of the sham_trial_dates that are higher than the first_real_trial_date
which_first_lowest <- function(real_trial_date, sham_trial_dates){
  if(sum(real_trial_date < sham_trial_dates) == 0){
    which_first_higher_lowest_value = NA
  } else {
    first_higher_lowest_value = min(sham_trial_dates[real_trial_date < sham_trial_dates])
    which_first_higher_lowest_value = which(first_higher_lowest_value == sham_trial_dates)
  }
  
  return(which_first_higher_lowest_value)
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

global_analysis_function = function(results_list, prob_of_significant_hits){

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
                                                                         p = prob_of_significant_hits, 
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

### Function that calculates the minimum number of successes that lead to a significant finding,
### given the number of trials

min_successes_for_significance <- function(n_rolls, p = 0.25, alpha = 0.05, alternative = "greater") {
  for (k in 0:n_rolls) {
    p_val <- binom.test(k, n_rolls, p = p, alternative = alternative)$p.value
    if (p_val < alpha) return(k)
  }
  return(NA)  # No significant result possible
}

### Wrapper function that does all preprocessing and data analysis at once to get the final confirmatory analysis data

do_GML_confirmatory_analysis = function(data){
  
  data_segmented_by_lab_preprocessed_1 = preprocessor_1(data)
  
  cycle_through_labs_for_n51 = c(
    names(data_segmented_by_lab_preprocessed_1)[sapply(data_segmented_by_lab_preprocessed_1, nrow) >= 51*1],
    names(data_segmented_by_lab_preprocessed_1)[sapply(data_segmented_by_lab_preprocessed_1, nrow) >= 51*2],
    names(data_segmented_by_lab_preprocessed_1)[sapply(data_segmented_by_lab_preprocessed_1, nrow) >= 51*3],
    names(data_segmented_by_lab_preprocessed_1)[sapply(data_segmented_by_lab_preprocessed_1, nrow) >= 51*4],
    names(data_segmented_by_lab_preprocessed_1)[sapply(data_segmented_by_lab_preprocessed_1, nrow) >= 51*5],
    names(data_segmented_by_lab_preprocessed_1)[sapply(data_segmented_by_lab_preprocessed_1, nrow) >= 51*6],
    names(data_segmented_by_lab_preprocessed_1)[sapply(data_segmented_by_lab_preprocessed_1, nrow) >= 51*7],
    names(data_segmented_by_lab_preprocessed_1)[sapply(data_segmented_by_lab_preprocessed_1, nrow) >= 51*8]
  )[1:8]

  
  results_list_preprocessed = NULL
  
  for(i in 1:length(data_segmented_by_lab_preprocessed_1)){
    n51s_for_this_study = rep(51, sum(cycle_through_labs_for_n51%in%names(data_segmented_by_lab_preprocessed_1)[[i]]))
    remainder = (nrow(data_segmented_by_lab_preprocessed_1[[i]])-sum(n51s_for_this_study)) %% 50
    n50s_for_this_lab = rep(50, (nrow(data_segmented_by_lab_preprocessed_1[[i]])-sum(n51s_for_this_study)-remainder)/50)
    n_of_studies_for_this_lab = c(n51s_for_this_study, n50s_for_this_lab, if(remainder == 0){NULL}else{remainder})

        
    results_list_preprocessed_lab = preprocessor_2(data_segmented_by_lab_preprocessed_1[[i]], n_of_studies_for_this_lab)
    results_list_preprocessed = c(results_list_preprocessed, results_list_preprocessed_lab)
  }
  
  trials_n50ormore = results_list_preprocessed[sapply(results_list_preprocessed, nrow) >= 50]
  
  remaining_trials=do.call(rbind, results_list_preprocessed[sapply(results_list_preprocessed, nrow) < 50])
  
  if(is.null(remaining_trials)){
    results_list_preprocessed_2 = trials_n50ormore
  } else {
    remainder = nrow(remaining_trials) %% 50
    n50s_for_remaining_trials = rep(50, (nrow(remaining_trials)-remainder)/50)
    n_of_studies_for_remaining_trials = c(n50s_for_remaining_trials, if(remainder == 0){NULL}else{remainder})
    
    remaining_studies_list = preprocessor_3(remaining_trials, n_of_studies_for_remaining_trials)
    
    results_list_preprocessed_2 = c(trials_n50ormore, remaining_studies_list)
    
  }

  results_list_preprocessed_3 = preprocessor_4(results_list_preprocessed_2)
  
  
  results_list = per_study_analysis_function(results_list_preprocessed_3)
  
  
  study_sizes = sapply(results_list, function(x) nrow(x$results))
  
  probability_grid = data.frame(study_size = NA, number_of_studies_with_this_study_size = NA, min_num_successes_to_significance = NA, p_value_at_min_num_successes_to_significance = NA)
  
  for(i in 1:length(unique(study_sizes))){
    probability_grid[i,"study_size"] = unique(study_sizes)[i]
    probability_grid[i,"number_of_studies_with_this_study_size"] = table(study_sizes)[names(table(study_sizes)) == unique(study_sizes)[i]]
    min_num_successes_to_significance = min_successes_for_significance(unique(study_sizes)[i])
    probability_grid[i,"min_num_successes_to_significance"] = min_num_successes_to_significance
    probability_grid[i,"p_value_at_min_num_successes_to_significance"] = pbinom((min_num_successes_to_significance-1), size = unique(study_sizes)[i], prob = 0.25, lower.tail = FALSE)
  }

  prob_of_significant_hits = sum(probability_grid[,"number_of_studies_with_this_study_size"] * probability_grid[,"p_value_at_min_num_successes_to_significance"])/sum(probability_grid[,"number_of_studies_with_this_study_size"])

  global_results_table = global_analysis_function(results_list, prob_of_significant_hits)
  
  return(global_results_table)
}

#########################################
#               Load data               #
#########################################


### import data into an object
### in the example below the analysis runs on simulated data, this needs to be replaced with the real data

hits_data_raw = read.csv("https://raw.githubusercontent.com/kekecsz/GML_project/refs/heads/main/simulated_hits_data.csv")
hits_data_raw = hits_data_raw[order(hits_data_raw$E_in_lab_experiment_start_time), ]

### Exclude invalid sessions

hits_data_valid = hits_data_raw[hits_data_raw$valid == 1,]

#########################################
#             Run analysis              #
#########################################

### run the analysis
results = do_GML_confirmatory_analysis(hits_data_valid)

results
