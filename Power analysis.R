### Simulation code to determine power and alpha for the GML project



### Packages 

library(pbapply) # for progress bar for long simulations

#########################################
#          Custom functions             #
#########################################


### functions for data simulation

lab_name_sampler <- function(n_trials, n_labs) {
  lab_post_fixes <- 1:n_labs
  urn_contents <- paste0("lab_", lab_post_fixes)
  n_complete_urns <- n_trials %/% length(urn_contents)
  n_remaining <- n_trials %% length(urn_contents)
  
  # Generate complete urns
  complete <- replicate(n_complete_urns, sample(urn_contents), simplify = FALSE)
  complete <- unlist(complete)
  
  # Add partial urn if needed
  if (n_remaining > 0) {
    partial <- sample(urn_contents, n_remaining)
    result <- c(complete, partial)
  } else {
    result <- complete
  }
  
  return(result)
}

data_simulation <- function(n_trials, n_labs, prob_X,prob_O,prob_SX,prob_SO){
  hits_data = data.frame(
    hits_X = rbinom(n_trials, 1, prob = prob_X),
    hits_O = rbinom(n_trials, 1, prob = prob_O),
    hits_SX = rbinom(n_trials, 1, prob = prob_SX),
    hits_SO = rbinom(n_trials, 1, prob = prob_SO),
    lab_ID = lab_name_sampler(n_trials, n_labs)
  )
  return(hits_data)
}

### Function to test success criteria

simulation_function = function(n_trials_first_research, n_trials_replication_research, n_labs, prob_X, prob_O, prob_SX, prob_SO){
  
  hits_data = data_simulation(n_trials = n_trials_first_research, n_labs = n_labs, prob_X = prob_X,prob_O = prob_O,prob_SX = prob_SX,prob_SO = prob_SO)
  
  global_results_table = do_GML_confirmatory_analysis(hits_data)
  
  global_results_first_research = NULL
  global_results_first_research["DR_triggered_in_previous_study"] = 0
  global_results_first_research["HR_triggered_in_previous_study"] = 0
  global_results_first_research["same_outcome_triggered_as_in_previous_study"] = 0
  
  if(global_results_table[global_results_table$outcome_name == "binom_of_num_sig_discoveries_p", "X"]<=0.05){global_results_first_research["DR_triggered"] = 1} else {global_results_first_research["DR_triggered"] = 0}
  if(global_results_table[global_results_table$outcome_name == "stouffer_z", "X"]>=1.64){global_results_first_research["HR_triggered"] = 1} else {global_results_first_research["HR_triggered"] = 0}
  
  effect_present_in_X = (global_results_first_research["DR_triggered"] == 1) | (global_results_first_research["HR_triggered"] == 1)
  global_results_first_research["effect_present_in_X"] = sum(effect_present_in_X)
  
  global_results_first_research["number_of_effects_in_control_experiments"] = sum(c(global_results_table[global_results_table$outcome_name == "binom_of_num_sig_discoveries_p", "O"]<=0.05, global_results_table[global_results_table$outcome_name == "stouffer_z", "O"]>=1.64, 
                                                                                    global_results_table[global_results_table$outcome_name == "binom_of_num_sig_discoveries_p", "SX"]<=0.05, global_results_table[global_results_table$outcome_name == "stouffer_z", "SX"]>=1.64, 
                                                                                    global_results_table[global_results_table$outcome_name == "binom_of_num_sig_discoveries_p", "SO"]<=0.05, global_results_table[global_results_table$outcome_name == "stouffer_z", "SO"]>=1.64))
  
  global_results_first_research["was_replication"] = 0
  global_results_first_research["replication_success"] = 0
  
  if(global_results_first_research["effect_present_in_X"] == 0){return(global_results_first_research)} else {
    
    hits_data_replication = data_simulation(n_trials = n_trials_replication_research, n_labs = n_labs, prob_X = prob_X,prob_O = prob_O,prob_SX = prob_SX,prob_SO = prob_SO)
    global_results_table_replication = do_GML_confirmatory_analysis(hits_data_replication)
    
    global_results_replication = NULL
    global_results_replication["DR_triggered_in_previous_study"] = 0
    global_results_replication["HR_triggered_in_previous_study"] = 0
    global_results_replication["same_outcome_triggered_as_in_previous_study"] = 0
    
    if(global_results_table_replication[global_results_table_replication$outcome_name == "binom_of_num_sig_discoveries_p", "X"]<=0.05){global_results_replication["DR_triggered"] = 1} else {global_results_replication["DR_triggered"] = 0}
    if(global_results_table_replication[global_results_table_replication$outcome_name == "stouffer_z", "X"]>=1.64){global_results_replication["HR_triggered"] = 1} else {global_results_replication["HR_triggered"] = 0}
    
    effect_present_in_X = (global_results_replication["DR_triggered"] == 1) | (global_results_replication["HR_triggered"] == 1)
    global_results_replication["effect_present_in_X"] = sum(effect_present_in_X)
    
    global_results_replication["number_of_effects_in_control_experiments"] = sum(c(global_results_table_replication[global_results_table_replication$outcome_name == "binom_of_num_sig_discoveries_p", "O"]<=0.05, global_results_table_replication[global_results_table_replication$outcome_name == "stouffer_z", "O"]>=1.64, 
                                                                                   global_results_table_replication[global_results_table_replication$outcome_name == "binom_of_num_sig_discoveries_p", "SX"]<=0.05, global_results_table_replication[global_results_table_replication$outcome_name == "stouffer_z", "SX"]>=1.64, 
                                                                                   global_results_table_replication[global_results_table_replication$outcome_name == "binom_of_num_sig_discoveries_p", "SO"]<=0.05, global_results_table_replication[global_results_table_replication$outcome_name == "stouffer_z", "SO"]>=1.64))
    
    global_results_replication["was_replication"] = 1
    if(global_results_first_research["DR_triggered"] == 1){global_results_replication["DR_triggered_in_previous_study"] = 1}
    if(global_results_first_research["HR_triggered"] == 1){global_results_replication["HR_triggered_in_previous_study"] = 1}
    
    if((global_results_replication["DR_triggered_in_previous_study"] == 1) & (global_results_replication["DR_triggered"] == 1)){global_results_replication["same_outcome_triggered_as_in_previous_study"] = 1}
    if((global_results_replication["HR_triggered_in_previous_study"] == 1) & (global_results_replication["HR_triggered"] == 1)){global_results_replication["same_outcome_triggered_as_in_previous_study"] = 1}
    
    if((global_results_replication["same_outcome_triggered_as_in_previous_study"] == 1) & (global_results_replication["number_of_effects_in_control_experiments"] < 2)){global_results_replication["replication_success"] = 1} else {global_results_replication["replication_success"] = 0}
    
    return(global_results_replication) # MUST BE global_results_replication
  }
}


#########################################
#           Run Simulation              #
#########################################

### Set parameters and run simulation

iter = 1000
sim_results = pbreplicate(iter, simulation_function(n_trials_first_research = 1408, 
                                                    n_trials_replication_research = 1408 * 1.5, 
                                                    n_labs = 7, 
                                                    prob_X = 0.28,
                                                    prob_O = 0.25,
                                                    prob_SX = 0.25,
                                                    prob_SO = 0.25))

sim_results_table_pre = data.frame(t(sim_results))

### power to detect the effect in the whole design
sum(sim_results_table_pre["replication_success"])/nrow(sim_results_table_pre)


iter = 1000
sim_results = pbreplicate(iter, simulation_function(n_trials_first_research = 1408, 
                                                    n_trials_replication_research = 1408 * 1.5, 
                                                    n_labs = 7, 
                                                    prob_X = 0.25,
                                                    prob_O = 0.25,
                                                    prob_SX = 0.25,
                                                    prob_SO = 0.25))

sim_results_table_pre = data.frame(t(sim_results))

### alpha, false positive rate in the whole design
sum(sim_results_table_pre["replication_success"])/nrow(sim_results_table_pre)

