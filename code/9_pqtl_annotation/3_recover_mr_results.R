##############
#INTRODUCTION#
##############

#This code recovers MR results and produces a table with them

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)
library(TwoSampleMR)

###################
#Loading functions#
###################

stats_compiler <- function(dir_, out_check){
  
  #STEP 1: do a conditional and assess:
  
  if(out_check == "outliers_removed"){
    
    mr_stats <- readRDS(paste(dir_, "/mr_res_after_outlier_extraction.RDS", sep = ""))
    
  } else {
    
    mr_stats <- readRDS(paste(dir_, "/mr_res_before_outlier_extraction.RDS", sep = ""))
    
  }
  
  return(mr_stats)
  
}

sensitivity_compiler <- function(dir_, out_check){
  
  #STEP 1: do a conditional and assess:
  
  if(out_check == "outliers_removed"){
    
    rucker <- readRDS(paste(dir_, "/sensitivity_test_after_outlier_extraction.RDS", sep = ""))
    het <- readRDS(paste(dir_, "/strict_het_test_after_outlier_extraction.RDS", sep = ""))
    
  } else {
    
    rucker <- readRDS(paste(dir_, "/sensitivity_test_before_outlier_extraction.RDS", sep = ""))
    het <- readRDS(paste(dir_, "/strict_het_test_before_outlier_extraction.RDS", sep = ""))
    
  }
  
  return(list(rucker, het))
  
}

# compiling_res <- function(list_of_dirs){
#   
#   #STEP 0:initiate loop: 
#   
#   for(index in seq(1, length(list_of_dirs))){
#     
#     #STEP 0: get the protein:
#     
#     prot <- as.character(unlist(str_split(list_of_dirs[index], "/")))
#     prot <- prot[length(prot)]
#     print(prot)
#     
#     #STEP 1: check data
#     
#     list_of_files <- list.files(list_of_dirs[index])
#     
#     analysis_check <- ifelse("mr_res_before_outlier_extraction.RDS"%in%list_of_files, "Performed", "other")
#     outlier_check <- ifelse("mr_res_after_outlier_extraction.RDS"%in%list_of_files, "outliers_removed", "standard")
#     nsnps_check <- ifelse("sensitivity_test_before_outlier_extraction.RDS"%in%list_of_files, analysis_check, "nSNPs < 3")
#     
#     if(nsnps_check == "nSNPs < 3"){
#       
#       next()
#       
#     }
#     
#     #STEP 2: get the stats!!!!
#     
#     if(analysis_check == "Performed"){
#       
#       #CAREFUL, we will start with an exception. 
#       #Maybe the analysis ran with 3 SNPs, but outlier removal took 1 out.
#       #We will allow it to pass, but we will check the sensitivity tests...
#       
#       nsnps_check <- ifelse("sensitivity_test_after_outlier_extraction.RDS"%in%list_of_files, analysis_check, "nSNPs < 3")
#       
#       if(nsnps_check == "nSNPs < 3"){
#         
#         outlier_check <- "standard"
#         
#         mr_stats <- stats_compiler(list_of_dirs[index], outlier_check)
#         mr_stats$protein <- prot
#         mr_stats$analyses <- outlier_check
#         
#         
#       }
#       
#       #If everything is normal, we are cool:
#       
#       mr_stats <- stats_compiler(list_of_dirs[index], outlier_check)
#       mr_stats$protein <- prot
#       mr_stats$analyses <- outlier_check
#       
#       
#     } else {
#       
#       mr_stats <- c("id.exposure", "id.outcome", "outcome", "exposure", "method", "nsnp", "b", "se", "pval")
#       mr_stats <- as.data.frame(rbind(mr_stats, rep(NA, 9)))
#       colnames(mr_stats) <- mr_stats[1,]
#       mr_stats <- mr_stats[-1,]
#       
#       mr_stats$protein <- prot
#       mr_stats$analyses <- analysis_check
#       
#     }
#     
#     #STEP 3: let's check now the sensitivity tests:
#     
#     sensitivity_list <- sensitivity_compiler(list_of_dirs[index], outlier_check)
#     
#     cochran_q_diff <- paste(sensitivity_list[[1]][[1]]$Q$Q[3], " P=", sensitivity_list[[1]][[1]]$Q$P[3], sep = "") #gets in the rucker list, extracts Q info for the Q_diff [3] between IVW and Egger!
#     intercept <- paste(sensitivity_list[[1]][[1]]$intercept$Estimate[1], " P=", sensitivity_list[[1]][[1]]$intercept$P[1], sep = "")
#     isq <- unlist(sensitivity_list[[2]])
#     
#     if(length(isq) != 1){ #this only happens when it is empty
#       
#       isq <- paste(isq$I2, "[", isq$lower.I2, ",", isq$upper.I2, "]", sep = "")
#       
#     }
#     
#     mr_stats$Q_diff <- cochran_q_diff
#     mr_stats$egger_intercept <- intercept
#     mr_stats$I2 <- isq
#     
#     #Finally, let's update this!!!
#     
#     if(!(exists("final_df"))){
#       
#       final_df <- mr_stats
#       
#     } else {
#       
#       final_df <- rbind(final_df, mr_stats)
#       
#     }
#     
#   }
#   
#   return(final_df)
#   
# }

compiling_res <- function(list_of_dirs){
  
  #STEP 0:initiate loop: 
  
  for(index in seq(1, length(list_of_dirs))){
    
    #STEP 0: get the outcome:
    
    outcome <- as.character(unlist(str_split(list_of_dirs[index], "/")))
    outcome <- outcome[length(outcome)]
    print(outcome)
    
    #STEP 1: check data
    
    list_of_files <- list.files(list_of_dirs[index])
    
    analysis_check <- ifelse("mr_res_before_outlier_extraction.RDS"%in%list_of_files, "Performed", "other")
    outlier_check <- ifelse("mr_res_after_outlier_extraction.RDS"%in%list_of_files, "outliers_removed", "standard")
    nsnps_check <- ifelse("sensitivity_test_before_outlier_extraction.RDS"%in%list_of_files, analysis_check, "nSNPs < 3")
    
    if(nsnps_check == "nSNPs < 3"){
      
      outlier_check <- "standard"
      
      mr_stats <- stats_compiler(list_of_dirs[index], outlier_check)
      mr_stats$outcome <- outcome
      mr_stats$analyses <- outlier_check
      
      mr_stats$Q_diff <- NA
      mr_stats$egger_intercept <- NA
      mr_stats$I2 <- NA
      
      if(!(exists("final_df"))){
        
        final_df <- mr_stats
        
      } else {
        
        final_df <- rbind(final_df, mr_stats)
        
      }
      
      next()
      
    }
    
    #STEP 2: get the stats!!!!
    
    if(analysis_check == "Performed"){
      
      #CAREFUL, we will start with an exception. 
      #Maybe the analysis ran with 3 SNPs, but outlier removal took 1 out.
      #We will allow it to pass, but we will check the sensitivity tests...
      
      nsnps_check <- ifelse("sensitivity_test_after_outlier_extraction.RDS"%in%list_of_files, analysis_check, "nSNPs < 3")
      
      if(nsnps_check == "nSNPs < 3"){
        
        outlier_check <- "standard"
        
        mr_stats <- stats_compiler(list_of_dirs[index], outlier_check)
        mr_stats$outcome <- outcome
        mr_stats$analyses <- outlier_check
        
        
      }
      
      #If everything is normal, we are cool:
      
      mr_stats <- stats_compiler(list_of_dirs[index], outlier_check)
      mr_stats$outcome <- outcome
      mr_stats$analyses <- outlier_check
      
      #Let's add the data if the sensitivity analyses were performed:
      
      sensitivity_list <- sensitivity_compiler(list_of_dirs[index], outlier_check)
      
      q_ <- format(sensitivity_list[[1]][[1]]$Q$Q[3], scientific = FALSE, digits = 2)
      q_p <- format(sensitivity_list[[1]][[1]]$Q$P[3], scientific = TRUE, digits = 3)
      cochran_q_diff <- paste(q_, " P=", q_p, sep = "") 
      
      #cochran_q_diff <- paste(sensitivity_list[[1]][[1]]$Q$Q[3], " P=", sensitivity_list[[1]][[1]]$Q$P[3], sep = "") #gets in the rucker list, extracts Q info for the Q_diff [3] between IVW and Egger!
      
      int_ <- format(sensitivity_list[[1]][[1]]$intercept$Estimate[1], scientific = TRUE, digits = 3)
      int_p <- format(sensitivity_list[[1]][[1]]$intercept$P[1], scientific = TRUE, digits = 3)
      
      #intercept <- paste(sensitivity_list[[1]][[1]]$intercept$Estimate[1], " P=", sensitivity_list[[1]][[1]]$intercept$P[1], sep = "")
      intercept <- paste(int_, " P=", int_p, sep = "") 
      
      
      isq <- unlist(sensitivity_list[[2]])
      
      if(length(isq) != 1){ #this only happens when it is empty
        
        isq <- paste(format(isq$I2, scientifiC=FALSE, digits=2), 
                     " [", format(isq$lower.I2, scientifiC=FALSE, digits=3), 
                     ",", format(isq$upper.I2, scientifiC=FALSE, digits=2), "]", sep = "")
        
      }
      
      mr_stats$Q_diff <- cochran_q_diff
      mr_stats$egger_intercept <- intercept
      mr_stats$I2 <- isq
      
      
    } else {
      
      mr_stats <- c("id.exposure", "id.outcome", "outcome", "exposure", "method", "nsnp", "b", "se", "pval")
      mr_stats <- as.data.frame(rbind(mr_stats, rep(NA, 9)))
      colnames(mr_stats) <- mr_stats[1,]
      mr_stats <- mr_stats[-1,]
      
      mr_stats$outcome <- outcome
      mr_stats$analyses <- analysis_check
      
      mr_stats$Q_diff <- NA
      mr_stats$egger_intercept <- NA
      mr_stats$I2 <- NA
      
    }
    
    #STEP 3: let's check now the sensitivity tests:
    
    if(!(exists("final_df"))){
      
      final_df <- mr_stats
      
    } else {
      
      final_df <- rbind(final_df, mr_stats)
      
    }
    
  }
  
  return(final_df)
  
}


protein_parser <- function(protein_id){
  
  #STEP 1: figure out what kind is it:
  
  check <- str_detect(protein_id, "UKB_PPP")
  
  if(check){
    
    prot <- unlist(str_split(protein_id, "UKB_PPP_EUR_"))[2]
    prot <- unlist(str_split(prot, "_"))[1]
    
    return(prot)
    
  } else {
    
    prot <- unlist(str_split(protein_id, "sun_2018_aptamer_plasma_"))[2]
    prot <- unlist(str_split(prot, "_"))[1]
    
    return(toupper(prot))
    
  }
  
}


lower_ci_simple <- function(beta_, se_){
  #A function that takes the beta and se from summary statistics
  #and generates the lower CI.
  
  lower_ci <- beta_ - qnorm(0.975)*se_
  
  return(lower_ci)
  
}

upper_ci_simple <- function(beta_, se_){
  #A function that takes the beta and se from summary statistics
  #and generates the upper CI.
  
  upper_ci <- beta_ + qnorm(0.975)*se_
  
  return(upper_ci)
  
}

# lower_ci <- function(beta_, se_, nsnp, method){
#   #This code makes confidence intervals according to the degrees of freedome for each of our analyses
#   
#   #STEP 0: make dummy variables to see that this works:
#   
#   # beta_=expected_res$b
#   # se_=expected_res$se
#   # nsnp=expected_res$nsnp
#   # method=expected_res$method
#   
#   #STEP 1: let's make a loop out of this:
#   
#   lower_ci_end <- c()
#   
#   for(index in seq(1, length(beta_))){
#     
#     print(index)
#     
#     #STEP 1.1: take degrees of freedom
#     
#     df <- switch(method[index],
#                  "MR Egger"           = nsnp[index] - 2,
#                  "Inverse variance weighted"             = nsnp[index] - 1,
#                  "Weighted median" = nsnp[index] - 1,
#                  "Weighted mode"   = nsnp[index] - 1,
#                  stop("Unknown method"))
#     
#     #STEP 1.2: compute lower ci
#     
#     lower_ci_ <- beta_[index] - qt(0.975, df) * se_[index]
#     
#     #STEP 1.3: return vector full of ci:
#     
#     lower_ci_end <- c(lower_ci_end, lower_ci_)
#     
#   }
#   
#   return(lower_ci_end)
# }
# 
# upper_ci <- function(beta_, se_, nsnp, method){
#   #This code makes confidence intervals according to the degrees of freedome for each of our analyses
#   
#   #STEP 0: make dummy variables to see that this works:
#   
#   #beta_=expected_res$b
#   #se_=expected_res$se
#   #nsnp=expected_res$nsnp
#   #method=expected_res$method
#   
#   #STEP 1: let's make a loop out of this:
#   
#   upper_ci_end <- c()
#   
#   for(index in seq(1, length(beta_))){
#     
#     #STEP 1.1: take degrees of freedom
#     
#     df <- switch(method[index],
#                  "MR Egger"           = nsnp[index] - 2,
#                  "Inverse variance weighted"             = nsnp[index] - 1,
#                  "Weighted median" = nsnp[index] - 1,
#                  "Weighted mode"   = nsnp[index] - 1,
#                  stop("Unknown method"))
#     
#     #STEP 1.2: compute lower ci
#     
#     upper_ci_ <- beta_[index] + qt(0.975, df) * se_[index]
#     
#     #STEP 1.3: return vector full of ci:
#     
#     upper_ci_end <- c(upper_ci_end, upper_ci_)
#     
#   }
#   
#   return(upper_ci_end)
# }


lower_ci <- function(beta_, se_, nsnp, method){
  lower_ci_end <- c()
  
  for(index in seq_along(beta_)){
    dist <- method[index]
    
    # Choose critical value based on method
    crit_val <- switch(dist,
                       "MR Egger"           = qt(0.975, df = nsnp[index] - 2),
                       "Weighted mode"      = qt(0.975, df = nsnp[index] - 1),
                       "Inverse variance weighted" = qnorm(0.975),
                       "Weighted median"    = qnorm(0.975),
                       stop(paste("Unknown method:", dist)))
    
    lower_ci_ <- beta_[index] - crit_val * se_[index]
    lower_ci_end <- c(lower_ci_end, lower_ci_)
  }
  
  return(lower_ci_end)
}

upper_ci <- function(beta_, se_, nsnp, method){
  upper_ci_end <- c()
  
  for(index in seq_along(beta_)){
    dist <- method[index]
    
    crit_val <- switch(dist,
                       "MR Egger"           = qt(0.975, df = nsnp[index] - 2),
                       "Weighted mode"      = qt(0.975, df = nsnp[index] - 1),
                       "Inverse variance weighted" = qnorm(0.975),
                       "Weighted median"    = qnorm(0.975),
                       stop(paste("Unknown method:", dist)))
    
    upper_ci_ <- beta_[index] + crit_val * se_[index]
    upper_ci_end <- c(upper_ci_end, upper_ci_)
  }
  
  return(upper_ci_end)
}


mr_plotter_first <- function(mr_df, title){
  #This plot takes the dataframe of result from MR and tries to plot it so that it looks amazing and beautiful
  
  #STEP 1: make dummy data
  
  betas_rounded <- format(mr_df$b, scientific = FALSE, digits = 4)
  lower_ci_rounded <- format(mr_df$lower_ci,scientific = FALSE, digits = 4)
  upper_ci_rounded <- format(mr_df$upper_ci, scientific = FALSE, digits = 4)
  
  mr_df$betas <- betas_rounded
  mr_df$lower_ci <- lower_ci_rounded
  mr_df$upper_ci <- upper_ci_rounded
  
  #STEP 2: change easy_for_code handles to beauty_informative_for_figure:
  
  mr_df$outcome <- ifelse(mr_df$outcome == "t2d", "T2D (FinnGen)", mr_df$outcome)
  mr_df$outcome <- ifelse(mr_df$outcome == "isiadjbmi", "ISIadjBMI (MAGIC)", mr_df$outcome)
  mr_df$outcome <- ifelse(mr_df$outcome == "fiadjbmi", "FIadjBMI (MAGIC)", mr_df$outcome)
  #mr_df$outcome <- ifelse(mr_df$outcome == "meta_analysis", "PCOS (Day et al)", mr_df$outcome)
  #mr_df$outcome <- ifelse(mr_df$outcome == "age_adjusted", "PCOS (Tyrmi et al)", mr_df$outcome)
  #mr_df$outcome <- ifelse(mr_df$outcome == "age_and_bmi_adjusted", "PCOSadjBMI (Tyrmi et al)", mr_df$outcome)
  
  #Let's convert the data into factors so that we can organize the results:
  
  mr_df$outcome <- factor(mr_df$outcome, levels = c("FIadjBMI (MAGIC)", "ISIadjBMI (MAGIC)", "T2D (FinnGen)"))
  mr_df$method  <- factor(mr_df$method, levels = c("Weighted mode", "Weighted median", "MR Egger", "Inverse variance weighted"))
  
  # Ensure pvals numeric and format
  mr_df$pvals_num <- as.numeric(mr_df$pval)
  mr_df$pval_formatted <- paste0("p=", format(mr_df$pvals_num, scientific = TRUE, digits = 2))
  
  # Create significance stars column based on p-value thresholds
  mr_df$signif_stars <- ""
  mr_df$signif_stars[mr_df$pvals_num < 0.05] <- "*"
  mr_df$signif_stars[mr_df$pvals_num < 5e-8] <- "**"
  
  # Position to place stars slightly above upper_ci
  star_y_pos <- as.numeric(mr_df$upper_ci) + 0.01  # tweak this offset if needed
  
  #Now we can make our plot!
  
  plotio <- ggplot(mr_df, aes(x = fct_rev(outcome), y = as.numeric(betas))) + #make a plot that takes into account outcome and betas
    
    geom_point(aes(color=method, shape=method), 
               size = 5, position = position_dodge(width = 0.75)) + #color by method, make them a bit separated to distinguish
    
    #A little change so that we can have a more beatiful plot:
    
    scale_shape_manual(values = c(
      "Inverse variance weighted" = 16,
      "Weighted median" = 16,
      "Weighted mode" = 16,
      "MR Egger" = 16
    )) +
    
    geom_errorbar(
      aes(
        ymin = as.numeric(lower_ci),
        ymax = as.numeric(upper_ci),
        group = method  # this is what was missing
      ),
      width = 0.15,
      position = position_dodge(width = 0.75),
      color = "black",
      linewidth = 1
    ) +
    
    geom_text(aes(y = star_y_pos, label = signif_stars, group = method),
              position = position_dodge(width = 0.75),
              size = 6, color = "black", fontface = "bold")+
    
    #facet_wrap(~outcome, strip.position = "left", nrow = 7, scales = "free_y") + #separate by outcome
    
    geom_hline(yintercept = 0, color = "red", size = 0.5, linetype = 2) +
    
    coord_flip() +
    xlab("Trait") +
    ylab("MR causal effect sizes") +
    ggtitle(title) +
    theme_bw() +
    theme(
      legend.position = "right",
      legend.text = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(),
      axis.text = element_text(size = 11),
      axis.title.x =  element_text(size = 11.5),
      axis.title.y =  element_blank(),
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5)
    ) +
    guides(
      color = guide_legend(reverse = TRUE),
      shape = guide_legend(reverse = TRUE)
    )
  
  plotio <- plotio + scale_y_continuous(limits = c(-0.10, 0.10), breaks = seq(as.numeric(format(-0.10, digits = 3)), as.numeric(format(0.10, digits = 3)), as.numeric(format(0.05, digits = 3)))) #for some reason I have the axis number showing me weird stuf
  
  return(plotio)
  
}

mr_plotter_last <- function(mr_df, title){
  #This plot takes the dataframe of result from MR and tries to plot it so that it looks amazing and beautiful
  
  #STEP 1: make dummy data
  
  betas_rounded <- format(mr_df$b, scientific = FALSE, digits = 4)
  lower_ci_rounded <- format(mr_df$lower_ci,scientific = FALSE, digits = 4)
  upper_ci_rounded <- format(mr_df$upper_ci, scientific = FALSE, digits = 4)
  
  mr_df$betas <- betas_rounded
  mr_df$lower_ci <- lower_ci_rounded
  mr_df$upper_ci <- upper_ci_rounded
  
  #STEP 2: change easy_for_code handles to beauty_informative_for_figure:
  
  mr_df$outcome <- ifelse(mr_df$outcome == "t2d", "T2D (FinnGen)", mr_df$outcome)
  mr_df$outcome <- ifelse(mr_df$outcome == "isiadjbmi", "ISIadjBMI (MAGIC)", mr_df$outcome)
  mr_df$outcome <- ifelse(mr_df$outcome == "fiadjbmi", "FIadjBMI (MAGIC)", mr_df$outcome)
  #mr_df$outcome <- ifelse(mr_df$outcome == "meta_analysis", "PCOS (Day et al)", mr_df$outcome)
  #mr_df$outcome <- ifelse(mr_df$outcome == "age_adjusted", "PCOS (Tyrmi et al)", mr_df$outcome)
  #mr_df$outcome <- ifelse(mr_df$outcome == "age_and_bmi_adjusted", "PCOSadjBMI (Tyrmi et al)", mr_df$outcome)
  
  #Let's convert the data into factors so that we can organize the results:
  
  mr_df$outcome <- factor(mr_df$outcome, levels = c("FIadjBMI (MAGIC)", "ISIadjBMI (MAGIC)", "T2D (FinnGen)"))
  mr_df$method  <- factor(mr_df$method, levels = c("Weighted mode", "Weighted median", "MR Egger", "Inverse variance weighted"))
  
  # Ensure pvals numeric and format
  mr_df$pvals_num <- as.numeric(mr_df$pval)
  mr_df$pval_formatted <- paste0("p=", format(mr_df$pvals_num, scientific = TRUE, digits = 2))
  
  # Create significance stars column based on p-value thresholds
  mr_df$signif_stars <- ""
  mr_df$signif_stars[mr_df$pvals_num < 0.05] <- "*"
  mr_df$signif_stars[mr_df$pvals_num < 5e-8] <- "**"
  
  # Position to place stars slightly above upper_ci
  star_y_pos <- as.numeric(mr_df$upper_ci) + 0.01  # tweak this offset if needed
  
  #Now we can make our plot!
  
  plotio <- ggplot(mr_df, aes(x = fct_rev(outcome), y = as.numeric(betas))) + #make a plot that takes into account outcome and betas
    
    geom_point(aes(color=method, shape=method), 
               size = 5, position = position_dodge(width = 0.75)) + #color by method, make them a bit separated to distinguish
    
    #A little change so that we can have a more beatiful plot:
    
    scale_shape_manual(values = c(
      "Inverse variance weighted" = 16,
      "Weighted median" = 16,
      "Weighted mode" = 16,
      "MR Egger" = 16
    )) +
    
    geom_errorbar(
      aes(
        ymin = as.numeric(lower_ci),
        ymax = as.numeric(upper_ci),
        group = method  # this is what was missing
      ),
      width = 0.15,
      position = position_dodge(width = 0.75),
      color = "black",
      linewidth = 1
    ) +
    
    geom_text(aes(y = star_y_pos, label = signif_stars, group = method),
              position = position_dodge(width = 0.75),
              size = 6, color = "black", fontface = "bold")+
    
    #facet_wrap(~outcome, strip.position = "left", nrow = 7, scales = "free_y") + #separate by outcome
    
    geom_hline(yintercept = 0, color = "red", size = 0.5, linetype = 2) +
    
    coord_flip() +
    xlab("Trait") +
    ylab("MR causal effect sizes") +
    ggtitle(title) +
    theme_bw() +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 11),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(),
      axis.text = element_text(size = 11),
      axis.title.x =  element_text(size = 11.5),
      axis.title.y =  element_blank(),
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5)
    ) +
    guides(
      color = guide_legend(reverse = TRUE),
      shape = guide_legend(reverse = TRUE)
    )
  
  plotio <- plotio + scale_y_continuous(limits = c(-0.10, 0.10), breaks = seq(as.numeric(format(-0.10, digits = 3)), as.numeric(format(0.10, digits = 3)), as.numeric(format(0.05, digits = 3)))) #for some reason I have the axis number showing me weird stuf
  
  return(plotio)
  
}



compute_p_value <- function(z) {
  return(2 * (1 - pnorm(abs(z))))
}


##############
#Loading data#
##############

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/")

fiadjbmi_path_2_res <- "output/9_pqtl_results/2_mr/1_mr_res/fiadjbmi_18032026/"
#isiadjbmi_path_2_res <- "output/9_pqtl_results/2_mr/1_mr_res/isiadjbmi/"
#t2d_path_2_res <- "output/9_pqtl_results/2_mr/1_mr_res/t2d/"

# bmi_path_2_res <- "output/9_pqtl_results/2_mr/1_mr_res/bmi/"
# whradjbmi_path_2_res <- "output/9_pqtl_results/2_mr/1_mr_res/whradjbmi/"
# hcadjbmi_path_2_res <- "output/9_pqtl_results/2_mr/1_mr_res/hcadjbmi/"

##############################################################################################################################
#The approach will be making a function that loops over the proteins in each folder that we give and searches for the results#
##############################################################################################################################

#We have been doing lots of analyses, let's focus on those that we finally ran: 

prot_og_df <- fread("output/9_pqtl_results/2_mr/0_data_4_input/all_associations_per_protein_w_more_than_3_credible_sets_18032026.txt") #131
prot_og_df$protein <- as.character(unlist(sapply(prot_og_df$studyId, protein_parser)))
prot_og_df$chr_pos <- paste("chr", prot_og_df$chr_hg19, ":", prot_og_df$pos_hg19, sep = "") #in build 37 to match associations from original GWAS for each protein

prot_list <- prot_og_df$prot_cred_set
prot_list <- prot_list[which(str_detect(prot_list, "UK"))]
prot_list <- as.character(unlist(sapply(prot_list, protein_parser)))
prot_list <- unique(prot_list) #20! Amazing
prot_list <- prot_list[order(prot_list)]

fiadjbmi_dirs <- paste("output/9_pqtl_results/2_mr/1_mr_res/fiadjbmi_18032026/" , prot_list, sep = "")
fiadjbmi_res <- compiling_res(fiadjbmi_dirs)

#Let's save this data because we are gonna use it as supplementary table 8:

fwrite(fiadjbmi_res, "manuscript/results/review_3/supplementary_tables/drafts/supplementary_table_19.txt")

# #######################################
# #Let's plot the results for KLK1 first#
# #######################################
# 
# klk1_fiadjbmi <- fiadjbmi_res[which(fiadjbmi_res$protein == "KLK1"),]
# #klk1_isiadjbmi <- isiadjbmi_res[which(isiadjbmi_res$protein == "KLK1"),]
# #klk1_t2d <- t2d_res[which(t2d_res$protein == "KLK1"),]
# 
# #Let's add the analyses they come from:
# 
# klk1_fiadjbmi$outcome <- "fiadjbmi"
# #klk1_isiadjbmi$outcome <- "isiadjbmi"
# #klk1_t2d$outcome <- "t2d"
# 
# #Let's combine this data:
# 
# klk1_df <- klk1_fiadjbmi
# #Let's remove simple method:
# 
# klk1_df <- klk1_df[which(klk1_df$method != "Simple mode"),]
# 
# #Let's add the confidence intervals:
# 
# klk1_df$lower_ci <- lower_ci(klk1_df$b, klk1_df$se, klk1_df$nsnp, klk1_df$method)
# klk1_df$upper_ci <- upper_ci(klk1_df$b, klk1_df$se, klk1_df$nsnp, klk1_df$method)
# 
# betas_rounded <- format(klk1_df$b, scientific = FALSE, digits = 4)
# lower_ci_rounded <- format(klk1_df$lower_ci,scientific = FALSE, digits = 4)
# upper_ci_rounded <- format(klk1_df$upper_ci, scientific = FALSE, digits = 4)
# 
# klk1_df$full_effect <- paste(betas_rounded,
#                                    " (",
#                                    lower_ci_rounded,
#                                    ",",
#                                    upper_ci_rounded,
#                                    ")",
#                                    sep=""
# )
# 
# #Can we load this in mr_plotter?
# 
# klk1_plot <- mr_plotter_first(klk1_df, "KLK1")
# 
# ########################################
# #Let's plot the results for GIMAP7 next#
# ########################################
# 
# gimap7_fiadjbmi <- fiadjbmi_res[which(fiadjbmi_res$protein == "GIMAP7"),]
# gimap7_isiadjbmi <- isiadjbmi_res[which(isiadjbmi_res$protein == "GIMAP7"),]
# gimap7_t2d <- t2d_res[which(t2d_res$protein == "GIMAP7"),]
# 
# #Let's add the analyses they come from:
# 
# gimap7_fiadjbmi$outcome <- "fiadjbmi"
# gimap7_isiadjbmi$outcome <- "isiadjbmi"
# gimap7_t2d$outcome <- "t2d"
# 
# #Let's combine this data:
# 
# gimap7_df <- rbind(gimap7_fiadjbmi, gimap7_isiadjbmi, gimap7_t2d)
# #Let's remove simple method:
# 
# gimap7_df <- gimap7_df[which(gimap7_df$method != "Simple mode"),]
# 
# #Let's add the confidence intervals:
# 
# gimap7_df$lower_ci <- lower_ci(gimap7_df$b, gimap7_df$se, gimap7_df$nsnp, gimap7_df$method)
# gimap7_df$upper_ci <- upper_ci(gimap7_df$b, gimap7_df$se, gimap7_df$nsnp, gimap7_df$method)
# 
# betas_rounded <- format(gimap7_df$b, scientific = FALSE, digits = 4)
# lower_ci_rounded <- format(gimap7_df$lower_ci,scientific = FALSE, digits = 4)
# upper_ci_rounded <- format(gimap7_df$upper_ci, scientific = FALSE, digits = 4)
# 
# gimap7_df$full_effect <- paste(betas_rounded,
#                              " (",
#                              lower_ci_rounded,
#                              ",",
#                              upper_ci_rounded,
#                              ")",
#                              sep=""
# )
# 
# #Can we load this in mr_plotter?
# 
# gimap7_plot <- mr_plotter_first(gimap7_df, "GIMAP7")
# 
# 
# # ################################################################################
# # #Let's retrieve the ones that are significant for, at least, two of the methods#
# # ################################################################################
# # 
# # sign_fiadjbmi <- fiadjbmi_res[which(as.numeric(fiadjbmi_res$pval) < 0.05),]
# # sign_isiadjbmi <- isiadjbmi_res[which(as.numeric(isiadjbmi_res$pval) < 0.05),]
# # sign_t2d <- t2d_res[which(as.numeric(t2d_res$pval) < 0.05),]
# # 
# # sign_fiadjbmi <- unique(sign_fiadjbmi$protein[which(duplicated(sign_fiadjbmi$protein))]) #GIMAP7, KLK1
# # sign_isiadjbmi <- unique(sign_isiadjbmi$protein[which(duplicated(sign_isiadjbmi$protein))]) #COL18A1, MEGF9
# # sign_t2d <- unique(sign_t2d$protein[which(duplicated(sign_t2d$protein))]) #KLK1, LPL, MEGF9 and NECTIN2
# # 
# # ###############################################################
# # #Let's check the SNPs for our the ones that replicate the most#
# # ###############################################################
# # 
# # pqtl_matched <- fread("output/9_pqtl_results/1_cis_pqtl_matches/cis_pqtl_matches.txt") #131
# # 
# # ###################################################
# # #Let's do the same for BMI, WHRadjBMI and HCadjBMI#
# # ###################################################
# # 
# # bmi_dirs <- paste("output/9_pqtl_results/2_mr/1_mr_res/bmi/" , prot_list, sep = "")
# # bmi_res <- compiling_res(bmi_dirs)
# # 
# # whradjbmi_dirs <- paste("output/9_pqtl_results/2_mr/1_mr_res/whradjbmi/" , prot_list, sep = "")
# # whradjbmi_res <- compiling_res(whradjbmi_dirs)
# # 
# # hcadjbmi_dirs <- paste("output/9_pqtl_results/2_mr/1_mr_res/hcadjbmi/" , prot_list, sep = "")
# # hcadjbmi_res <- compiling_res(hcadjbmi_dirs)
