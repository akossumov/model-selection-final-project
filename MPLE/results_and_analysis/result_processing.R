############ CREATE TABLES OF COUNTS OF INFORMTION CRITERIA

library(dplyr)
library(xtable)
source("result_processing_functions.R")

# File path (edit as needed)
files_path <- "all_used_logs_1000_replications"
MODEL_SELECTION_output_filename <- "MODEL_SELECTION.tex"
print(files_path)

#tau_list <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.5, 0.75)
tau_list <- c(0.05, 0.1, 0.15, 0.2, 0.25)
n_list <- c(100, 200, 400)

cop_TYPES <-  c("clayton", "gumbel", "joe", "frank", "gaussian")
all_criteria <- c("AIC", "xv1", "xv_CIC", "xvK_2", "xvK_5", "xvK_10",  "xvK_n:2", "xvK_n:logn", "xvK_n05", "xvK_n05:logn")

# all outputs  are given in "count_output_filename" '
MODEL_SELECTION(files_path, n_list, tau_list, cop_TYPES, all_criteria, MODEL_SELECTION_output_filename)

############ COMPUTATION OF COINCIDENCE PERCENTAGES

num_repl <- 1000

(coinc_perc_AIC_xv1_STRING <- COINC_PERC(files_path, n_list, tau_list, cop_TYPES, all_criteria, "AIC", "xv1"))
(coinc_perc_AIC_xvK_n05_STRING <- COINC_PERC(files_path, n_list, tau_list, cop_TYPES, all_criteria, "AIC", "xvK_n05"))
(coinc_perc_AIC_xv_CIC_STRING <- COINC_PERC(files_path, n_list, tau_list, cop_TYPES, all_criteria, "AIC", "xv_CIC"))

#colnames(coinc_perc_AIC_xv1_STRING) <- c("$\\tau = 0.05$", "$\\tau =0.1$", "$\\tau=0.15$", "$\\tau=0.2$", "$\\tau=0.25$", "$\\tau=0.5$", "$\\tau=0.75$", "All" )
#colnames(coinc_perc_AIC_xvK_n05_STRING) <- c("$\\tau = 0.05$", "$\\tau =0.1$", "$\\tau=0.15$", "$\\tau=0.2$", "$\\tau=0.25$", "$\\tau=0.5$", "$\\tau=0.75$", "All" )
#colnames(coinc_perc_AIC_xv_CIC_STRING) <- c("$\\tau = 0.05$", "$\\tau =0.1$", "$\\tau=0.15$", "$\\tau=0.2$", "$\\tau=0.25$", "$\\tau=0.5$", "$\\tau=0.75$", "All" )

colnames(coinc_perc_AIC_xv1_STRING) <- c("$\\tau = 0.05$", "$\\tau =0.1$", "$\\tau=0.15$", "$\\tau=0.2$", "$\\tau=0.25$", "All")
colnames(coinc_perc_AIC_xvK_n05_STRING) <- c("$\\tau = 0.05$", "$\\tau =0.1$", "$\\tau=0.15$", "$\\tau=0.2$", "$\\tau=0.25$", "All")
colnames(coinc_perc_AIC_xv_CIC_STRING) <- c("$\\tau = 0.05$", "$\\tau =0.1$", "$\\tau=0.15$", "$\\tau=0.2$", "$\\tau=0.25$","All")


(coinc_perc_AIC_xv1_NUMERIC <- round(as.data.frame(
  apply(coinc_perc_AIC_xv1_STRING, 2, function(x) as.numeric(as.character(x)))
), digits = 2))

(coinc_perc_AIC_xvK_n05_NUMERIC <- round(as.data.frame(
  apply(coinc_perc_AIC_xvK_n05_STRING, 2, function(x) as.numeric(as.character(x)))
), digits = 2))

(coinc_perc_AIC_xv_CIC_NUMERIC <- round(as.data.frame(
  apply(coinc_perc_AIC_xv_CIC_STRING, 2, function(x) as.numeric(as.character(x)))
), digits = 2))

### standard errors of coincidence percentages

(se_coinc_perc_AIC_xv1_STRING <- SE_COINC_PERC(files_path, n_list, tau_list, cop_TYPES, all_criteria, "AIC", "xv1", num_repl))
#(se_coinc_perc_AIC_xvK_n05_STRING <- SE_COINC_PERC(files_path, n_list, tau_list, cop_TYPES, all_criteria, "AIC", "xvK_n05", num_repl))
(se_coinc_perc_AIC_xv_CIC_STRING <- SE_COINC_PERC(files_path, n_list, tau_list, cop_TYPES, all_criteria, "AIC", "xv_CIC", num_repl))

(se_coinc_perc_AIC_xv1_NUMERIC <- round(as.data.frame(
  apply(se_coinc_perc_AIC_xv1_STRING, 2, function(x) as.numeric(as.character(x)))
), digits = 3))

#(se_coinc_perc_AIC_xvK_n05_NUMERIC <- round(as.data.frame(
#  apply(se_coinc_perc_AIC_xvK_n05_STRING, 2, function(x) as.numeric(as.character(x)))
#), digits = 3))

(se_coinc_perc_AIC_xv_CIC_NUMERIC <- round(as.data.frame(
  apply(se_coinc_perc_AIC_xv_CIC_STRING, 2, function(x) as.numeric(as.character(x)))
), digits = 3))

caption_text_AIC_xv1 <- sprintf("Coincedence of $\\text{AIC}$ and $\\text{xv}_{1}$, with 95 \\%% confidence intervals (all values multiplied by 100).")
#caption_text_AIC_xvK_n05 <- sprintf("Coincedence of $\\text{AIC}$ and $\\text{xv}_{K_{n}}$ with $K_{n} = \\sqrt{n}$, with 95 \\%% confidence intervals (all values multiplied by 100).")
caption_text_AIC_xv_CIC <- sprintf("Coincedence of $\\text{AIC}$ and $\\text{xv}\\_\\text{CIC}$, with 95 \\%% confidence intervals (all values multiplied by 100).")


AIC_xv1 <- matrix(
  mapply(make_latex_entry, as.matrix(coinc_perc_AIC_xv1_NUMERIC), as.matrix(se_coinc_perc_AIC_xv1_NUMERIC)),
  nrow = nrow(as.matrix(coinc_perc_AIC_xv1_STRING)),
  ncol = ncol(as.matrix(se_coinc_perc_AIC_xv1_STRING)),
  byrow = FALSE,
  dimnames = dimnames(coinc_perc_AIC_xv1_STRING)
)

print(xtable(AIC_xv1, align = rep("c", ncol(AIC_xv1) + 1), caption = caption_text_AIC_xv1), 
      sanitize.text.function = identity, 
      caption.placement = "bottom",
      na.string = "")

#AIC_xvK_n05 <- matrix(
#  mapply(make_latex_entry, as.matrix(coinc_perc_AIC_xvK_n05_NUMERIC), as.matrix(se_coinc_perc_AIC_xvK_n05_NUMERIC)),
#  nrow = nrow(as.matrix(coinc_perc_AIC_xvK_n05_STRING)),
#  ncol = ncol(as.matrix(se_coinc_perc_AIC_xvK_n05_STRING)),
#  byrow = FALSE,
#  dimnames = dimnames(coinc_perc_AIC_xvK_n05_STRING)
#)

#print(xtable(AIC_xvK_n05, align = rep("c", ncol(AIC_xvK_n05) + 1), caption = caption_text_AIC_xvK_n05),
#      sanitize.text.function = identity, 
#      na.string = "")

AIC_xv_CIC <- matrix(
  mapply(make_latex_entry, as.matrix(coinc_perc_AIC_xv_CIC_NUMERIC), as.matrix(se_coinc_perc_AIC_xv_CIC_NUMERIC)),
  nrow = nrow(as.matrix(coinc_perc_AIC_xv_CIC_STRING)),
  ncol = ncol(as.matrix(se_coinc_perc_AIC_xv_CIC_STRING)),
  byrow = FALSE,
  dimnames = dimnames(coinc_perc_AIC_xv_CIC_STRING)
)

print(xtable(AIC_xv_CIC, align = rep("c", ncol(AIC_xv_CIC) + 1), caption = caption_text_AIC_xv_CIC), 
      sanitize.text.function = identity, 
      na.string = "")

############ HIT RATES

HIT_RATES_output_filename <- "HIT_RATES.tex"
all_criteria <- c("AIC", "xv1", "xv_CIC", "xvK_2", "xvK_5", "xvK_10",  "xvK_n.2", "xvK_n.logn", "xvK_n05", "xvK_n05.logn")
HIT_RATES(files_path, n_list, tau_list, cop_TYPES, all_criteria, HIT_RATES_output_filename)

