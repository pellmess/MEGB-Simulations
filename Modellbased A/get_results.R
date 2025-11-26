library(checkmate)
library(tidyverse)
source("auxiliar/QualityMeasure.R")
get_estimations <- function(estimator = "bhf", mod_name){
  
  pfad <- paste0(mod_name, "/results/", estimator)
    
    csv_index <-
      list.files(pfad) %>% grep(pattern = "es.csv")

    csv_path <- list.files(pfad, full.names = T)[csv_index]
    if(length(csv_path) == 0){
      stop(paste("no files found for:", pfad))
    }
    if(length(csv_path) > 1) {
      tmp_csv <- lapply(csv_path, function(path){
        read.csv(path, sep = " ") 
      }) %>% bind_rows()
    } else {
      tmp_csv <- read.csv(csv_path, sep = " ")
    }
    assert_data_frame(tmp_csv, any.missing = F, min.rows = 1, min.cols = 1)

    tmp_results <- tmp_csv %>% arrange(run, domain) %>%
      dplyr::select(c("run", "Mean", "domain")) %>%
      reshape(idvar = "run",
              timevar = "domain",
              direction = "wide") %>%
      dplyr::select(-c("run")) %>%
      t() %>%
      `row.names<-`(NULL) %>%
      `colnames<-`(NULL)
    
    tmp_results
}







get_quality <- function(results, mod_name){
  
  domain.pop.mean <- read_csv(paste0(mod_name, "/results/pop.value.csv"))$mean
  

  n <- length(as.vector(t(domain.pop.mean))) / 500 
  

  vector <- as.vector(t(domain.pop.mean))
  

  transposed_matrix <- matrix(vector, nrow = n, ncol = 500)
  

  original_matrix <- t(transposed_matrix)
  

  lapply(results, function(x){
    QualityMeasure(
      True_mean = t(original_matrix), 
      Est_mean = x,
      MSETF = F
    )
  })
  
}


plot_estimations <- function(estimations, qual = c("Bias", "RB", "RRMSE", "True.RMSE")){
  
  lapply(qual, function(x) {
    bquote(ggplot(estimations, aes(x = estimator, y =  .(as.symbol(x)), fill = estimator)) +
      geom_boxplot(col = "gray") +
        ggtitle(x) +
        ylab(x) +
        scale_fill_manual(values = viridis::viridis_pal(begin = 0.2, end = 0.8)(11)[1:11]) +
        theme_bw(base_size = 20) +
        theme(
          axis.title.x = element_blank(),
          legend.position = "none",
          legend.title = element_blank(),
          axis.text.x = element_text(
            angle = 25,
            vjust = 0.5,
            hjust = 1,
            size = 8
          )
        )) |> eval()
    
  })
} 


report <- function(estimators, model_name, qual = c("Bias", "RB", "RRMSE", "True.RMSE")) {
  results <- setNames(lapply(estimators, function(x) {
    get_estimations(x, model_name)
  }), estimators)
  
  

  
  est_quality <- get_quality(results, model_name) %>%
    data.table::rbindlist(idcol = "estimator")
  
  # Ensure columns are numeric
  est_quality <- est_quality %>%
    mutate(across(c(Bias, RB, RRMSE, True.RMSE), as.numeric))
  
 
    
      res_tabelle <- est_quality %>%
        mutate(estimator = as.character(estimator)) %>%
        arrange(toupper(estimator)) %>%
        group_by(estimator) %>%
        dplyr::summarise(across(where(is.numeric), list(mean = mean, median = median), .names = "{.col}_{.fn}")) %>%
        mutate(across(where(is.numeric), round, 3)) %>%
        dplyr::select(estimator, Bias_mean, Bias_median, RB_mean, RB_median, RRMSE_mean, RRMSE_median, True.RMSE_mean, True.RMSE_median) %>%
        kableExtra::kbl(
          col.names = c("Estimator", "Mean", "Median", "Mean", "Median", "Mean", "Median", "Mean)", "Median"),
          format = "latex",
          booktabs = TRUE,
          digits = 3,
          linesep = ""
        ) %>%
        kableExtra::kable_classic(full_width = F, latex_options = c("striped", "hold_position")) %>%
        kableExtra::add_header_above(c(" " = 1, "Bias" = 2, "RB" = 2, "RRMSE" = 2, "True RMSE" = 2), 
                                     escape = FALSE, line = TRUE)
      tabelle_df <- est_quality %>%
        mutate(estimator = as.character(estimator)) %>%
        arrange(toupper(estimator)) %>%
        group_by(estimator) %>%
        dplyr::summarise(across(where(is.numeric), list(mean = mean, median = median), .names = "{.col}_{.fn}")) %>%
        mutate(across(where(is.numeric), ~ .x * 100)) %>%  
        mutate(across(where(is.numeric), round, 3)) %>%  
        dplyr::select(estimator, Bias_mean, Bias_median, RB_mean, RB_median, RRMSE_mean, RRMSE_median, True.RMSE_mean, True.RMSE_median)
      
      
      res <- plot_estimations(estimations = est_quality)
      list(plots = res, tabelle = res_tabelle, tabelle_df = tabelle_df)
}



color_mapping <- c(
  "BHF"          = "#440154FF",
  "EBP"       = "#3B528BFF",
  "GB"          = "#2A788EFF",
  "MEGB"         = "#21908CFF",
  "MERF"         = "#5DC863FF",
  "MBOOST"    = "#75D054FF",
  "MBOOST_T"  = "#8FD744FF",
  "MERTBoosting"    = "#ABDC31FF",
  "REEMBoosting"    = "#C8E020FF",
  "RF"           = "#E3E418FF",
  "XGBoost" = "#FDE725FF",
  "XGBoost-IDdummy"      = "#FFF079FF",
  "HT"           = "#ADD8E6"
)



estimators <- c(
  "BHF",
  "EBP",
  "MEGB",
  "MERF",
  "XGBoost",
  "XGBoost-IDdummy",
  "MERTBoosting",
  "REEMBoosting",
  "RF"
)


generate_results <- function(estimators, model_name_suffix = "") {
  # Define models
  mods <- c("mod1", "mod2", "mod3", "mod4")
  
  # Initialize lists to store reports, bias plots, and RMSE boxplots
  report_list <- list()
  bias_plots <- list()
  rmse_boxplots <- list()
  size <- 24
  # Generate reports, bias plots, and RMSE boxplots for each model
  for (i in mods) {
    # Generate report for the current model
    report_list[[i]] <- report(estimators = estimators, model_name = i)
    
    # Generate bias plot for the current model
    bias_plots[[i]] <- report_list[[i]]$plots[[2]] + 
      scale_fill_manual(values = color_mapping) +
      theme_bw() +
      xlab("") +
      theme(
        text = element_text(size = size),
        plot.title = element_text(size = size),
        axis.title = element_text(size = size),
        legend.text = element_text(size = size),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      geom_hline(
        yintercept = 0,
        linetype = "dashed",
        col = "red",
        size = 1.2
      )
    
    # Generate RMSE boxplot for the current model
    rmse_boxplots[[i]] <- report_list[[i]]$plots[[3]] + 
      scale_fill_manual(values = color_mapping) +
      theme_bw() +
      xlab("") +
      theme(
        text = element_text(size = size),
        plot.title = element_text(size = size),
        axis.title = element_text(size = size),
        legend.text = element_text(size = size),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    
    # Optionally save each combined plot
    combined_plot <- grid.arrange(bias_plots[[i]], rmse_boxplots[[i]], ncol = 2)
    ggsave(filename = paste0("Stat Woche/combined_mod", i, model_name_suffix, "_testskript_2024.jpg"), combined_plot, height = height, width = width)
  }
  
  
  arrange_bias <- ggarrange(
    bias_plots[[1]] + ggtitle("(a) Linear-Normal"),
    bias_plots[[2]] + ggtitle("(b) Complex-Normal"),
    bias_plots[[3]] + ggtitle("(c) Linear-Pareto"),
    bias_plots[[4]] + ggtitle("(d) Complex-Pareto"),
    ncol = 2,
    nrow = 2,
    common.legend = TRUE
  )
  
  
  arrange_boxrmse <- ggarrange(
    rmse_boxplots[[1]] + ggtitle("(a) Linear-Normal"),
    rmse_boxplots[[2]] + ggtitle("(b) Complex-Normal"),
    rmse_boxplots[[3]] + ggtitle("(c) Linear-Pareto"),
    rmse_boxplots[[4]] + ggtitle("(d) Complex-Pareto"),
    ncol = 2,
    nrow = 2,
    common.legend = TRUE
  )
  
  
  ggsave(paste0("plots/mb_boxplot_rb", model_name_suffix, ".jpg"), plot = arrange_bias, width = 16, height = 12)
  ggsave(paste0("plots/mb_boxplot_rmse", model_name_suffix, ".jpg"), plot = arrange_boxrmse, width = 16, height = 12)
}

generate_results(estimators = estimators)
