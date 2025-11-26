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
  "XGBOOST"      = "#FDE725FF",
  "HT"           = "#ADD8E6"
)


library(pacman)
p_load(pbapply)
library(tidyverse)
res_files <- list.files("results/", pattern = "rds", full.names = T)

res_sae <- pblapply(res_files, function(x){
  readRDS(x)
})


keys <- sapply(res_sae, function(x) {
  s <- x$sim_obj$setting
  paste0("ICC=", s$ICC, "_n_domains=", s$n_domains, "_ni_domain=", s$ni_domain)
})


res_grouped <- split(res_sae, keys)


names(res_grouped)   
lengths(res_grouped)


setting_name <- names(res_grouped)[1]
qm_results <- lapply(names(res_grouped), function(setting_name){
  area_est_L2 <- lapply(res_grouped[[setting_name]], function(x) x$area_est_L2$mw) %>% do.call(cbind, .)
  area_est_BoostMERT_L2 <- lapply(res_grouped[[setting_name]], function(x) x$area_est_BoostMERT_L2$mw) %>% do.call(cbind, .)
  area_est_BoostREEM_L2 <- lapply(res_grouped[[setting_name]], function(x) x$area_est_BoostREEM_L2$mw) %>% do.call(cbind, .)
  area_est_MBOOST_L2 <- lapply(res_grouped[[setting_name]], function(x) x$area_est_MBOOST_L2$mw) %>% do.call(cbind, .)
  area_est_MBOOST_T_L2 <- lapply(res_grouped[[setting_name]], function(x) x$area_est_MBOOST_T_L2$mw) %>% do.call(cbind, .)
  area_est_RF <- lapply(res_grouped[[setting_name]], function(x) x$area_est_RF$mw) %>% do.call(cbind, .)
  area_est_MEGB <- lapply(res_grouped[[setting_name]], function(x) x$area_est_MEGB$mw) %>% do.call(cbind, .)
  area_est_EBP <- lapply(res_grouped[[setting_name]], function(x) x$area_est_ebp$mw) %>% do.call(cbind, .)
  area_est_MERF <- lapply(res_grouped[[setting_name]], function(x) x$area_est_MERF$mw) %>% do.call(cbind, .)
  true_values <- sapply(res_grouped[[setting_name]], function(x){
    
    x$sim_obj$pop_data |> group_by(idd) |> dplyr::summarise(true_mean = mean(Y)) %>% pull(true_mean)
    
  })
  
  QM_BoostMERT_L2   <- QualityMeasure(true_values, area_est_BoostMERT_L2,  MSETF = FALSE)
  QM_BoostREEM_L2   <- QualityMeasure(true_values, area_est_BoostREEM_L2,  MSETF = FALSE)
  QM_L2             <- QualityMeasure(true_values, area_est_L2,            MSETF = FALSE)
  QM_MBOOST_L2      <- QualityMeasure(true_values, area_est_MBOOST_L2,     MSETF = FALSE)
  QM_MBOOST_T_L2    <- QualityMeasure(true_values, area_est_MBOOST_T_L2,   MSETF = FALSE)
  QM_MEGB           <- QualityMeasure(true_values, area_est_MEGB,          MSETF = FALSE)
  QM_RF             <- QualityMeasure(true_values, area_est_RF,            MSETF = FALSE)
  QM_EBP            <- QualityMeasure(true_values, area_est_EBP,           MSETF = FALSE)
  QM_MERF           <- QualityMeasure(true_values, area_est_MERF,          MSETF = FALSE)

  QM_all <- list(
    EBP = QM_EBP,
    MERF = QM_MERF,
    MERTBoosting = QM_BoostMERT_L2,
    REEMBoosting = QM_BoostREEM_L2,
    GB = QM_L2,
    MBOOST = QM_MBOOST_L2,
    MBOOST_T = QM_MBOOST_T_L2,
    MEGB = QM_MEGB,
    RF = QM_RF
  )
  QM_all
  
})


names(qm_results) <- names(res_grouped)

# plots --------------------------------------------------------

library(data.table)
library(ggplot2)

plot_results <- function(result_list, quality_measure = "RRMSE") {
  

  dt <- rbindlist(
    lapply(names(result_list), function(nm) {
      if (!is.null(result_list[[nm]][[quality_measure]])) {
        data.table(
          estimator = nm,
          value = as.numeric(result_list[[nm]][[quality_measure]])
        )
      } else {
        NULL
      }
    }),
    use.names = TRUE, fill = TRUE
  )
  
  # Boxplot 
  ggplot(dt, aes(x = estimator, y = value, fill = estimator)) +
    geom_boxplot(alpha = 0.6) +
    theme_minimal(base_size = 14) +
    labs(
      x = "estimator",
      y = quality_measure
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


library(data.table)



out_dir <- "plots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


safe_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)


measures <- c("RB", "RRMSE")


if (is.null(names(qm_results)) || any(names(qm_results) == "")) {
  names(qm_results) <- paste0("result_", seq_along(qm_results))
}


for (nm in names(qm_results)) {
  result_list <- qm_results[[nm]]
  
  for (m in measures) {
   
    p <- plot_results(result_list, quality_measure = m)
    
    
    file_name <- sprintf("%s_%s.jpg", tolower(m), safe_name(nm))
    out_path <- file.path(out_dir, file_name)
    
    ggsave(filename = out_path, plot = p, width = 8, height = 5, dpi = 300)
  }
}






library(data.table)

extract_quality_df <- function(qm_results, quality_measure = "RRMSE") {

  all_data <- list()
  
  for (nm in names(qm_results)) {

    icc <- as.numeric(sub(".*ICC=([0-9.]+)_.*", "\\1", nm))
    D   <- as.numeric(sub(".*n_domains=([0-9]+)_.*", "\\1", nm))
    ni  <- as.numeric(sub(".*ni_domain=([0-9]+).*", "\\1", nm))
    

    for (est in names(qm_results[[nm]])) {
      vals <- qm_results[[nm]][[est]][[quality_measure]]
      
      if (!is.null(vals)) {
        all_data[[length(all_data) + 1]] <- data.table(
          ICC        = icc,
          n_domains  = D,
          ni_domain  = ni,
          estimator  = est,
          value      = as.numeric(vals)
        )
      }
    }
  }
  

  dt <- rbindlist(all_data)
  dt[]
}


df_rrmse <- extract_quality_df(qm_results, "RRMSE")
df_rb    <- extract_quality_df(qm_results, "RB")


library(ggplot2)
library(data.table)
library(viridis)
plot_quality_measure <- function(df, quality_measure = "RRMSE") {

  df <- as.data.table(df)
  

  out_dir <- "plots"
  dir.create(out_dir, showWarnings = FALSE)
  

  combos <- unique(df[, .(n_domains, ni_domain)])

  for (i in seq_len(nrow(combos))) {
    D_val  <- combos$n_domains[i]
    ni_val <- combos$ni_domain[i]
    
    d <- df[n_domains == D_val & ni_domain == ni_val]
    
    p <- ggplot(d, aes(x = factor(ICC), y = value, fill = estimator)) +
      geom_boxplot(alpha = 0.8, outlier.shape = 21, outlier.size = 0.8) +
      theme_bw(base_size = 13) +
      scale_fill_manual(values = color_mapping, drop = FALSE) +
      labs(
        x = "ICC",
        y = quality_measure,
        fill = "estimator"
      ) +
      theme(
        legend.position = "bottom",
        plot.title = element_text(face = "bold", hjust = 0.5),

        axis.text.x = element_text(
          vjust = 1.5,   
          margin = margin(t = 5) 
        )
      )
    

    ggsave(
      filename = sprintf("%s/%s_D%d_ni%d.jpg", out_dir, tolower(quality_measure), D_val, ni_val),
      plot = p, width = 10, height = 6, dpi = 300
    )
  }
}

df_rrmse <- extract_quality_df(qm_results, "RRMSE")
plot_quality_measure(df_rrmse, "RRMSE")

df_rb <- extract_quality_df(qm_results, "RB")
plot_quality_measure(df_rb, "RB")


# Tables ----------------------------------------------------------------





library(kableExtra)

for (nm in names(qm_results)) {
  df <- qm_results[[nm]]
  cap <- paste0("results for ", gsub("_", "\\\\_", nm, fixed = TRUE)) 
  
  print(
    kbl(
      df,
      format   = "latex",
      booktabs = TRUE,
      digits   = 3,
      caption  = cap,
      row.names = TRUE               
    ) |>
      kable_styling(full_width = FALSE, position = "center")
  )
}

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(kableExtra)


parse_key <- function(nm) {
  tibble(raw = nm) |>
    tidyr::extract(
      raw,
      into = c("ICC", "n_domains", "ni_domain"),
      regex = "^ICC=([^_]+)_n_domains=([^_]+)_ni_domain=([^_]+)$",
      remove = FALSE
    ) |>
    mutate(
      ICC_num = as.numeric(ICC),
      n_domains = as.integer(n_domains),
      ni_domain = as.integer(ni_domain)
    )
}


make_panel_table <- function(qm_results, n_dom, ni_dom, metrics = c("mw_rrmse","median_rrmse")) {
  

  meta <- map_chr(names(qm_results), identity) |> parse_key()
  keep_names <- meta %>%
    filter(n_domains == n_dom, ni_domain == ni_dom) %>%
    arrange(ICC_num) %>%
    pull(raw)
  
  stopifnot(length(keep_names) > 0)
  

  icc_vals <- keep_names %>% parse_key() %>% arrange(ICC_num) %>% pull(ICC)
  

  base_df <- qm_results[[keep_names[[1]]]] %>% as.data.frame %>%
    tibble::rownames_to_column("Model")
  

  wide <- base_df %>% select(Model)
  for (nm in keep_names) {
    df_icc <- qm_results[[nm]] %>% as.data.frame %>%
      tibble::rownames_to_column("Model") %>%
      select(Model, all_of(metrics))

    icc_lab <- parse_key(nm)$ICC[1]
    names(df_icc)[-1] <- paste0(icc_lab, "__", names(df_icc)[-1])
    wide <- wide %>%
      left_join(df_icc, by = "Model")
  }
  

  wide$Model <- kableExtra::escape_latex(wide$Model)
  

  col_order <- c("Model", unlist(lapply(icc_vals, function(iv) paste0(iv, "__", metrics))))
  wide <- wide[, col_order]
  

  num_cols <- setdiff(names(wide), "Model")
  for (cl in num_cols) {
    x <- wide[[cl]]
    mn <- min(x, na.rm = TRUE)
    as_chr <- sprintf("%.3f", x)
    idx <- which(x == mn)
    as_chr[idx] <- cell_spec(as_chr[idx], bold = TRUE)
    wide[[cl]] <- as_chr
  }
  

  header <- setNames(
    rep(length(metrics), length(icc_vals)),
    paste0("ICC = ", icc_vals)
  )
  

  cap <- paste0(
    "Results for n\\_domains=", n_dom,
    ", ni\\_domain=", ni_dom,
    " — ", paste(metrics, collapse = " & ")
  )
  
  # kable
  kbl(
    wide,
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    row.names = FALSE,
    caption = cap,
    align = c("l", rep("r", length(num_cols)))
  ) |>
    add_header_above(c(" " = 1, header), escape = FALSE) |>
    kable_styling(full_width = FALSE, position = "center")
}



all_meta <- map_chr(names(qm_results), identity) |> parse_key()
pairs <- all_meta %>%
  distinct(n_domains, ni_domain) %>%
  arrange(n_domains, ni_domain)


for (i in seq_len(nrow(pairs))) {
  nd <- pairs$n_domains[i]
  ni <- pairs$ni_domain[i]
  
  # RMSE-Panel
  print(make_panel_table(
    qm_results,
    n_dom = nd,
    ni_dom = ni,
    metrics = c("mw_rrmse", "median_rrmse")
  ))
  
  # RB-Panel
  print(make_panel_table(
    qm_results,
    n_dom = nd,
    ni_dom = ni,
    metrics = c("mw_rb", "median_rb")
  ))
}


runtime_summary <- lapply(names(res_grouped), function(setting_name) {

  time_elements <- grep("^time_", names(res_grouped[[setting_name]][[1]]), value = TRUE)
  

  sapply(time_elements, function(t_name) {
    t_values <- sapply(res_grouped[[setting_name]], function(x) abs(x[[t_name]]))
    stats <- quantile(t_values, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = F)
    c(
      min = stats[1],
      q25 = stats[2],
      median = stats[3],
      mean = mean(t_values, na.rm = F),
      q75 = stats[4],
      max = stats[5]
    ) 
  })
})

names(runtime_summary) <- names(res_grouped)
runtime_summary
