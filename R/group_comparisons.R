#' @title Group Statistics Function
#' @description Performs statistical analysis on grouped data, including normality tests, homogeneity tests, and post-hoc comparisons.
#' @param group_df Data frame containing the data to be analyzed, including grouping variables.
#' @param metrics Character vector of metric names to be analyzed.
#' @param comparisons List of group comparisons to be performed.
#' @param p_adjust_method Method for p-value adjustment (default is "BH").
#' @param test_type Type of test to be used (parametric, nonparametric, or permutation). If NULL, it will be determined based on data.
#' @param num_permutations Number of permutations for permutation test (default is 1000).
#' @param analysis_type Type of analysis: "between" for between-group comparisons, "within" for within-subject comparisons.
#' @param id_var Identifier variable for within-subject comparisons.
#' @param time_var Time variable for within-subject comparisons.
#' @importFrom stats TukeyHSD as.formula friedman.test kruskal.test t.test wilcox.test
#' @importFrom dunn.test dunn.test
#' @importFrom utils combn
#' @return List containing test results, post-hoc comparisons, and normality/homogeneity test results.

group_statistics <- function(group_df, metrics, comparisons, p_adjust_method = "BH", test_type = NULL, num_permutations = 1000, analysis_type = "between", id_var = NULL, time_var = NULL) {

  determine_test_type <- function(metric_data, group_variable) {
    group_variable <- as.factor(group_variable)
    shapiro_results <- by(metric_data, group_variable, function(x) stats::shapiro.test(x))
    shapiro_p <- sapply(shapiro_results, function(x) x$p.value)
    shapiro_overall_result <- all(shapiro_p >= 0.05)
    levene_test <- car::leveneTest(metric_data ~ group_variable)
    levene_overall_result <- levene_test$`Pr(>F)`[1] >= 0.05
    test_type <- if (shapiro_overall_result && levene_overall_result) "parametric" else "nonparametric"
    return(list(test_type = test_type, shapiro_p_values = shapiro_p, levene_test = levene_test))
  }

  permutation_test <- function(formula, data, test_function, num_permutations = 1000) {
    observed_result <- test_function(formula, data = data)
    observed_stat <- as.numeric(tidy(observed_result)$statistic)[1]

    perm_stats <- numeric(num_permutations)
    for (i in 1:num_permutations) {
      permuted_data <- data
      permuted_data[[all.vars(formula)[2]]] <- sample(permuted_data[[all.vars(formula)[2]]])
      perm_result <- test_function(formula, permuted_data)
      perm_stats[i] <- as.numeric(tidy(perm_result)$statistic)[1]
    }
    p_value <- mean(perm_stats >= observed_stat)
    return(p_value)
  }

  filter_complete_cases <- function(df, id_var, time_var) {
    complete_subjects <- df %>%
      group_by(across(all_of(id_var))) %>%
      summarize(n_timepoints = n_distinct(!!sym(time_var))) %>%
      filter(n_timepoints == length(unique(df[[time_var]]))) %>%
      pull(!!sym(id_var))
    return(df %>% filter(!!sym(id_var) %in% complete_subjects))
  }

  within_subject_analysis <- function(metric, df, id_var, time_var) {
    df <- filter_complete_cases(df, id_var, time_var)

    if (nrow(df) == 0) {
      stop("No complete cases found. Check your data or filtering criteria.")
    }

    df[[id_var]] <- as.factor(df[[id_var]])
    df[[time_var]] <- as.factor(df[[time_var]])
    df[[id_var]] <- droplevels(df[[id_var]])
    df[[time_var]] <- droplevels(df[[time_var]])

    # Repeated Measures ANOVA
    formula_str <- paste(metric, "~", time_var, "+ Error(", id_var, "/", time_var, ")")
    aov_result <- aov(as.formula(formula_str), data = df)
    tidy_aov <- tidy(aov_result)

    if (is.null(tidy_aov$p.value) || all(is.na(tidy_aov$p.value))) {
      # Friedman Test if ANOVA is not applicable
      formula_friedman <- as.formula(paste(metric, "~", time_var, "|", id_var))
      friedman_result <- friedman.test(formula_friedman, data = df)
      tidy_friedman <- tidy(friedman_result)
      return(tidy_friedman)
    } else {
      return(tidy_aov)
    }
  }

  post_hoc_comparisons <- function(data, metric, id_var = NULL, time_var = NULL, test_type, analysis_type, p_adjust_method = "BH") {
    results <- list()

    if (analysis_type == "within") {
      levels_time <- factor(levels(data[[time_var]]))
      comb <- combn(levels_time, 2)

      for (i in 1:ncol(comb)) {
        comparison_label <- paste(comb[, i], collapse = "_vs_")
        comparison_data <- data %>% filter(!!sym(time_var)%in% comb[, i])

        if (test_type == "parametric") {
          # Paired t-test for parametric data
          pairwise_test_result <- t.test(as.formula(paste(metric, "~", time_var)), paired = TRUE, data = comparison_data) %>% tidy()
        } else if (test_type == "nonparametric") {
          # Wilcoxon signed-rank test for non-parametric data
          pairwise_test_result <- wilcox.test(as.formula(paste(metric, "~", time_var)), paired = TRUE, data = comparison_data) %>% tidy()
        } else {
          stop("Invalid test_type. Use 'parametric' or 'nonparametric'.")
        }

        results[[comparison_label]] <- pairwise_test_result
      }

    } else if (analysis_type == "between") {
      for (comparison in comparisons) {
        comparison_label <- paste(comparison, collapse = "_vs_")
        comparison_data <- data %>% filter(group %in% comparison)

        if (length(unique(comparison_data$group)) < 2) {
          warning(paste("Comparison:", comparison_label, "- requires at least two group levels. Skipping..."))
          next
        }

        if (test_type == "parametric") {
          # Tukey HSD for parametric data
          aov_model <- aov(as.formula(paste(metric, "~ group")), data = comparison_data)
          tukey_result <- TukeyHSD(aov_model)
          tidy_tukey <- broom::tidy(tukey_result)
          results[[comparison_label]] <- tidy_tukey
        } else if (test_type == "nonparametric") {
          # Dunn's test for nonparametric data
          dunn_test_result <- dunn.test::dunn.test(
            x = comparison_data[[metric]],
            g = comparison_data$group,
            altp = TRUE
          )

          tidy_dunn <- data.frame(
            comparison = dunn_test_result$comparisons,
            p.value = dunn_test_result$altP
          )
          results[[comparison_label]] <- tidy_dunn
        } else {
          stop("Invalid test_type. Use 'parametric' or 'nonparametric'.")
        }
      }
    }

    # # Adjust p-values
    # p_values <- sapply(results, function(x) x$p.value)
    # adjusted_p_values <- p.adjust(p_values, method = p_adjust_method)
    #
    # for (i in seq_along(results)) {
    #   results[[i]]$p.adjusted <- adjusted_p_values[i]
    # }

    return(results)
  }

  if(is.null(group_df$group)){
    stop('dataframe must have variable: group')
  }

  missing_metrics <- metrics[!metrics %in% colnames(group_df)]

  if (length(missing_metrics) > 0) {
    stop(paste('The dataframe is missing the following metrics:', paste(missing_metrics, collapse = ", ")))
  }

  if (analysis_type == "within") {
    if (is.null(id_var) || is.null(time_var)) {
      stop("Both id_var and time_var must be supplied when performing within-subject comparisons.")
    }
    # If not already done, remove incomplete observations  without all time points and those not in comparisons
    all_comparison_levels <- unique(unlist(comparisons))
    group_df <- group_df %>%
      filter(group %in% all_comparison_levels) %>%
      filter_complete_cases(id_var, time_var)
    if (length(unique(group_df[[time_var]]))<2){
      stop("Less than 2 unique time_var groups in data")
    }
  }
  if (analysis_type == "between") {
    if(!is.null(id_var) || !is.null(time_var)){
      stop("When performing between-subject comparisons, id_var and time_var cannot be supplied.")
    }
    # If not already done, remove observations not within comparisons
    group_df <- group_df %>% filter(group %in% unique(unlist(comparisons)))
    if (length(unique(group_df$group))<2){
      stop("Less than 2 unique groups in data")
    }
  }

  test_results <- list()
  post_hoc_results <- list()
  normality_homogeneity_results <- list()

  for (metric in metrics) {
    formula <- as.formula(paste(metric, "~ group"))

    metric_group_df <- group_df %>% mutate_if(is.factor, droplevels)


    if (analysis_type == "within") {
      test_info <- determine_test_type(metric_group_df[[metric]], metric_group_df[[time_var]])
    }
    else if (analysis_type == "between") {
      test_info <- determine_test_type(metric_group_df[[metric]], metric_group_df$group)
    }
    normality_homogeneity_results[[metric]] <- test_info
    if (is.null(test_type)) {
      test_type <- test_info$test_type
    }


    if (analysis_type == "between") {
      if (test_type == "parametric") {
        overall_test_result <- aov(formula, data = metric_group_df) %>% tidy()
        post_hoc <- post_hoc_comparisons(metric_group_df, metric, NULL, NULL, test_type, "between", p_adjust_method)
      } else if (test_type == "nonparametric") {
        overall_test_result <- kruskal.test(formula, data = metric_group_df) %>% tidy()
        post_hoc <- post_hoc_comparisons(metric_group_df, metric, NULL, NULL, test_type, "between", p_adjust_method)
      } else if (test_type == "permutation") {
        p_value <- permutation_test(formula, metric_group_df, kruskal.test, num_permutations)
        overall_test_result <- data.frame(term = metric, statistic = NA, p.value = p_value, method = "Permutation Test")
        post_hoc <- NULL
      } else {
        stop("Invalid test_type. Use 'parametric', 'nonparametric', or 'permutation'.")
      }
    } else if (analysis_type == "within") {
      overall_test_result <- within_subject_analysis(metric, metric_group_df, id_var, time_var)
      post_hoc <- NULL
      if (any(overall_test_result$p.value < 0.05, na.rm = TRUE)) {
        post_hoc <- post_hoc_comparisons(metric_group_df, metric, id_var, time_var, test_type, "within", p_adjust_method)
      }
    } else {
      stop("Invalid analysis_type. Use 'between' or 'within'.")
    }

    test_results[[metric]] <- overall_test_result
    if(!is.null(post_hoc)){
      post_hoc_results[[metric]] <- list_rbind(post_hoc)
    }
  }

  return(list(test_results = test_results, post_hoc = post_hoc_results, normality_homogeneity = normality_homogeneity_results))
}

#' @title Plot Group Statistics with Significance Annotations
#' @description Generates box plots with jitter for specified metrics across groups, with significance annotations based on adjusted p-values from pairwise comparisons.
#' @param group_df Data frame containing 'subject', 'group', and metric columns.
#' @param metrics Character vector of metric column names to plot.
#' @param comparisons List of pairwise group comparisons for significance annotations.
#' @param t_test_results_list List of t-test or Wilcoxon test results with adjusted p-values.
#' @param group_order Optional character vector specifying the order of groups.
#' @return A combined plot with significance annotations for each metric.
#' @importFrom ggplot2 ggplot
#' @importFrom ggsignif geom_signif
#' @importFrom cowplot ggdraw draw_plot draw_label
#' @export
group_statistics_plots <- function(group_df, metrics, comparisons, stats_results, group_order = NULL) {
  plot_list <- list()

  transform_metric_name <- function(name) {
    name <- gsub("_wei$", "", name)  # Remove '_wei' at the end
    name <- gsub("_", " ", name)     # Replace underscores with spaces
    name <- str_to_title(name)       # Capitalize
    return(name)
  }

  for (index in seq_along(metrics)) {
    metric <- metrics[index]
    filtered_df <- group_df
    y_label <- transform_metric_name(metric)

    y_max <- max(filtered_df[[metric]], na.rm = TRUE)
    y_min <- min(filtered_df[[metric]], na.rm = TRUE)

    signif_y_position <- y_max + ((y_max - y_min)* 0.02)

    plot <- ggplot(filtered_df, aes(x = group, y = .data[[metric]], color = group)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.6) +
      labs(title = NULL, y = y_label, x = NULL) +
      scale_x_discrete(labels = function(x) str_to_title(x)) +
      theme_minimal() +
      theme(legend.position = "none",
            plot.margin = unit(c(2, 0.5, 0.5, 0.5), "lines"))

    for (i in seq_along(comparisons)) {
      comparison <- comparisons[[i]]
      comparison_label <- paste(comparison, collapse = "_vs_")
      if (metric %in% names(stats_results$post_hoc)) {
        post_hoc_results <- stats_results$post_hoc[[metric]]
        p_value <- post_hoc_results$adj.p.value[[i]]
        if (!is.null(p_value) && p_value < 0.05) {
          plot <- plot + ggsignif::geom_signif(
              comparisons = list(comparison),
              annotations = sprintf("p = %.3f", p_value),
              map_signif_level = FALSE,
              textsize = 3,
              color = "slategray",
              y_position = signif_y_position
            )
        }
      }
    }

    plot_with_title <- cowplot::ggdraw() +
      cowplot::draw_plot(plot) +
      cowplot::draw_label(LETTERS[index], fontface = 'bold', size = 18, x = 0, y = 1, hjust = -0.1, vjust = 1.1)

    plot_list[[metric]] <- plot_with_title
  }

  combined_plot <- do.call(grid.arrange, c(plot_list, ncol = length(metrics)))

  return(combined_plot)
}

#' @title Group Comparisons for Metrics with Dynamic Testing and Visualization
#' @description Performs statistical tests on specified metrics across groups, adjusting p-values and generating combined plots. Supports dynamic determination of test type.
#' @param group_df Data frame containing 'subject', 'group', and metric columns.
#' @param metrics Character vector of metric column names to analyze.
#' @param comparisons List of pairwise group comparisons for tests or 'all' to include all possible pairwise comparisons.
#' @param p_adjust_method Character string specifying the p-value adjustment method (default is "BH").
#' @param group_order Optional character vector specifying the order of groups.
#' @param test_type Character string specifying "parametric" or "nonparametric" (optional). If NULL, test type is determined by normality and homogeneity tests.
#' @return A list with overall test results, pairwise test results, normality and homogeneity results, and a combined plot.
#' @export
group_comparisons <- function(group_df, metrics, comparisons, p_adjust_method = "BH", group_order = NULL, test_type = NULL, num_permutations = 1000, analysis_type = "between", id_var = NULL, time_var = NULL) {
  # All unique levels of comparison
  comparison_groups <- unique(unlist(comparisons))

  # If 'all' is specified, all possible combinations are used.
  if(is.character(comparisons)){
    if(comparisons == 'all'){
      comparisons <- utils::combn(comparison_groups, 2, simplify = FALSE)
    } else {
      stop("Comparisons must be list of comparisons or 'all'")
    }
  } else if(!is.list(comparisons)){
    if (length(comparison_groups)<1){
      stop("At least 1 comparison needs to be specified")
    }
  }

  # Filter data to only include unique levels of comparisons
  group_df <- group_df %>% filter(group %in% comparison_groups)

  # Format group as factor
  group_df$group <- as.factor(group_df$group)

  required_columns <- c("subject", "group")
  if (!all(required_columns %in% colnames(group_df))) {
    stop("Group dataframe must have 'subject' and 'group' columns.")
  }

  missing_metrics <- setdiff(metrics, colnames(group_df))
  if (length(missing_metrics) > 0) {
    stop("The following metrics are missing in the data frame: ", paste(missing_metrics, collapse = ", "))
  }

  comparison_groups <- unique(unlist(comparisons))
  missing_groups <- setdiff(comparison_groups, levels(group_df$group))
  if (length(missing_groups) > 0) {
    stop("The following groups are missing in the 'group' column: ", paste(missing_groups, collapse = ", "))
  }

  if (!is.null(group_order)) {
    if (!all(group_order %in% comparison_groups)) {
      stop("All specified group_order values must be present in the comparisons")
    }
    group_df$group <- factor(group_df$group, levels = group_order)
  }

  stats_results <- group_statistics(group_df, metrics, comparisons, p_adjust_method, test_type, num_permutations, analysis_type, id_var, time_var)
  combined_plot <- group_statistics_plots(group_df, metrics, comparisons, stats_results, group_order)

  return(list(
    test_results = stats_results$test_results,
    t_tests = stats_results$t_tests,
    normality_homogeneity = stats_results$normality_homogeneity,
    combined_plot = combined_plot
  ))
}
