#' @title Group Statistics for Metrics with Parametric and Nonparametric Tests
#' @description Performs ANOVA or Kruskal-Wallis tests on specified metrics across groups, adjusting p-values and considering normality and homogeneity.
#' @param group_df Data frame containing 'subject', 'group', and metric columns.
#' @param metrics Character vector of metric column names to analyze.
#' @param comparisons List of pairwise group comparisons for t-tests or Wilcoxon tests.
#' @param p_adjust_method Character string specifying the p-value adjustment method (default is "BH").
#' @param test_type Character string specifying "parametric" or "nonparametric" (optional). If NULL, test type is determined by Shapiro-Wilk and Levene's tests.
#' @importFrom stats as.formula kruskal.test t.test wilcox.test shapiro.test
#' @importFrom car leveneTest
#' @return A list with test results, pairwise test results, and normality and homogeneity test results.
group_statistics <- function(group_df, metrics, comparisons, p_adjust_method = "BH", test_type = NULL) {
  determine_test_type <- function(metric_data, group_variable) {
    # Ensure group is factor
    group_variable <- as.factor(group_variable)

    # Shapiro-Wilk test for normality
    shapiro_results <- by(metric_data, group_variable, function(x) stats::shapiro.test(x))
    shapiro_p <- sapply(shapiro_results, function(x) x$p.value)
    shapiro_overall_result <- all(shapiro_p >= 0.05)

    # Levene's test for homogeneity of variance
    levene_test <- car::leveneTest(metric_data ~ group_variable)
    levene_overall_result <- levene_test$`Pr(>F)`[1] >= 0.05

    # Determine test type
    test_type <- if (shapiro_overall_result && levene_overall_result) "parametric" else "nonparametric"

    return(list(test_type = test_type, shapiro_p_values = shapiro_p, levene_test = levene_test))
  }

  test_results <- list()
  t_test_results_list <- list()
  normality_homogeneity_results <- list()  # Store results of normality and homogeneity tests

  for (metric in metrics) {
    formula <- as.formula(paste(metric, "~ group"))
    filtered_df <- group_df

    if (is.null(test_type)) {
      test_info <- determine_test_type(filtered_df[[metric]], filtered_df$group)
      test_type <- test_info$test_type
      normality_homogeneity_results[[metric]] <- test_info
    }

    if (test_type == "parametric") {
      test_result <- stats::aov(formula, data = filtered_df) %>% tidy()
    } else if (test_type == "nonparametric") {
      test_result <- stats::kruskal.test(formula, data = filtered_df) %>% tidy()
    } else {
      stop("Invalid test_type. Use 'parametric' or 'nonparametric'.")
    }
    test_results[[metric]] <- test_result

    comparison_results <- list()
    p_values <- numeric()

    for (comparison in comparisons) {
      comparison_label <- paste(comparison, collapse = "_vs_")
      comparison_df <- filtered_df %>% filter(group %in% comparison)

      if (length(unique(comparison_df$group)) < 2) {
        warning(paste("Comparison:", comparison_label, "- requires at least two group levels. Skipping..."))
        next
      }

      if (test_type == "parametric") {
        pairwise_test_result <- t.test(as.formula(paste(metric, "~ group")), data = comparison_df) %>% tidy()
      } else if (test_type == "nonparametric") {
        pairwise_test_result <- wilcox.test(as.formula(paste(metric, "~ group")), data = comparison_df) %>% tidy()
      }

      p_values <- c(p_values, pairwise_test_result$p.value)
      comparison_results[[comparison_label]] <- pairwise_test_result
    }

    adjusted_p_values <- p.adjust(p_values, method = p_adjust_method)

    for (i in seq_along(comparison_results)) {
      comparison_results[[i]]$p.adjusted <- adjusted_p_values[i]
    }

    t_test_results_list[[metric]] <- comparison_results
  }

  return(list(test_results = test_results, t_tests = t_test_results_list, normality_homogeneity = normality_homogeneity_results))
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
group_statistics_plots <- function(group_df, metrics, comparisons, t_test_results_list, group_order = NULL) {
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
    signif_y_position <- y_max + 0.05 * y_max

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
      if (t_test_results_list[[metric]][[comparison_label]]$p.adjusted < 0.05) {
        plot <- plot +
          ggsignif::geom_signif(
            comparisons = list(comparison),
            annotations = sprintf("p = %.3f", t_test_results_list[[metric]][[comparison_label]]$p.adjusted),
            map_signif_level = FALSE,
            textsize = 3,
            color = "slategray",
            y_position = signif_y_position + (signif_y_position * ((i - 1) / 75))
          )
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
group_comparisons <- function(group_df, metrics, comparisons, p_adjust_method = "BH", group_order = NULL, test_type = NULL) {
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

  stats_results <- group_statistics(group_df, metrics, comparisons, p_adjust_method, test_type)
  combined_plot <- group_statistics_plots(group_df, metrics, comparisons, stats_results$t_tests, group_order)

  return(list(
    test_results = stats_results$test_results,
    t_tests = stats_results$t_tests,
    normality_homogeneity = stats_results$normality_homogeneity,
    combined_plot = combined_plot
  ))
}
