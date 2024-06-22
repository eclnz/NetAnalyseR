#' @title Group Comparisons for Metrics with Visualization
#' @description Performs ANOVA and pairwise t-tests for specified metrics across groups, adjusting p-values and generating combined plots.
#' @param group_df Data frame containing 'subject', 'group', and metric columns.
#' @param metrics Character vector of metric column names to analyze.
#' @param comparisons List of pairwise group comparisons for t-tests.
#' @param p_adjust_method Character string specifying the p-value adjustment method (default is "BH").
#' @param group_order Optional character vector specifying the order of groups.
#' @param test_type Optional character vector specifying whether the tests should be parametric or nonparametric
#' @return A list with ANOVA results, t-test results, and a combined plot.
#' @importFrom stats as.formula kruskal.test t.test wilcox.test
#' @importFrom ggplot2 ggplot
#' @importFrom ggsignif geom_signif
#' @importFrom cowplot ggdraw draw_plot draw_label
#' @export
group_comparisons <- function(group_df, metrics, comparisons, p_adjust_method = "BH", group_order = NULL, test_type = "parametric") {

  # Internal function to transform metric names
  transform_metric_name <- function(name) {
    name <- gsub("_wei$", "", name)  # Remove '_wei' at the end
    name <- gsub("_", " ", name)     # Replace underscores with spaces
    name <- str_to_title(name)       # Capitalize
    return(name)
  }

  # Check for 'subject' and 'group' columns
  required_columns <- c("subject", "group")
  if (!all(required_columns %in% colnames(group_df))) {
    stop("Group dataframe must have 'subject' and 'group' columns.")
  }

  group_df$group <- as.factor(group_df$group)
  if (!is.null(group_order)) {
    if (!all(group_order %in% levels(group_df$group))) {
      stop("All specified group_order values must be present in the 'group' column.")
    }
    group_df$group <- factor(group_df$group, levels = group_order)
  } else {
    group_df$group <- factor(group_df$group)
  }

  # Check if all metrics exist in the data frame
  missing_metrics <- setdiff(metrics, colnames(group_df))
  if (length(missing_metrics) > 0) {
    stop("The following metrics are missing in the data frame: ", paste(missing_metrics, collapse = ", "))
  }

  # Check if all groups in comparisons exist in the data frame
  all_groups <- unique(unlist(comparisons))
  missing_groups <- setdiff(all_groups, levels(group_df$group))
  if (length(missing_groups) > 0) {
    stop("The following groups are missing in the 'group' column: ", paste(missing_groups, collapse = ", "))
  }

  # Initialize lists to store results
  test_results <- list()
  t_test_results_list <- list()
  plot_list <- list()  # List to store individual plots

  # Iterate over each metric
  for (index in seq_along(metrics)) {
    metric <- metrics[index]
    # Define the formula dynamically
    formula <- as.formula(paste(metric, "~ group"))

    # Filter the data
    filtered_df <- group_df %>% filter(group %in% all_groups)

    # Perform statistical tests based on test_type
    if (test_type == "parametric") {
      # Perform ANOVA for the current metric
      test_result <- aov(formula, data = filtered_df) %>% tidy()
    } else if (test_type == "nonparametric") {
      # Perform Kruskal-Wallis test for the current metric
      test_result <- kruskal.test(formula, data = filtered_df) %>% tidy()
    } else {
      stop("Invalid test_type. Use 'parametric' or 'nonparametric'.")
    }
    test_results[[metric]] <- test_result

    # Initialize list to store t-test or Wilcoxon test results for this metric
    comparison_results <- list()

    # Gather p-values for adjustment
    p_values <- numeric()

    # Perform pairwise tests for each comparison
    for (i in seq_along(comparisons)) {
      comparison <- comparisons[[i]]
      comparison_label <- paste(comparison, collapse = "_vs_")

      # Filter the data for the comparison
      comparison_df <- filtered_df %>% filter(group %in% comparison)

      # Perform pairwise test based on test_type
      if (test_type == "parametric") {
        pairwise_test_result <- t.test(as.formula(paste(metric, "~ group")), data = comparison_df) %>% tidy()
      } else if (test_type == "nonparametric") {
        pairwise_test_result <- wilcox.test(as.formula(paste(metric, "~ group")), data = comparison_df) %>% tidy()
      }

      # Add the p-value to results
      p_values <- c(p_values, pairwise_test_result$p.value)
      comparison_results[[comparison_label]] <- pairwise_test_result
    }

    # Adjust p-values
    adjusted_p_values <- p.adjust(p_values, method = p_adjust_method)

    # Incorporate adjusted p-values into results
    for (i in seq_along(comparison_results)) {
      comparison_results[[i]]$p.adjusted <- adjusted_p_values[i]
    }

    t_test_results_list[[metric]] <- comparison_results

    # Calculate y position for significance annotations
    y_max <- max(filtered_df[[metric]], na.rm = TRUE)
    signif_y_position <- y_max + -0.01 * y_max  # Closer to the data points

    # Create a plot for each metric
    y_label <- transform_metric_name(metric)
    plot <- ggplot(filtered_df, aes(x = group, y = .data[[metric]], color = group)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.6) +
      labs(title = NULL, y = y_label, x = NULL) +
      scale_x_discrete(labels = function(x) str_to_title(x)) +
      theme_minimal() +
      theme(legend.position = "none",
            plot.margin = unit(c(2, 0.5, 0.5, 0.5), "lines"))  # Add margin to the top

    # Add significance annotations
    for (i in seq_along(comparisons)) {
      comparison <- comparisons[[i]]
      comparison_label <- paste(comparison, collapse = "_vs_")
      if (comparison_results[[comparison_label]]$p.adjusted < 0.05) {
        plot <- plot +
          ggsignif::geom_signif(
            comparisons = list(comparison),
            annotations = sprintf("p = %.3f", comparison_results[[comparison_label]]$p.adjusted),
            map_signif_level = FALSE,
            textsize = 3,
            color = "slategray",
            y_position = signif_y_position + (signif_y_position* ((i-1)/75))
          )
      }
    }

    # Add the title annotation above the plot in the top left corner
    plot_with_title <- cowplot::ggdraw() +
      cowplot::draw_plot(plot) +
      cowplot::draw_label(LETTERS[index], fontface = 'bold', size = 18, x = 0, y = 1, hjust = -0.1, vjust = 1.1)

    plot_list[[metric]] <- plot_with_title
  }

  # Combine plots using grid.arrange
  combined_plot <- do.call(grid.arrange, c(plot_list, ncol = length(metrics)))

  # Return results as a list
  return(list(test_results = test_results, t_tests = t_test_results_list, combined_plot = combined_plot))
}
