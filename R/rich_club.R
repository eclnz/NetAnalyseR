rich_club <- function(W, klevel = NULL) {
  validate_matrix(W)
  NofNodes <- ncol(W)
  NodeDegree <- colSums(W != 0)

  if (length(unique(NodeDegree))<=1){
    if (network_density(W)==1){
      stop("Network is fully connected, a rich club organisation cannot exist")
    } else{

    }
    stop("There are less than two unique node degrees, a rich club organisation cannot exist")
  }

  if (is.null(klevel)) {
    klevel <- max(NodeDegree)
  } else if (length(klevel) == 1) {
    klevel <- klevel
  } else {
    stop("Number of inputs incorrect. Should be [W], or [W, klevel]")
  }

  wrank <- sort(W, decreasing = TRUE)

  Rw <- numeric(klevel)

  for (kk in 1:klevel) {
    SmallNodes <- which(NodeDegree < kk)

    if (length(SmallNodes) == 0) {
      Rw[kk] <- NA
      next
    }

    CutoutW <- W[-SmallNodes, -SmallNodes]

    Wr <- sum(CutoutW)
    Er <- sum(CutoutW != 0)

    wrank_r <- wrank[1:Er]

    Rw[kk] <- Wr / sum(wrank_r)
  }
  while(length(Rw)<NofNodes){
    Rw <- c(Rw,NA)
  }
  return(Rw)
}

norm_rich_club <- function(W, n_rand=100, rand_array = NULL){
  # n_rand = 100
  # Generate randomized matrices if not provided
  if (is.null(rand_array)) {
    rand_array <- generateRewiredMatrices(W, n_rand) %>%
      abind::abind(along = 3)
  }
  rich <- rich_club(W)

  rand_rich <- apply(rand_array, FUN = rich_club, MARGIN = 3)
  # rand_rich[is.na(rand_rich)] <- Inf

  num_exceeding <- rowSums(as.matrix(rand_rich>rich)+0)
  rand_rich_sum <- apply(rand_rich, FUN = mean, MARGIN = 1)
  p_value = num_exceeding/n_rand
  norm_rich = rich/rand_rich_sum


  list(rich_club = rich/rand_rich_sum, p_value = p_value, node = 1:length(rich) )

}


grouping_list<- list(mTBI = "RUGLONG", rugby = c("Pre", "mid", "post"))
grouping_list<- list(Preseason = "Pre", Postseason = "post")

#' @title Rich Club Coefficient Comparison
#' @description Computes and compares rich club coefficients for multiple subjects grouped by given criteria.
#' @param mat_array An array of connectivity matrices where each slice corresponds to a different subject.
#' @param subjects A vector of subject identifiers.
#' @param grouping_list A list specifying how subjects should be grouped for analysis.
#' @param plot A boolean indicating whether to generate plots of the results.
#' @param permutation A boolean to choose whether permutation tests should be used instead of ANOVA for statistical comparison.
#' @param n_rand The number of random permutations to use if permutation is TRUE.
#' @param p_adj A boolean indicating whether significant differences between groups on the graph should be controlled for multiple comparisons.
#' @return A list containing statistical analysis results and optionally a plot, depending on the input parameters.
#' @examples
#' data_dir <- system.file("extdata", package = "NetAnalyseR")
#' subjects <- c("A", "B", "C", "D")
#' grouping <- list(AB = c("A" , "B"), CD = c("C", "D"))
#' file_convention <- ".csv"
#' output <- process_matrices(data_dir, subjects, file_convention)
#' compare_rich_club(mat_array = output$matrices,
#'                   subjects = output$subjects,
#'                   n_rand = 100,
#'                   grouping_list = grouping,
#'                   plot = T)
#'
#'
#'
#' @importFrom tidyr pivot_longer drop_na
#' @importFrom broom tidy
#' @importFrom stats aov
#' @importFrom dplyr summarise group_by n cur_data pull mutate n_distinct
#' @importFrom ggplot2 ggplot geom_line geom_errorbar geom_point theme_light ylab xlab
#' @importFrom purrr list_rbind
#' @export

compare_rich_club <- function(mat_array,
                              subjects,
                              grouping_list,
                              plot = FALSE,
                              permutation = FALSE,
                              n_rand = 100,
                              p_adj = TRUE) {
  subjects <- data.frame(subject = subjects)

  # Determining subjects to be analysed from grouping list
  grouped_subjects <- allocate_groups(subjects, grouping_list)
  include_index <- !is.na(grouped_subjects$group)

  if (!is.list(grouping_list)) {
    stop("Input 'grouping_list' must be a list.")
  }

  cat("Calculating Normalised Rich Club\n")

  norm_rich_club <- apply(mat_array[,,include_index], MARGIN = 3, FUN = function(mat) norm_rich_club(mat, n_rand))

  cat("Formulating Rich Club Data Statistics\n")
  long_norm_rich_club <- as.data.frame(norm_rich_club) %>%
    pivot_longer(
      cols = everything(),
      names_to = c("subject", ".value"),
      names_pattern = "(.*)\\.(.*)"
    ) %>%
    allocate_groups(grouping_list) %>%
    drop_na()

  permute_test <- function(data, n_perm = 100) {
    if (dplyr::n_distinct(data$group) <= 1) {
      stop("Not enough distinct groups for permutation testing.")
    }
    original_data <- data
    original_fit <- tidy(aov(rich_club ~ group, data = original_data))
    original_stat <- original_fit$statistic[1]  # assuming the first row is the group effect

    # Permutation logic
    perm_stats <- replicate(n_perm, {
      permuted_data <- data
      permuted_data$group <- sample(permuted_data$group)  # shuffle group labels
      permuted_fit <- tidy(aov(rich_club ~ group, data = permuted_data))
      permuted_stat <- permuted_fit$statistic[1]  # F-statistic for the group effect
    })

    p_value <- mean(perm_stats >= original_stat)
    return(data.frame(original_statistic = original_stat, p_value = p_value))
  }

  # Calculate results based on permutation or ANOVA
  results <- if (permutation) {
    long_norm_rich_club %>%
      group_by(node) %>%
      dplyr::filter(p_value<0.05) %>%
      dplyr::filter(n_distinct(group) > 1) %>%
      dplyr::filter(n()/n_distinct(group) > 1) %>% # Must be more than 1 observation in each group
      summarise(perm_results = list(permute_test(pick(rich_club))), .groups = 'drop') %>%
      ungroup() %>%
      mutate(p_value = sapply(perm_results, function(x) x$p_value),
             adjusted_p_value = p.adjust(p_value, method = "BH"))
  } else {
    long_norm_rich_club %>%
      group_by(node) %>%
      dplyr::filter(p_value<0.05) %>%
      dplyr::filter(n_distinct(group) > 1) %>% # Must be more than one group
      dplyr::filter(n()/n_distinct(group) > 1) %>% # Must be more than 1 observation in each group
      summarise(aov_results = list(tidy(aov(rich_club ~ group, data = pick(rich_club)))), .groups = 'drop') %>%
      ungroup() %>%
      mutate(p_value = sapply(aov_results, function(x) x$p.value[1]),  # Assuming first p-value is the group effect
             adjusted_p_value = p.adjust(p_value, method = "BH"))
  }
  # Prepare data for plotting and identify significant nodes
  data_for_plotting <- long_norm_rich_club %>%
    group_by(node, group) %>%
    summarise(
      mean_rich_club = mean(rich_club, na.rm = TRUE),
      sd = sd(rich_club, na.rm = TRUE),
      n = n(),
      sem = sd / sqrt(n),
      sd_lower = mean_rich_club - sd,
      sd_upper = mean_rich_club + sd,
      mean_p_value = mean(p_value)
    ) %>%
    dplyr::filter(dplyr::n()>1) %>% # Must be more than 1 observation in each group
    dplyr::filter(mean_p_value<0.05 & mean_rich_club > 1) %>%
    dplyr::filter(sd_lower>1) # Prevents the last few nodes from being significant sometimes.

  sig_nodes <- if (permutation) {
    # Extract node numbers for significant permutation results
    results %>%
      mutate(is_significant = if (p_adj) adjusted_p_value < 0.05 else p_value < 0.05) %>%
      dplyr::filter(is_significant) %>%
      pull(node)  # Extract node numbers
  } else {
    # Extract node numbers for significant ANOVA results
    results %>%
      mutate(is_significant = if (p_adj) adjusted_p_value < 0.05 else p_value < 0.05) %>%
      dplyr::filter(is_significant) %>%
      pull(node)  # Extract node numbers
  }

  if (plot) {
    rich_plot <- long_norm_rich_club %>%
      group_by(node, group) %>%
      summarise(mean_rich_club = mean(rich_club, na.rm = TRUE)) %>%
      ggplot(aes(x = as.numeric(node), y = mean_rich_club, color = group)) +
      geom_line(linewidth = 1) +
      geom_errorbar(data = data_for_plotting, aes(ymin = sd_lower, ymax = sd_upper), width = 0.8) +
      geom_point(data = long_norm_rich_club[long_norm_rich_club$node %in% sig_nodes,] %>% group_by(node) %>% summarize(mean_rich_club = mean(rich_club)), aes(x = as.numeric(node), y = mean_rich_club+0.04*mean_rich_club), color = "red") +
      theme_light() +
      ylab("Normalised Rich Club Coefficient") +
      xlab("Node Degree")
  }
  return(
    if(permutation==FALSE){
      list(results$aov_results %>%
        list_rbind() %>%
          cbind(node = results$node),
      rich_plot)
    } else{
      list(results$perm_results %>%
        list_rbind() %>%
          cbind(node = results$node),
        rich_plot)
    }
  )
}
