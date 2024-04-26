W <- matrix(c(0, 2, 1, 0, 2, 0, 3, 5, 1, 3, 0, 6, 0, 5, 6, 0), nrow = 4, byrow = TRUE)

# Define test networks
single_node <- matrix(0, nrow = 1, ncol = 1)
empty_network <- matrix(0, nrow = 3, ncol = 3)
disconnected_network <- matrix(c(0, 2, 0, 0, 2, 0, 0, 5, 0, 0, 0, 0, 0, 5, 0, 0), nrow = 4, byrow = TRUE)
fully_connected <- matrix(1, nrow = 4, ncol = 4) - diag(4)
network_with_loops <- matrix(c(1,0,0, 0,1,0, 0,0,1), nrow = 3, byrow = TRUE)

test_network_metric <- function(metric_func, expected, networks, description) {
  test_that(paste("Test", description), {
    for (i in seq_along(networks)) {
      result <- metric_func(networks[[i]])
      expect_equal(result, expected[[i]])
    }
  })
}

expected_betweenness <- list(c(0), c(0, 0, 0), c(0, 0, 0, 0, 0), c(2, 2, 2, 2), c(1, 1, 1), c(0, 0.5, 0))

test_network_metric(betweenness_wei, expected_betweenness,
                    list(single_node, empty_network, disconnected_network, fully_connected, network_with_loops, weighted_network),
                    "betweenness centrality")


test_that("global_efficiency_wei", {
  expect_equal(
    global_efficiency_wei(W), 3.10476190)
  })

test_that("global_clustering_coefficient_wei", {
  expect_equal(
    global_clustering_coefficient_wei(W), 2.62438556)
  })
