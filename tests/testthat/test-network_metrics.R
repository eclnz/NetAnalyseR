# ---------------------------------------------------------------------------
# Test data shared across test blocks
# ---------------------------------------------------------------------------

W <- matrix(c(0, 2, 1, 0,
              2, 0, 3, 5,
              1, 3, 0, 6,
              0, 5, 6, 0), nrow = 4, byrow = TRUE)

# 3-node path graph: 1-2-3, all weights = 1
# distances: d(1,2)=1, d(2,3)=1, d(1,3)=2 => CPL = (1+2+1)/3 = 4/3
path3 <- matrix(c(0, 1, 0,
                  1, 0, 1,
                  0, 1, 0), nrow = 3, byrow = TRUE)

# Fully connected 4-node graph, unit weights (no self-loops)
full4 <- matrix(1, 4, 4)
diag(full4) <- 0

# Networks used in legacy betweenness tests
single_node        <- matrix(0, nrow = 1, ncol = 1)
empty_network      <- matrix(0, nrow = 3, ncol = 3)
disconnected_network <- matrix(c(0, 2, 0, 0,
                                  2, 0, 0, 5,
                                  0, 0, 0, 0,
                                  0, 5, 0, 0), nrow = 4, byrow = TRUE)
fully_connected    <- matrix(1, nrow = 4, ncol = 4) - diag(4)
network_with_loops <- matrix(c(1, 0, 0,
                                0, 1, 0,
                                0, 0, 1), nrow = 3, byrow = TRUE)
# 3-node weighted network used in legacy betweenness test
weighted_network <- matrix(c(0, 1, 0,
                              1, 0, 1,
                              0, 1, 0), nrow = 3, byrow = TRUE)

# ---------------------------------------------------------------------------
# 1. validate_matrix
# ---------------------------------------------------------------------------

test_that("validate_matrix errors on non-matrix input", {
  expect_error(validate_matrix(1:4), "must be a matrix")
})

test_that("validate_matrix errors on non-square matrix", {
  expect_error(validate_matrix(matrix(1:6, nrow = 2)), "square")
})

test_that("validate_matrix errors on matrix with NA values", {
  m <- matrix(c(0, 1, 1, NA), nrow = 2)
  expect_error(validate_matrix(m), "NA")
})

test_that("validate_matrix warns and symmetrizes an asymmetric matrix", {
  m <- matrix(c(0, 2, 1, 0,
                3, 0, 4, 0,
                1, 4, 0, 5,
                0, 0, 5, 0), nrow = 4, byrow = TRUE)
  expect_warning(out <- validate_matrix(m), "symmetric")
  expect_true(isSymmetric(out))
})

test_that("validate_matrix passes a valid symmetric matrix unchanged", {
  out <- validate_matrix(W)
  expect_equal(out, W)
})

# ---------------------------------------------------------------------------
# 2. network_density
# ---------------------------------------------------------------------------

test_that("network_density of fully connected 4-node graph is 1", {
  expect_equal(network_density(full4), 1)
})

test_that("network_density of path3 is 2/3", {
  # Lower triangle has 3 pairs; 2 are connected
  expect_equal(network_density(path3), 2 / 3)
})

test_that("network_density of W is between 0 and 1", {
  d <- network_density(W)
  expect_true(d >= 0 && d <= 1)
})

# ---------------------------------------------------------------------------
# 3. node_strength
# ---------------------------------------------------------------------------

test_that("node_strength equals column sums with diagonal zeroed", {
  expected <- colSums(W) - diag(W)   # W has zero diagonal so same as colSums
  expect_equal(node_strength(W), expected)
})

test_that("node_strength of fully connected 4-node unit graph is 3 for all nodes", {
  expect_equal(node_strength(full4), rep(3, 4))
})

# ---------------------------------------------------------------------------
# 4. characteristic_path_length — 3-node path graph
# ---------------------------------------------------------------------------

test_that("characteristic_path_length of path3 is 4/3", {
  # d(1,2)=1, d(2,3)=1, d(1,3)=2 => mean of lower triangle = (1+2+1)/3
  expect_equal(characteristic_path_length(path3), 4 / 3)
})

# ---------------------------------------------------------------------------
# 5. global_efficiency_wei (legacy test preserved)
# ---------------------------------------------------------------------------

test_that("global_efficiency_wei", {
  expect_equal(
    global_efficiency_wei(W), 3.10476190)
})

# ---------------------------------------------------------------------------
# 6. global_clustering_coefficient_wei (legacy test preserved)
# ---------------------------------------------------------------------------

test_that("global_clustering_coefficient_wei", {
  expect_equal(
    global_clustering_coefficient_wei(W), 2.62438556)
})

# ---------------------------------------------------------------------------
# 7. Internal / public consistency
# ---------------------------------------------------------------------------

test_that("public characteristic_path_length equals internal _ version on valid matrix", {
  expect_equal(
    characteristic_path_length(W),
    characteristic_path_length_(W)
  )
})

test_that("public network_density equals internal _ version on valid matrix", {
  expect_equal(
    network_density(W),
    network_density_(W)
  )
})

test_that("public node_strength equals internal _ version on valid matrix", {
  expect_equal(
    node_strength(W),
    node_strength_(W)
  )
})

# ---------------------------------------------------------------------------
# 8. compute_global_metrics smoke test
# ---------------------------------------------------------------------------

test_that("compute_global_metrics returns a data.frame with expected columns", {
  arr <- array(W, dim = c(4, 4, 1))
  result <- compute_global_metrics(
    arr,
    c("characteristic_path_length", "network_density", "global_efficiency_wei"),
    subject_names = "S1"
  )
  expect_s3_class(result, "data.frame")
  expect_true("subject" %in% names(result))
  expect_true("characteristic_path_length" %in% names(result))
  expect_true("network_density" %in% names(result))
  expect_true("global_efficiency_wei" %in% names(result))
  expect_equal(nrow(result), 1L)
  expect_equal(result$subject, "S1")
})

# ---------------------------------------------------------------------------
# Betweenness centrality tests (fixed weighted_network reference)
# ---------------------------------------------------------------------------

test_that("betweenness_wei on fully connected 4-node graph", {
  expect_equal(betweenness_wei(fully_connected), c(2, 2, 2, 2))
})

test_that("betweenness_wei on 3-node path graph (weighted_network)", {
  # path: 1-2-3, node 2 lies on the only path between 1 and 3
  expect_equal(betweenness_wei(weighted_network), c(0, 0.5, 0))
})
