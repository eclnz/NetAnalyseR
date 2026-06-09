# NetAnalyseR — Refactor TODO

## In progress

## Done
- [x] **1. Two-layer validation** — strip `validate` arg from public API; create unexported
  internal helpers that skip validation; enforce the boundary by function identity not a boolean flag
- [x] **2. Replace `get()` dispatch** — replace `get(metric_function)(...)` in
  `compute_global_metrics` and `compute_nodal_metrics` with explicit named function registries
- [x] **3. Parallelise subject loops** — wrap independent per-subject iterations in
  `compute_global_metrics` and `compute_nodal_metrics` with `parallel::mclapply` via `workers` param
- [x] **4. Expose magic numbers** — `0.35` edge-prevalence threshold in
  `z_structural_deviation` exposed as `edge_prevalence_threshold` parameter; `0.55` density cutoff
  in `shortest_distance_` documented with a comment explaining the algorithm crossover
- [x] **5. Test suite** — built out `tests/testthat/test-network_metrics.R` with regression tests
  for validate_matrix, network_density, node_strength, characteristic_path_length,
  global_efficiency_wei, global_clustering_coefficient_wei, internal/public consistency,
  compute_global_metrics smoke test, and fixed broken `weighted_network` reference
- [x] validate_matrix: square check before upper.tri
- [x] All metric functions: capture validate_matrix() return value
- [x] compare_rich_club: fix crash when plot=FALSE
- [x] z_structural_deviation: fix vacuous control-group guard
- [x] threshold_density: fix 2x density bug in local network_density closure
- [x] compute_global_metrics: fix small_worldness column name; fix array guard order
- [x] process_matrices: validate file_convention before use; fix O(n^2) vector growth
- [x] rich_club: remove debug global assignments; fix O(n) while loop
- [x] dijkstraAllPairs: DBL_MAX → infinity() for unreachable nodes
- [x] localEfficiencyWei: rewrite C++ with correct BCT formula; wire up R wrapper
