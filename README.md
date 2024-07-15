# NetAnalyseR

An R package for producing graph theory metrics and other metrics from networks represented as adjacency matrices. Networks are processed in an automated pipeline format and are outputted as a .dataframe object to facilitate statistical analysis and graphing.
### Features 
- Streamlined adjacency matrix import with scan name parsing.
- Automated global and nodal network metric calculation.
- Automated statistical comparisons of network metrics between groups.

## Installation

The package devtools must first be installed:

```{R}
install.packages("devtools")
devtools::install_github("eclnz/NetAnalyseR")
library(NetAnalyseR)
```

## Usage
### Import adjacency matrices

```{r}
library(NetAnalyseR)

# Input parameters
data_dir <- system.file("extdata", package = "NetAnalyseR")
file_convention <- ".csv" 

# Import adjacency matrices
output <- process_matrices(
  directory = data_dir,
  file_convention = file_convention
  )

attributes(output)
```

### Calculate global metrics

```{r}
global_metrics <- c('global_efficiency_wei', 'global_clustering_coefficient_wei')
global_df <- compute_global_metrics(
  matrices_array = output$matrices, 
  global_metrics = global_metrics, 
  subject_names = output$subjects) %>% 
  # Allocate subjects to group
  allocate_groups(list(sub1 = 'sub1', sub2 = 'sub2'))
head(global_df)
```

### Statistical analysis of group differences

```{r}
global_df %>% group_comparisons(global_metrics, list(c('sub1','sub2')))
```

## Licence

This project is licensed under the GNU General Public License. See the LICENSE file for details.

## Credits

-   The logic from many functions were translated from code in the Brain Connectivity Toolbox (<https://sites.google.com/site/bctnet/>)
