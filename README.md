# NetAnalyseR

An R package for producing graph theory metrics and other metrics from networks represented as adjacency matrices. Networks are processed in an automated pipeline format and are outputted as a .dataframe object to facilitate statistical analysis and graphing.

## Table of Contents

-   [Installation](#installation)
-   [Usage](#usage)
-   [Features](#features)
-   [Contributing](#contributing)
-   [License](#licence)
-   [Credits](#credits)

## Installation

The package devtools must first be installed:

```{R}
install.packages("devtools")
devtools::install_github("eclnz/NetAnalyseR")
library(NetAnalyseR)
```

## Usage

```{R}
library(NetAnalyseR)
library(magrittr)
library(ggplot2)

# Set directory where adjacencey matrices are located. Location of exemplar matrices below. 
data_dir <- system.file("extdata", package = "NetAnalyseR")

# Set names of all subjects you wish to analyse. 
subjects <- c("A", "B", "C", "D")

# The file naming convention of these files. This is any characters following the name. 
file_convention <- ".csv"

# Initial import of matrices into matrix array and edge data frame.
output <- process_matrices(data_dir, subjects, file_convention)

# Specify network metrics of interest.
global_metrics <- c("characteristic_path_length","network_density")
nodal_metrics <- c("node_strength", "local_efficiency_wei")

# Compute specified global metrics and allocate groups
global_df <- compute_global_metrics(output$matrices, global_metrics, output$subjects) %>% 
  allocate_groups(list(control = c("A","B"), case = c("C", "D")))

# Compute specified nodal metrics and allocate groups
nodal_df <- compute_nodal_metrics(output$matrices, nodal_metrics, output$subjects) %>% 
  allocate_groups(list(control = c("A","B"), case = c("C", "D")))

# Allocate edges to groups
edge_df <- classify_connections(output$edge_df, cortical_allocation = c(1,2), subcortical_allocation = c(3,4))

# Example plots
nodal_df %>% 
  ggplot(aes(x = group, y = local_efficiency_wei))+geom_point()
edge_df %>% 
  filter(self_connectivity=="Network-Connected") %>% 
  ggplot(aes(x = connection_type, y = edge_strength))+geom_point()

# Example statistical analysis
lm <- lm(characteristic_path_length ~ group, data = global_df )
summary(lm)

edge_df %>% 
  filter(self_connectivity=="Network-Connected") %>% 
  ggplot(aes(x = connection_type, y = edge_strength))+geom_point()

```

## Features 

-   Import adjacency matrices

## Contributing

Add the package folder

```{bash}
cd YOUR_DIRECTORY
git clone https://github.com/eclnz/NetAnalyseR.git
```

Open the file NetAnalyseR.Rproj in RStudio

```{R}
library(devtools)
document()
```

## Licence

This project is licensed under the GNU General Public License. See the LICENSE file for details.

## Credits

-   The logic from many functions were translated from code in the Brain Connectivity Toolbox (<https://sites.google.com/site/bctnet/>)
