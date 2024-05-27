##################
# STATS2TABLE.R #
##################
# Written by Olga Veth - s2067579 - University of Leiden
# Created on 30-09-2019
# Edited by Lara M Wierenga on 22-02-21
# Version 4.0
# Edited for FreeSurfer7
# From the Qoala-T QC repo
#   https://github.com/Qoala-T/QC/blob/master/Scripts/Stats2Table/Stats2Table_fs7.R
# Edited by Ed Clarkson - University of Auckland - 25-04-2024
##################
# BSD 3-Clause License
#
# Copyright (c) 2017-2019,
# Developed and created by Lara Wierenga and the Qoala-T team
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#   1. Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#          SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##################
#' Process MRI data from FreeSurfer output files
#'
#' @description This function reads, processes, and returns volumetric data from FreeSurfer output files for all subjects in a specified directory.
#' It reads volumetric and structural data from aseg.stats files, metadata from these files, area and thickness measurements from aparc.stats files,
#' and combines all these into a single dataset.
#' @param datasetDir The directory containing the data files for all subjects.
#' @param subject_names  A vector specifying the subject names to be analysed within the `datasetDir`.
#' @return Returns a dataframe containing FreeSurfer brain volumes
#' @importFrom purrr map
#' @importFrom purrr list_rbind
#' @importFrom dplyr select
#' @importFrom utils read.table
#' @importFrom tidyselect everything

processFSdata <- function(datasetDir, subject_names = NULL) {

  # Helper function to read aseg.stats files
  readAseg <- function(subjectDir) {
    aseg_file <- data.frame(utils::read.table(file.path(subjectDir,"stats","aseg.stats"), row.names=1))[,c(3,4)]
    t(data.frame(aseg_file[,1], row.names = aseg_file[,2]))
  }

  # Helper function to read metadata from aseg.stats
  readMetaAseg <- function(subjectDir) {
    aseg_meta <- readLines(file.path(subjectDir,"stats","aseg.stats"), n=35)[15:34]
    meta1 <- gsub("# ", "", aseg_meta)
    meta <- t(data.frame(strsplit(meta1, ",")))[,c(2,4)]
    metaTable <- t(data.frame(meta[,2]))
    colnames(metaTable) <- meta[,1]
    return(metaTable)
  }

  # Helper function to read and process aparc.stats files
  readAparc <- function(subjectDir, value) {
    sides <- c("lh", "rh")
    aparcTable <- NULL

    for (side in sides) {
      file_path <- file.path(subjectDir, "stats", paste0(side, ".aparc.stats"))
      areaThickness <- as.data.frame(read.table(file_path, row.names=1))[, c(2,4)]
      rowValues <- rownames(areaThickness)

      meta <- readLines(file_path)[c(20, 21)]
      meta1 <- gsub("# ", "", meta)
      meta2 <- t(data.frame(strsplit(meta1, ",")))[, c(2,4)]
      meta3 <- data.frame(meta2[ifelse(value == "area", 1, 2),2])
      value2 <- gsub(" ", "", meta2[ifelse(value == "area", 1, 2),1])

      colnames(meta3) <- paste(side, "_", value2, "_", value, sep="")
      extra <- t(matrix(areaThickness[,ifelse(value == "area", 1, 2)]))
      colnames(extra) <- paste(side, "_", rowValues, "_", value, sep="")
      aparcTable <- ifelse(is.null(aparcTable), cbind(extra, meta3), cbind(aparcTable, extra, meta3))
    }
    return(aparcTable)
  }

  # Read and merge all data
  readFiles <- function(subjectDir) {
    asegTable <- readAseg(subjectDir)
    metaTable <- readMetaAseg(subjectDir)
    areaAparc <- readAparc(subjectDir, "area")
    thickAparc <- readAparc(subjectDir, "thickness")
    cbind(asegTable, metaTable, areaAparc, thickAparc)
  }

  subjectDirs <- unique(list.dirs(datasetDir, recursive=FALSE))
  if(!is.null(subject_names)){
    include_index <- basename(subjectDirs) %in% subject_names
    subjectDirs <- subjectDirs[include_index]
  }
  statsExistIndex <- file.exists(file.path(subjectDirs, "stats", "aseg.stats"))
  warning("The following subjects are missing stats files: ", paste(basename(subjectDirs[!statsExistIndex]), collapse =", "),call. = FALSE)
  subjectDirs <- subjectDirs[statsExistIndex]
  # Function to read and process individual subject directories
  read_subject_data <- function(subjectDir) {
    cat("\t", subjectDir, "\n")
    stats_file <- file.path(subjectDir, "stats", "aseg.stats")

    if (file.exists(stats_file)) {
      subjectTable <- readFiles(subjectDir) %>% as.data.frame()  # Reading the data
      subjectTable$subject <- basename(subjectDir)  # Add subject identifier
      subjectTable <- subjectTable %>% dplyr::select(subject, tidyselect::everything())
      return(subjectTable)
    } else {
      warning(paste("No volumetric stats for", subjectDir), call. = FALSE)
    }
  }
  cat("Reading FS Volumes from:\n")
  stats2Table2 <- purrr::map(subjectDirs, read_subject_data) %>%
    purrr::list_rbind()

  # Replace "-" with "_" in column names
  names(stats2Table2) <- gsub("-", "_", names(stats2Table2))
  names(stats2Table2) <- gsub(" ", "", names(stats2Table2))
  rownames(stats2Table2) <- NULL

  # Make sure all volume entries in dataframe are numeric
  stats2Table2[,-1]  <- lapply(stats2Table2[,-1], as.numeric)

  return(stats2Table2)
}


