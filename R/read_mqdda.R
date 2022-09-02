#' Read and filter the DDA quantification data produced by the MaxQuant
#' @description This functions expects that the directory contains proteinGroups.txt and peptides.txt produced by the MaxQuant without any post-analysis modification.
#'     Potential contaminants, reverse and only identified by site will be removed.
#'     Additionally, conditions file will be written in the directory which should be modified based on experimental conditions.
#'     Note that first you need to copy the conditions.txt file and rename as conditions_modified.txt and change only the second column. do not attempt to rename column names.
#' @param exclude_samples if not empty, excludes specified sample/s from further analysis (only if necessary, e.g. after inspecting PCA)
#' @param lfq if non-labelled data is loaded, lfq must be set to true if labelling was performed (e.g. TMT) lfq should be set to false. For TMT Reporter.intensity.corrected is taken for quantification
#' @return The list of three elements, the first is the filtered peptides file, the second is the filtered protein groups file and the last is the character that stores the type of experiment
#' @export
#' @import dplyr utils stringr
#' @importFrom magrittr %>%
#'
#' @examples read_mqdda(exclude_samples=c("samplename"), lfq = TRUE)
read_mqdda <- function(exclude_samples=c(), lfq = TRUE) {


  # check if the necessary files exist in the directory
  stopifnot("proteinGroups.txt and/or peptides.txt file does not exist in the directory" = file.exists("proteinGroups.txt") &
              file.exists("peptides.txt"))

  # read the files
  data_pg         <- read.delim("proteinGroups.txt", sep = "\t", header = T)
  data_peptide    <- read.delim("peptides.txt",      sep = "\t", header = T)

  # check if the required columns exist
  stopifnot("Protein.IDs column was not found in the protein groups file"     = "Protein.IDs" %in% colnames(data_pg))
  stopifnot("Leading.razor.protein column was not found in the peptides file" = "Leading.razor.protein" %in% colnames(data_peptide))


  # filter by site, reverse and contaminants in the protein groups file
  data_pg <- data_pg %>%
    dplyr::filter(.data$Potential.contaminant   != "+") %>%
    dplyr::filter(.data$Only.identified.by.site != "+") %>%
    dplyr::filter(.data$Reverse != "+")

  # add only identified by site to peptides from the protein groups file (as they are missing there)
  data_peptide$Only.identified.by.site = data_pg$Only.identified.by.site[charmatch(data_peptide$Leading.razor.protein,
                                                                                   data_pg$Protein.IDs)]

  # filter by site, reverse and contaminants, in the peptides file
  data_peptide <- data_peptide %>%
    dplyr::filter(.data$Potential.contaminant  != "+") %>%
    dplyr::filter(.data$Reverse != "+") %>%
    dplyr::filter(.data$Only.identified.by.site != "+")


  # add protein grouping column to the peptide (from the protein groups file)
  data_peptide$id = data_pg$Protein.IDs[charmatch(data_peptide$Leading.razor.protein,
                                                  data_pg$Protein.IDs)]


  # distinct peptides from the same proteins will have the same ids
  # convert id to an unique id by adding .number (1-nrow) (sequence)
  # depending if its labelled or lfq data different columns will be selected

  if (lfq == TRUE) {

  stopifnot("columns starting with Intensity. do not exist, did you set lfq to true for the labelled data?" = any(grepl("Intensity.", colnames(data_peptide))))

     data_peptide <- data_peptide %>%
      dplyr::mutate(unique_id = paste(id, seq(1: nrow(data_peptide)), sep = "."), .keep= "all") %>%
      dplyr::select(.data$id, .data$unique_id, starts_with("Intensity.")) %>%
      dplyr::select(-all_of(exclude_samples))

  }


  else {

  stopifnot("columns starting with Reporter.intensity.corrected. do not exist, did you set lfq to false for the lfq data?" = any(grepl("Reporter.intensity.corrected", colnames(data_peptide))))

     data_peptide <- data_peptide %>%
      dplyr::mutate(unique_id = paste(id, seq(1: nrow(data_peptide)), sep = "."), .keep="all") %>%
      dplyr::select(.data$id, .data$unique_id, starts_with("Reporter.intensity.corrected.")) %>%
      dplyr::select(-all_of(exclude_samples))

  # for convention Reporter.intensity.corrected will be replaced by Intensity.
      names(data_peptide)<-gsub("Reporter.intensity.corrected.", "Intensity.", colnames(data_peptide))

  }


  # write conditition file
  Bioreplicate <- colnames(data_peptide)[! colnames(data_peptide) %in% c('id', 'unique_id')]
  Condition    <- stringr::str_remove(Bioreplicate, "Intensity.")
  Conditions   <- data.frame(Bioreplicate, Condition) %>%
    write.table("conditions.txt", row.names = F, sep = "\t")


  # print information
  cat("conditions file was generated. First rename file as: conditions_modified.txt. afterwards,
         modify ONLY the second column according to the experimental conditions. Do not change the column headers")


  # store the results in the list
  filtered_data      <- list()
  filtered_data[[1]] <- data_peptide
  filtered_data[[2]] <- data_pg
  filtered_data[[3]] <- "DDA"


  # return the filtered peptides and protein groups as well as the information about the type of experiment
  return(filtered_data)
}
