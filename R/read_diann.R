#' Read and filter the DIA quantification data produced by the DIA-NN
#'
#' @param Q_Val refer to https://github.com/vdemichev/DiaNN
#' @param Global_Q_Val refer to https://github.com/vdemichev/DiaNN
#' @param Global_PG_Q_Val refer to https://github.com/vdemichev/DiaNN
#' @param Lib_Q_Val refer to https://github.com/vdemichev/DiaNN
#' @param Lib_PG_Q_Val refer to https://github.com/vdemichev/DiaNN
#' @param experimental_library set true if you use empirical libraries (e.g. prefractionation or GPF), false in case of lib free search with mbr enabled
#' @param unique_peptides_only TRUE only unique peptides will be used for quantification (recommended)
#' @param Quant_Qual  refer to https://github.com/vdemichev/DiaNN; pepquantify by default sets it to 0.5
#' @param id_column default "Genes"
#' @param exclude_samples if not empty, excludes specified sample/s from further analysis (only if necessary, e.g. after inspecting PCA)
#' @param quantity_column default "Genes.MaxLFQ.Unique", not important for MS-EmpiRe
#' @param sum_charge how precursor charge states will be aggregated to peptide level, True means the sum will be taken, in case of false, precursor with the highest intensity will be kept, default false
#' @param save_supplementary default TRUE, output is peptide and protein level data which can be used as a supplement
#' @import utils dplyr stringr tidyr
#' @importFrom magrittr %>%
#' @importFrom seqinr read.fasta getName
#'
#'
#'
#' @return The list of three elements, the first is the filtered peptides file, the second is the filtered protein groups file and the last is the character that stores the type of experiment
#' @export
#'
#' @examples read_diann(exclude_samples=c("samplename"), experimental_library = TRUE)


read_diann <- function(Q_Val = 0.01, Global_Q_Val = 0.01,
                                     Global_PG_Q_Val = 0.01, Lib_Q_Val = 0.01,
                                     Lib_PG_Q_Val = 0.01,
                                     experimental_library,
                                     unique_peptides_only = TRUE,
                                     Quant_Qual = 0.5,
                                     id_column = "Genes", quantity_column = "Genes.MaxLFQ.Unique",
                                     sum_charge = TRUE, save_supplementary = TRUE, exclude_samples=c()) {

  stopifnot("working directory does not contain .tsv file; make sure your working directory contains main output of the DIA-NN" =
              any(stringr::str_ends(list.files(getwd()), ".tsv")))



  # find the dia nn main output (naivly assumed that the largest file with extension .tsv is the one)
  raw_diann_path <- data.frame(list.files(getwd()), file.size(list.files(getwd()))) %>%
    dplyr::rename(files = 1, size = 2) %>%
    dplyr::filter(stringr::str_ends(.data$files, ".tsv") & .data$size == max(.data$size)) %>%
    dplyr::select(.data$files) %>%
    as.character()


  # read the DIA-NN output
  data <- read.delim(raw_diann_path, header = T, sep = "\t")


  # experimental library based analysis (e.g. GPF)
  if (experimental_library == T) {

    data <- data %>%
    dplyr::filter(.data$Q.Value           <= Q_Val) %>%
    dplyr::filter(.data$Global.Q.Value    <= Global_Q_Val) %>%
    dplyr::filter(.data$Global.PG.Q.Value <= Global_PG_Q_Val)

  }

  # library free analysis (with mbr enabled)
  else {

    data <- data %>%
      dplyr::filter(.data$Q.Value         <= Q_Val) %>%
      dplyr::filter(.data$Lib.Q.Value     <= Lib_Q_Val) %>%
      dplyr::filter(.data$Lib.PG.Q.Value  <= Lib_PG_Q_Val)


  }

  # filter for unique peptides and signal quality
  if (unique_peptides_only == TRUE) {

    unique <- 1

  }

  else {unique <- 0}

  data <- data %>%
    dplyr::filter(.data$Proteotypic      >= unique) %>%
    dplyr::filter(.data$Quantity.Quality >= Quant_Qual)


  # save filtered data
  write.table(data, "diannoutput_filtered.tsv", sep = "\t", row.names = F, quote = F)



  # subset for necessary columns
  data_filtered <- data %>%
    dplyr::select("Precursor.Quantity", "Precursor.Normalised", "Run", dplyr::all_of(quantity_column),
                  "Stripped.Sequence", dplyr::all_of(id_column))



  # aggregate charge states
  data_peptide <- data_filtered %>%
    dplyr::select(-dplyr::all_of(quantity_column))


  if (sum_charge == TRUE) {

    # summerise charge states by taking the sum
    data_peptide <- data_peptide %>%
      dplyr::group_by(.data$Run, .data$Stripped.Sequence, !!as.symbol(id_column)) %>%
      dplyr::summarise_all(sum, na.rm = TRUE) %>%
      dplyr::ungroup()

    }

  else {

    # aggregate charge states by taking the precursor with the highest intensity
    data_peptide <- data_peptide %>%
      dplyr::group_by(.data$Run, .data$Stripped.Sequence, !!as.symbol(id_column)) %>%
      dplyr::summarise(Precursor.Quantity   = max(.data$Precursor.Quantity),
                       Precursor.Normalised = max(.data$Precursor.Normalised)) %>%
      dplyr::ungroup()

  }

  # from long to wide also import additional columns from the original data
  data_peptide <- data_peptide %>%
    tidyr::pivot_wider(names_from = "Run",
                       values_from = c("Precursor.Quantity", "Precursor.Normalised"),
                       c(dplyr::all_of(id_column), "Stripped.Sequence"))



  data_peptide$peptide_Q_Val <- data$Global.Q.Value[match(data_peptide[[id_column]], data[[id_column]])]



  # re-order
  data_peptide <- data_peptide %>%
    dplyr::select(dplyr::all_of(id_column), .data$Stripped.Sequence,
                  starts_with("Precursor.Quantity"), starts_with("Precursor.Quantity"),
                  starts_with("Precursor.Normalised"), .data$peptide_Q_Val)


  # remove entries without gene names
  if ("Genes" %in% colnames(data_peptide)) {

    data_peptide      <- data_peptide[data_peptide$Genes != "", ]

  }


  # write results
  if (save_supplementary == TRUE) {

    write.table(data_peptide, "peptides.txt", row.names = F, sep = "\t")

  }


  # protein level data
  data_pg <- data_filtered %>%
    dplyr::select(.data$Run, dplyr::all_of(id_column), dplyr::all_of(quantity_column)) %>%
    dplyr::distinct(.data$Run, !!as.symbol(id_column), .keep_all = T) %>%
    dplyr::mutate(Run = stringr::str_c("LFQ.intensity_", .data$Run)) %>%
    tidyr::pivot_wider(names_from = "Run", values_from = dplyr::all_of(quantity_column), dplyr::all_of(id_column))


  # count number of peptides
  n_pep <- data_peptide %>%
    dplyr::group_by(!!as.symbol(id_column)) %>%
    dplyr::summarise(n_pep = dplyr::n_distinct(.data$Stripped.Sequence))


  # match other columns
  n_pep$pg_Q_Val                  <- data$Global.PG.Q.Value[match(n_pep[[id_column]], data[[id_column]])]


  # remove entries without gene names
  if ("Genes" %in% colnames(data_peptide)) {

    data_pg      <- data_pg[data_pg$Genes != "", ]

  }

  # join number of peptides and q-values
  data_pg <- data_pg %>%
    dplyr::left_join(n_pep)



  # write results
  if (save_supplementary == TRUE) {

    write.table(data_pg, "proteingroups.txt", row.names = F, sep = "\t")

  }


  # prepare final matrix for MS-EmpiRe
  data_peptide <- data_peptide %>%
  dplyr::select(dplyr::all_of(id_column), tidyr::starts_with("Precursor.Quantity")) %>%
  dplyr::rename(id=1) %>%
  dplyr::mutate_all(~replace(., is.na(.), 0)) %>%
  dplyr::mutate(unique_id = paste(.data$id, seq(1: nrow(data_peptide)), sep = "."), .keep="all") %>%
  dplyr::select(-dplyr::all_of(exclude_samples))


  # for convention Precursor.Normalised will be replaced by Intensity.
  names(data_peptide)<-gsub("Precursor.Quantity_", "Intensity.", colnames(data_peptide))


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
  filtered_data[[3]] <- "DIA"
  filtered_data[[4]] <- id_column

  # return the filtered peptides and protein groups as well as the information about the type of experiment
  return(filtered_data)

}


