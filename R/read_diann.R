#' Read and filter DIA precursor-level data produced by the DIA-NN
#'
#' @param Q_Val refer to https://github.com/vdemichev/DiaNN
#' @param Global_Q_Val refer to https://github.com/vdemichev/DiaNN
#' @param Global_PG_Q_Val refer to https://github.com/vdemichev/DiaNN
#' @param Lib_Q_Val refer to https://github.com/vdemichev/DiaNN
#' @param Lib_PG_Q_Val refer to https://github.com/vdemichev/DiaNN
#' @param experimental_library set true if you use empirical libraries (e.g. prefractionation or GPF), false in case of lib free search with mbr enabled
#' @param unique_peptides_only TRUE only unique peptides will be used for quantification (recommended)
#' @param Quant_Qual  refer to https://github.com/vdemichev/DiaNN; pepquantify by default sets it to 0.5
#' @param id_column default "Genes" (can be switched to Protein.Group)
#' @param second_id_column default "Protein.Group" (can be switched to Genes)
#' @param exclude_samples if not empty, excludes specified sample/s from further analysis (only if necessary, e.g. after inspecting PCA)
#' @param quantity_column default "Genes.MaxLFQ.Unique", not important for MS-EmpiRe
#' @param diann_file_name by default function will assume that largest .tsv file is the DIA-output, however specific file name can be specified diann_file_name=filename.tsv
#' @param directory default is a current directory
#' @param sum_charge how precursor charge states will be aggregated to peptide level, True means the sum will be taken, in case of false, precursor with the highest intensity will be kept, default false
#' @param save_supplementary default TRUE, output is peptide and protein level data which can be used as a supplement
#' @param for_msempire default TRUE so pepquantify will prepare data for MS-EmpiRe quantification, if purpose is to filter DIA-NN output and generate peptides and protein groups file set to false
#' @param include_mod_in_pepreport default FALSE, if true includes modifications in the output peptide file (currently only Carbamidomethyl (C))
#' @import utils dplyr stringr tidyr
#' @importFrom magrittr %>%
#'
#'
#'
#' @return The list of three elements, the first is the filtered peptides file, the second is the filtered protein groups file and the last is the character that stores the type of experiment
#' @export
#'
#' @examples read_diann(exclude_samples=c("samplename"), experimental_library = TRUE)


read_diann <- function(Q_Val = 0.01,
                       Global_Q_Val             = 0.01,
                       Global_PG_Q_Val          = 0.01,
                       Lib_Q_Val                = 0.01,
                       Lib_PG_Q_Val             = 0.01,
                       experimental_library,
                       unique_peptides_only     = TRUE,
                       Quant_Qual               = 0.5,
                       diann_file_name          = "",
                       directory                = ".",
                       id_column                = "Genes",
                       second_id_column         ="Protein.Group",
                       quantity_column          = "Genes.MaxLFQ.Unique",
                       for_msempire             = T,
                       sum_charge               = TRUE,
                       save_supplementary       = TRUE,
                       exclude_samples          = c(),
                       include_mod_in_pepreport = F) {


  # set working directory
  setwd(directory)
  print(paste0("R is located in ", getwd(), ", if this path is wrong, change it from R studio, or by specifying the correct path using the directory command e.g. directory = name_of_the_path"))


  # read main output of DIA-NN (if dia_nn_file name is not specified largest .tsv file will be loaded)
  if (!grepl(pattern = '.tsv|.txt|.csv', diann_file_name)) {

    stopifnot("pepquantify attempted to automatically read DIA-NN main output but working directory does not seem to contain .tsv file; make sure working directory (printed above) is correct;
              also consider explicit reading using: diann_file_name = filename.extension" =
                any(stringr::str_ends(list.files(getwd()), ".tsv")))


    # find the dia nn main output (naivly assumed that the largest file with extension .tsv is the one)
    raw_diann_path <- data.frame(list.files(getwd()), file.size(list.files(getwd()))) %>%
      dplyr::rename(files = 1, size = 2) %>%
      dplyr::filter(stringr::str_ends(files, ".tsv") & size == max(size, na.rm = T)) %>%
      dplyr::select(files) %>%
      as.character()


    # read the DIA-NN output
    print(paste0("file, with the name ", raw_diann_path, " is currently loading"))
    data <- read.delim(raw_diann_path, header = T, sep = "\t")


  }

  if (grepl(pattern = '.tsv|.txt|.csv', diann_file_name)) {

    # stop if specified file does not exist in the directory
    stopifnot("such file does not exist in the working directory, make sure to type the exact file name with an extension e.g. diann_file_name=myfilename.tsv" =
                any(stringr::str_detect(list.files(getwd()), diann_file_name)))


    # read the DIA-NN output
    print(paste0("file, with the name ", diann_file_name, " is currently loading"))
    data <- read.delim(diann_file_name, header = T, sep = "\t")


  }


  # experimental library based analysis (e.g. GPF)
  if (experimental_library == T) {

    data <- data %>%
    dplyr::filter(Q.Value           <= Q_Val) %>%
    dplyr::filter(Global.Q.Value    <= Global_Q_Val) %>%
    dplyr::filter(Global.PG.Q.Value <= Global_PG_Q_Val)

  }

  # library free analysis (with mbr enabled)
  if (experimental_library == F) {

    data <- data %>%
      dplyr::filter(Q.Value         <= Q_Val) %>%
      dplyr::filter(Lib.Q.Value     <= Lib_Q_Val) %>%
      dplyr::filter(Lib.PG.Q.Value  <= Lib_PG_Q_Val)


  }

  # filter for unique peptides and signal quality
  if (unique_peptides_only == TRUE) {

    unique <- 1

  }

  else {unique <- 0}

  data <- data %>%
    dplyr::filter(Proteotypic      >= unique) %>%
    dplyr::filter(Quantity.Quality >= Quant_Qual)


  # create folder where outputs will be saved
  if (exists("data", inherits = F) & !dir.exists("txt")) {

    dir.create("txt")

  }

  # save filtered data
  write.table(data, "txt/diannoutput_filtered.tsv", sep = "\t", row.names = F, quote = F)


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
      dplyr::group_by(Run, Stripped.Sequence, !!as.symbol(id_column)) %>%
      dplyr::summarise_all(sum, na.rm = TRUE) %>%
      dplyr::ungroup()

    }

  else {

  # aggregate charge states by taking the precursor with the highest intensity
  data_peptide <- data_peptide %>%
      dplyr::group_by(Run, Stripped.Sequence, !!as.symbol(id_column)) %>%
      dplyr::summarise(Precursor.Quantity   = max(Precursor.Quantity),
                       Precursor.Normalised = max(Precursor.Normalised)) %>%
      dplyr::ungroup()

  }


  # from long to wide also import additional columns from the original data
  data_peptide <- data_peptide %>%
    tidyr::pivot_wider(names_from = "Run",
                       values_from = c("Precursor.Quantity", "Precursor.Normalised"),
                       c(dplyr::all_of(id_column), "Stripped.Sequence"))


  # match other data (charge, modifications, q-values)
  if (include_mod_in_pepreport == T) {

    print("in the peptide output modifications will be included (so far (v.2.1.2) only works for Carbamidomethyl(C))")

    precursors_data_modification <- data %>%
      dplyr::select(Modified.Sequence, Stripped.Sequence) %>%
      dplyr::mutate(Modification = str_extract_all(Modified.Sequence,"(?<=\\().+(?=\\))")) %>%
      dplyr::mutate(Modification = case_when(Modification == "UniMod:4" ~ "Carbamidomethyl (C)", TRUE ~ "")) %>%
      dplyr::distinct(Stripped.Sequence, Modification, .keep_all = T) %>%
      dplyr::group_by(Stripped.Sequence) %>%
      dplyr::summarize(Modification=paste(Modification,collapse=";")) %>%
      dplyr::ungroup()

  }

  # get unique protein descriptions (distinct protein names for the same id_column will be aggregated in one row separated by semicolon)
  protein_description <- data %>%
    dplyr::select(!!as.symbol(id_column), First.Protein.Description) %>%
    dplyr::distinct(!!as.symbol(id_column), First.Protein.Description, .keep_all = T) %>%
    dplyr::group_by(!!as.symbol(id_column)) %>%
    dplyr::summarize(First.Protein.Description=paste(First.Protein.Description,collapse=";")) %>%
    dplyr::ungroup()


  # get unique second_id
  second_id <- data %>%
    dplyr::select(!!as.symbol(id_column), !!as.symbol(second_id_column)) %>%
    dplyr::distinct(!!as.symbol(id_column), !!as.symbol(second_id_column), .keep_all = T) %>%
    group_by(!!as.symbol(id_column)) %>%
    dplyr::summarize(second_ids = paste(!!as.symbol(second_id_column), collapse=",")) %>%
    dplyr::ungroup()


  # combine
  combined_additional <- protein_description %>%
    dplyr::left_join(second_id)


  # q values separately for each charge state (experimental library based)
  precursors_data_qvalue <- data %>%
    dplyr::select(Modified.Sequence, Stripped.Sequence, Global.Q.Value) %>%
    dplyr::distinct(Stripped.Sequence, Global.Q.Value, .keep_all = T) %>%
    dplyr::group_by(Stripped.Sequence) %>%
    dplyr::summarize(Global.Q.Value=paste(Global.Q.Value,collapse=";")) %>%
    dplyr::rename(Q.value = Global.Q.Value) %>%
    dplyr::ungroup()


  # q values separately for each charge state (lib free mbr enabled)
  precursors_data_qvalue_mbr <- data %>%
    dplyr::select(Modified.Sequence, Stripped.Sequence, Lib.Q.Value) %>%
    dplyr::distinct(Stripped.Sequence, Lib.Q.Value, .keep_all = T) %>%
    dplyr::group_by(Stripped.Sequence) %>%
    dplyr::summarize(Lib.Q.Value=paste(Lib.Q.Value,collapse=";")) %>%
    dplyr::rename(Q.value = Lib.Q.Value) %>%
    dplyr::ungroup()


  # charges
  precursors_data_charge <- data %>%
    dplyr::select(Modified.Sequence, Stripped.Sequence, Precursor.Charge) %>%
    dplyr::mutate(Precursor.Charge = str_c("+", Precursor.Charge)) %>%
    dplyr::distinct(Stripped.Sequence, Precursor.Charge, .keep_all = T) %>%
    dplyr::group_by(Stripped.Sequence) %>%
    dplyr::summarize(Precursor.Charge=paste(Precursor.Charge,collapse=";")) %>%
    dplyr::rename(Charge = Precursor.Charge) %>%
    dplyr::ungroup()


  # include q values
  if (experimental_library == T) {

    prec_uggregated <- precursors_data_qvalue %>%
      dplyr::left_join(precursors_data_charge)

  }


  if (experimental_library == F) {

    prec_uggregated <- precursors_data_qvalue_mbr %>%
      dplyr::left_join(precursors_data_charge)

  }


  # include modification
  if (include_mod_in_pepreport == T) {

   prec_uggregated <- prec_uggregated %>%
      dplyr::left_join(precursors_data_modification)

      }

  else {

      prec_uggregated <- prec_uggregated
  }


  # add data to peptides
  data_peptide <- data_peptide %>%
    dplyr::left_join(prec_uggregated) %>%
    dplyr::left_join(data %>% dplyr::select(Stripped.Sequence, First.Protein.Description, !!as.symbol(second_id_column)) %>%
                       dplyr::distinct(Stripped.Sequence, !!as.symbol(second_id_column), .keep_all = T))



  # remove entries without gene names
  if ("Genes" == colnames(data_peptide)[1]) {

    data_peptide      <- data_peptide[data_peptide$Genes != "", ]

  }


  # write results
  if (save_supplementary == TRUE) {

    write.table(data_peptide, "txt/peptides.txt", row.names = F, sep = "\t", quote = F)

  }


  # protein level data
  data_pg <- data_filtered %>%
    dplyr::select(Run, dplyr::all_of(id_column), dplyr::all_of(quantity_column)) %>%
    dplyr::distinct(Run, !!as.symbol(id_column), .keep_all = T) %>%
    dplyr::mutate(Run = stringr::str_c("LFQ.intensity_", Run)) %>%
    tidyr::pivot_wider(names_from = "Run", values_from = dplyr::all_of(quantity_column), dplyr::all_of(id_column))


  # count number of peptides
  n_pep <- data_peptide %>%
    dplyr::group_by(!!as.symbol(id_column)) %>%
    dplyr::summarise(n_pep = dplyr::n_distinct(Stripped.Sequence))


  # match q values
  if (experimental_library == T) {

    n_pep$pg_Q_Val    <- data$Global.PG.Q.Value[match(n_pep[[id_column]], data[[id_column]])]

  }


  if (experimental_library == F) {

    n_pep$pg_Q_Val    <- data$Lib.PG.Q.Value[match(n_pep[[id_column]], data[[id_column]])]

  }



  # remove entries without gene names
  if ("Genes" == colnames(data_pg)[1]) {

    data_pg      <- data_pg[data_pg$Genes != "", ]

  }


  # join number of peptides and q-values
  data_pg <- data_pg %>%
    dplyr::left_join(n_pep) %>%
    dplyr::left_join(combined_additional)



  # write results
  if (save_supplementary == TRUE) {

    write.table(data_pg, "txt/proteingroups.txt", row.names = F, sep = "\t", quote = F)

  }


  # print message what has been saved
  if (save_supplementary == TRUE) {

    cat("the following files were saved in the txt folder:
    1 - diann_output_filtered;
    2 - peptides;
    3 - proteingroups.")

  }

  if (save_supplementary == FALSE) {

    cat("the following file were saved in the txt folder:
    diann_output_filtered.")

    }


  if (for_msempire == T) {

    # prepare final matrix for MS-EmpiRe
    data_peptide <- data_peptide %>%
      dplyr::select(dplyr::all_of(id_column), tidyr::starts_with("Precursor.Quantity")) %>%
      dplyr::rename(id=1) %>%
      dplyr::mutate_all(~replace(., is.na(.), 0)) %>%
      dplyr::mutate(unique_id = paste(id, seq(1: nrow(data_peptide)), sep = "."), .keep="all") %>%
      dplyr::rename_all(~str_replace(.,"Precursor.Quantity_", "Intensity.")) %>%
      dplyr::select(-dplyr::all_of(exclude_samples))


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

  if (for_msempire == F) {

    filtered_data      <- list()
    filtered_data[[1]] <- data_peptide
    filtered_data[[2]] <- data_pg
    return(filtered_data)
  }


}


