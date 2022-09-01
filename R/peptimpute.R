#' Performs missing value imputation (Perseus-like approach)
#'
#' @param data this is the data that is returned from the preppetide function
#' @param downshift see the Perseus documentation "Replace missing values from normal distribution"
#' @param width     see the Perseus documentation "Replace missing values from normal distribution"
#' @param n_ko_like minimum number of peptides that should have missing and valid value pattern (all valid in one condition, less than 2 in the second or otherwise by user defined criteria) default 2
#' @param fraction_valid between 0-1. 1 means that imputed peptides are taken into account if they are present in all samples of one of the conditions, 0.5 means if they are present in the half of the samples of one of the conditions. default 1
#' @param second_condition maximum acceptable number of valid values in other condition when fraction valid is met in the other, default 1
#' @param seed set seed as values for imputation are derived randomly, seed makes sure the reproducibility. default 1234
#' @description By default MS-EmpiRe does not require imputation, but requires at least two valid values in each condition.
#'    In extreme cases it could be that the peptide is consistently detected in one of the conditions but is not found in the second condition,
#'    such peptides might be interesting to futher explore or even visualize on the volcano plot. This function founds such peptides and imputes missing values. However after quantitative analysis generated p-values and fold-changes must be taken into
#'    account with caution. This function tries to take those imputed peptides for quantification that are chosen with very stringend criterias. Additionally, if protein already have
#'    peptides more than min_pep (see the preppetide function), no imputation will be performed even if there are peptides for this protein that meet imputation criteria.
#'    For the criterias see the options of this function.
#' @import utils dplyr tibble tidyr
#' @importFrom stats sd median rnorm
#' @importFrom magrittr %>%
#' @export
#'
#' @examples peptimpute(data)
peptimpute <- function(data, downshift = 1.8, width = 0.3, n_ko_like = 2, fraction_valid = 1,
                       second_condition = 1, seed = 1234) {


  # get the statistics measured from the valid data. this will be used to build a random distribution from which the values will be
  #drawn
  valid_data_descriptives <- data[[2]] %>%
    rownames_to_column("unique_id") %>%
    tidyr::pivot_longer(names_to = "Bioreplicate", values_to = "Intensity", -.data$unique_id) %>%
    dplyr::group_by(.data$Bioreplicate) %>%
    dplyr::summarise(mean_valid = mean(.data$Intensity, na.rm = T),
              median_valid = median(.data$Intensity, na.rm = T),
              sd_valid   = sd(.data$Intensity,   na.rm = T),
              n_valid    = sum(!is.na(.data$Intensity)),
              n_missing  = nrow(data[[1]]) - .data$n_valid) %>%
    ungroup()


  # imputation
  column_indices <- 1:nrow(valid_data_descriptives)
  random_data <- list()
  # makes the list which contains as many elements as the samples are and as many random values as the missing values are in each
  #sample
  for (i in column_indices) {

    set.seed(seed = seed)
    random_data[[i]] <- rnorm(n = valid_data_descriptives$n_missing[i],
                              mean = valid_data_descriptives$median_valid[i] - (downshift * valid_data_descriptives$sd_valid[i]),
                              sd = valid_data_descriptives$sd_valid[i]   * width)

    # impute the missing values
    data[[2]][is.na(data[[2]])[, valid_data_descriptives$Bioreplicate[i]],
              valid_data_descriptives$Bioreplicate[i]] <- random_data[[i]]

  }


  # get once more information about the conditions
  conditions    <- read.delim("conditions_modified.txt")
  condition1    <- data[[4]][1]
  condition2    <- data[[4]][2]
  length_cond_1 <- length(conditions[conditions == data[[4]][1]])
  length_cond_2 <- length(conditions[conditions == data[[4]][2]])



  # is the fraction valid between 0 and 1?
  stopifnot("fraction_valid should be between 0 and 1" = fraction_valid <= 1 & fraction_valid > 0)


  # stop if fraction valid is so low that it includes peptides that are already quantifiable by msempire
  if (fraction_valid * length_cond_1 < data[[5]][1] | fraction_valid * length_cond_2 < data[[5]][2]) {

    fraction_valid <- 1

    cat(" fraction valid was low enough to include peptides that are already quantifiable by msempire, therefore it was set to 1.
            continue with this or increase the value until you do not get this warning ")
  }

  else {fraction_valid <- fraction_valid}



  # imputed matrix. imputation was performed on entire dataset. But finally only imputed peptides which were present in all samples
  # of any condition and maximum two samples of other condition will be left for the msempire
  imputed_matrix <- data[[2]] %>%
    tibble::rownames_to_column("unique_id") %>%
    dplyr::left_join(data[[1]]) %>%
    dplyr::filter(!!as.symbol(condition1)  >= fraction_valid * length_cond_1  & !!as.symbol(condition2)  <= second_condition
           |!!as.symbol(condition1) <=  second_condition & !! as.symbol(condition2) >= fraction_valid * length_cond_2) %>%
    dplyr::select(.data$id, .data$unique_id, starts_with("Intensity.")) %>%
    dplyr::mutate_if(is.numeric, ~2^(.)) %>%
    dplyr::group_by(.data$id) %>%
    dplyr::mutate(n_ko = dplyr::n()) %>%
    dplyr::filter(.data$n_ko >= n_ko_like) %>%
    dplyr::ungroup()



  # this is unimputed matrix which contains potentially msempire quantifiable peptides
  peptides_for_msempire_min_2_pep <- data[[3]] %>%
    dplyr::group_by(.data$id) %>%
    dplyr::mutate(n_pep = dplyr::n()) %>%
    dplyr::filter(.data$n_pep >= data[[6]]) %>%
    dplyr::ungroup()



  # detect peptides in imputed matrix that can anyway be quantified by msempire so are not necessary to include
  imputed_matrix <- imputed_matrix %>%
    dplyr::filter(!id %in% peptides_for_msempire_min_2_pep$id)


  # print information about imputed peptides
  if (nrow(imputed_matrix) > 0) {
    cat(paste0(" there are ", nrow(imputed_matrix), " peptides with the missing values imputed, they will be included in the quantification "))
  }

  else {cat(" no peptides have met to the imputation criteria. you can relax the parameters or continue without the imputation ")}


  # add imputed peptides (if any)
    peptides_for_msempire <- data[[3]] %>%
      dplyr::bind_rows(imputed_matrix) %>%
      dplyr::group_by(.data$id) %>%
      dplyr::mutate(n_pep = dplyr::n()) %>%
      dplyr::filter(.data$n_pep >= data[[6]]) %>%
      dplyr::ungroup() %>%
      dplyr::select(-.data$id) %>%
      write.table(paste0(condition1, "_vs_", condition2, "/peptidesformsempire.txt"), sep = "\t", row.names = F)

}
