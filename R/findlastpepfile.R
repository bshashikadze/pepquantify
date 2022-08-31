#' returns the location of the most recent peptides and modified conditions file
#'
#' @return last folder
#' @export
#' @import dplyr
#' @importFrom magrittr %>%
#' @examples findlastpepfile()
findlastpepfile <- function() {

  # which is the last folder
  mypath <- getwd()
  dir_info <- file.info(list.files(mypath, full.names = T)) %>%
    dplyr::filter(.data$isdir == "TRUE")
  last_folder <- rownames(dir_info)[which.max(dir_info$ctime)]

  # peptides and the labels file from that folder
  peptides               <- paste0(last_folder, "/peptidesformsempire.txt")
  labels                 <- paste0(last_folder, "/conditions_dual.txt")

  # save locations in the list to read it during msempire analysis
  locations_last <- list()
  locations_last[[1]] <- paste0(last_folder, "/peptidesformsempire.txt")
  locations_last[[2]] <- paste0(last_folder, "/conditions_dual.txt")
  locations_last[[3]] <- last_folder
  return(locations_last)
}
