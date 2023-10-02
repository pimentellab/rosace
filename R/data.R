#' OCT1 drug cytotoxicity screen: raw count
#' 
#' A dataset containing three replicates of OCT1 drug cytotoxicity screen. 
#' Paper link: https://www.biorxiv.org/content/10.1101/2023.06.06.543963v1.full
#'  
#' @name oct1
#' @usage data(oct1)
#'
#' @format A data frame with variant as rows. The first column is the variant 
#' label, and the rest of them are raw counts at different time points in order. 
#' @keywords datasets
NULL

#' @rdname oct1
"oct1_rep1"
#' @rdname oct1
"oct1_rep2"
#' @rdname oct1
"oct1_rep3"

#' OCT1 drug cytotoxicity screen: processed rosace object
#' 
#' @name oct1_rosace
#' @usage data(oct1_rosace)
#'
#' @format A rosace object with three assays (replicates) and one assayset.  
#' @keywords datasets
"oct1_rosace"


