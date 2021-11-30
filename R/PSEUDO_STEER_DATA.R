#' Generate a pseudo dataset to reproduce results in the analysis study
#'
#' A pseudo dataset containing subect's id, group name, time points, and response variables.
#'
#'
#' @format A dataframe including 6 subjects with 4 variables has the repeated measured at 6 different time points for each subject :
#' \describe{
#'   \item{Steer}{identification number of a subject (steer).}
#'   \item{Group}{name of the group each subject belongs to}
#'     \item{Time}{6 time points at which the response variables are measured.}
#'   \item{citrulline}{response variables (amino acids) measured at 6 different time points from the same steer}
#' }
#'e
#' @source see data-raw/generate_data.R
#'
#'
#'
#'
#'
"pseudo_steer"

