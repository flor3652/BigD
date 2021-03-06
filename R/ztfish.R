#' The zero truncated number of fish caught
#' 
#' This datset was pulled offline (see site below), and contains the amount of fish caught by certain people. The meanings of the columns are made up by me to give it context (though I think they are reasonable interpretations of the column headings). The zero counts were removed from this dataset.
#' 
#' @format A data frame with 108 rows and 8 columns
#' \describe{
#'   \item{nofish}{whether they went fishing, binary (1 meaning no fishing)}
#'   \item{livebait}{whether live bait was used, binary (1 meaning live bait was used)}
#'   \item{camper}{whether a camper was used, binary (1 meaning camper was used)}
#'   \item{persons}{how many persons attended, count}
#'   \item{child}{how many children attended, count}
#'   \item{xb}{unknown, apparently continuous}
#'   \item{zg}{unknown, apparently conitnuous}
#'   \item{count}{the amount of fish caught, count}
#' }
#' 
#' @source \url{http://www.ats.ucla.edu/stat/data/fish.csv}
"ztfish"