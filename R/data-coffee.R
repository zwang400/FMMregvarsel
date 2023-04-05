#' Coffee Data
#'
#' Data on the chemical composition of coffee samples collected from around the world, comprising
#' 43 samples from 29 countries. Each sample is either of the Arabica or Robusta variety. Twelve
#' of the thirteen chemical constituents reported in the study are given. The omitted variable is total
#' chlorogenic acid; it is generally the sum of the chlorogenic, neochlorogenic and isochlorogenic acid
#' values.
#'
#' @format A data frame with 43 observations and 14 variables. The first two columns contain Variety
#' and Country, respectively, while the remaining 12 columns contain the chemical properties. The
#' Variety is either (1) Arabica or (2) Robusta
#'
#' @note The German to English translations of the variable names were carried out by Dr. Sharon M. McNicholas.
#' @source Streuli, H. (1973). Der heutige stand der kaffeechemie. In Association Scientifique International
#' du Cafe, 6th International Colloquium on Coffee Chemisrty, Bogata, Columbia, pp. 61–72.
#'
#' @examples
#' coffee
#' variety_coffee <- coffee$Variety
#'
"coffee"
