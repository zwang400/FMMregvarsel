#' Happy Score Data
#'
#' Happy scores and six life quality factors for 149 countries in 2021.
#'
#' @format A data frame with 149 observations and 8 variables:
#' \describe{
#' \item{country}{country name}
#' \item{Happy.Score}{happy score for countries in 2021, it is the national average response to the question of life evaluations.}
#' \item{GDP.per.cap}{The GDP per capita. Calculated as the Total GDP divided by the total population.}
#' \item{Socl.spprt}{Social support (or having someone to count on in times of trouble) is the national
#' average of the binary responses (either 0 or 1) to the GWP question “If you
#' were in trouble, do you have relatives or friends you can count on to help you
#' whenever you need them, or not?”.}
#' \item{Life.expt}{Healthy Life Expectancy (HLE), based on the data extracted from the World Health Organization’s (WHO) Global Health Observatory data repository.}
#' \item{Freedom}{Freedom to make life choices is the national average of responses to the GWP
#' question “Are you satisfied or dissatisfied with your freedom to choose what
#' you do with your life?”.}
#' \item{Generosity}{Generosity is the residual of regressing national average of response to the GWP
#' question “Have you donated money to a charity in the past month?” on GDP
#' per capita.}
#' \item{Corruption}{Corruption Perception. The measure is the national average of the survey responses to two questions in the GWP: “Is corruption widespread throughout
#' the government or not” and “Is corruption widespread within businesses or
#' not?”.}
#' }
#'

#' @source World Happiness Report 2021, https://worldhappiness.report/ed/2021/.
#'
#' @examples
#' happy
#' name_country <- happy$Country
#'
"happy"
