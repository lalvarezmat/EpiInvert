#' A dataset containing daily incidence of COVID-19 for France, 
#' Germany, UK and the USA
#'
#'
#' @format A data frame with 5 variables:
#' \describe{
#'   \item{date}{date of the incidence.}
#'   \item{FRA}{incidence of France.}
#'   \item{DEU}{incidence of Germany.}
#'   \item{UK}{incidence of UK.}
#'   \item{USA}{incidence of USA.}
#' }
#' 
#' An updated version of this dataset can be found at 
#' https://www.ctim.es/covid19/incidence.csv
#' 
#' @source \url{https://github.com/owid/covid-19-data/tree/master/public/data}
#' @source \url{https://www.santepubliquefrance.fr/dossiers/coronavirus-covid-19/coronavirus-chiffres-cles-et-evolution-de-la-covid-19-en-france-et-dans-le-monde}
#' @source \url{https://experience.arcgis.com/experience/478220a4c454480e823b17327b2bf1d4/page/page_1/}
"incidence"


#' A dataset containing weekly aggregated incidence of COVID-19 for France, 
#' Germany, UK and the USA
#'
#'
#' @format A data frame with 5 variables:
#' \describe{
#'   \item{date}{date of the weekly aggregated incidence.}
#'   \item{FRA}{weekly aggregated incidence of France.}
#'   \item{DEU}{weekly aggregated incidence of Germany.}
#'   \item{UK}{weekly aggregated incidence of UK.}
#'   \item{USA}{weekly aggregated incidence of USA.}
#' }
#' 
#' An updated version of this dataset can be found at 
#' https://www.ctim.es/covid19/incidence.csv
#' 
#' @source \url{https://github.com/owid/covid-19-data/tree/master/public/data}
#' @source \url{https://www.santepubliquefrance.fr/dossiers/coronavirus-covid-19/coronavirus-chiffres-cles-et-evolution-de-la-covid-19-en-france-et-dans-le-monde}
#' @source \url{https://experience.arcgis.com/experience/478220a4c454480e823b17327b2bf1d4/page/page_1/}
"incidence_weekly_aggregated"


#' A dataset containing festive days in  France, 
#' Germany, UK and the USA
#'
#'
#' @format A list with 4 variables:
#' \describe{
#'   \item{FRA}{festive day of France}
#'   \item{DEU}{festive day  of Germany}
#'   \item{UK}{festive day  of UK}
#'   \item{USA}{festive day  of USA}
#'  
#' }
"festives"

#' A dataset containing the values of a serial interval
#'
#' @format A numeric vector with 1 variable representing 
#' the serial interval 
#' 
"si_distr_data"

#' A dataset of restored daily incidence trend curves
#' @description A dataset including 27,418 samples of different
#' restored incidence curves computed by EpiInvert using real data. Each 
#' restored incidence curve includes the last 56 values of the sequence. 
#'
#' @format A 27,418 X 56 numeric matrix
#' 
"restored_incidence_database"

#' A datased with COVID-19 indicators
#' @description A dataset containing COVID-19 epidemiological indicators for Canada, France, 
#' Germany, Italy, UK and the USA from Our World in data organization
#' https://github.com/owid/covid-19-data/tree/master/public/data up to 2022-11-28. 
#' In the case a data value is not avalaible for a given day we assign the value 0 
#' to the indicator
#'
#'
#' @format A dataframe with 13 variables:
#' \describe{
#'   \item{iso_code}{iso code of the country}
#'   \item{location}{country name}
#'   \item{date}{date of the indicator value }
#'   \item{new_cases}{new confirmed cases}
#'   \item{new_cases_smoothed}{new confirmed cases smoothed}
#'   \item{new_cases_restored_EpiInvert}{new confirmed cases restored using EpiInvert}
#'   \item{new_deaths}{new deaths attributed to COVID-19}
#'   \item{new_deaths_smoothed}{new deaths smoothed}
#'   \item{new_deaths_restored_EpiInvert}{new deaths restored using EpiInvert}
#'   \item{icu_patients}{number of COVID-19 patients in intensive care units (ICUs) 
#'   on a given day}
#'   \item{hosp_patients}{number of COVID-19 patients in hospital on a given day}
#'   \item{weekly_icu_admissions}{number of COVID-19 patients newly admitted to 
#'   intensive care units (ICUs) in a given week (reporting date and the preceeding 
#'   6 days)}
#'   \item{weekly_hosp_admissions}{number of COVID-19 patients newly admitted to 
#'hospitals in a given week (reporting date and the preceeding 6 days)}
#'  
#' }
"owid_data"


