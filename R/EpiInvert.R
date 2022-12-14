#' @title 
#' \code{EpiInvert} estimates the reproduction number Rt and a restored 
#' incidence curve from the original daily incidence curve and the serial 
#' interval distribution. EpiInvert also corrects the festive and weekly biases 
#' present in the registered  daily incidence.  
#'
#' @param incid The original daily incidence curve (a numeric vector).
#'
#' @param last_incidence_date The date of the last value of the incidence curve
#' in the format YYYY-MM-DD. EpiInvert does not 
#' allow missing values. On days when a country does not report data, a zero must
#' be registered as the value associated with the incidence of that day.
#' 
#' @param festive_days The festive or anomalous dates in the format YYYY-MM-DD 
#' (a character vector). In these dates we "a priori" expect that the incidence
#' value is perturbed.
#'
#' @param config An object of class \code{estimate_R_config}, as returned by 
#' function \code{select_params}. The element of config are : 
#' 
#' \itemize{
#'  \item{si_distr}{: a numeric vector with the distribution of the serial 
#'  interval (the default value is an empty vector). If this vector is empty, 
#'  the serial interval is estimated using a parametric shifted log-normal. }
#'  \item{shift_si_distr}{: shift of the above user provided serial interval. 
#'  This shift can be negative, which means that secondary cases can show 
#'  symptoms before the primary cases (the default value is 0). }
#'  \item{max_time_interval}{: Maximum number of days used by 
#'  EpiInvert (the default value is 150, which means that EpiInvert uses the 
#'  last 150 days. The computational cost strongly depends on this value). }
#'  \item{mean_si}{: mean of the parametric shifted log-normal to approximate
#'   the serial interval (in the case the above si_distr vector is empty),
#'  (the default value is 12.267893). }
#'  \item{sd_si}{: standard deviation of the parametric shifted log-normal to 
#'  approximate the serial interval (in the case the above si_distr vector is empty),
#'   (the default value is  5.667547). }
#'  \item{shift_si=}{: shift of the parametric shifted log-normal to approximate
#'   the serial interval (in the case the above si_distr vector is empty),
#'  (the default value is -5). }
#'  \item{Rt_regularization_weight}{: regularization weight for Rt in the 
#'  variational model used by EpiInvert (the default value is 5). }
#'  \item{seasonality_regularization_weight}{: regularization weight for the 
#'  weekly bias correction factors in the variational model used by EpiInvert (the default value is 5).}
#'   \item{incidence_weekly_aggregated}{: a boolean value which determines if we use weekly aggregated 
#'  incidence. In such case  every week a single data is stored with the accumulated incidence in the last 
#'  7 days (the default value is FALSE).}
#'   \item{NweeksToKeepIncidenceSum}{: number of weeks to keep the value of the incidence accumulation. 
#'    The default is 2.}
#' }
#' 
#'
#' @return {
#'   an object of class \code{estimate_EpiInvert}, given by a list with 
#'   components:
#'   \itemize{
#'
#'   \item{i_original}{: the original daily incidence curve. In the case of weekly aggregated incidence, 
#'   we initialize the original curve assigning each day of the week 1/7 of the weekly aggregated value.}
#'
#'   \item{i_festive}{: the incidence after correction of the festive days bias.}
#'
#'   \item{i_bias_free}{: the incidence after correction of the festive and 
#'   weekly biases.}
#'
#'   \item{i_restored}{: the restored incidence obtained using the renewal 
#'   equation.}
#'
#'   \item{Rt}{: the reproduction number Rt obtained by inverting the renewal 
#'   equation.}
#' 
#'   \item{Rt_CI95}{: 95\% confidence interval radius for the value of Rt taking 
#'   into account the variation of Rt when more days are added to the estimation. 
#'   }
#' 
#'   \item{seasonality}{: the weekly bias correction factors.}
#' 
#'   \item{dates}{: a vector of dates corresponding to the incidence curve. }
#' 
#'   \item{festive}{: boolean associated to each incidence value to check if it 
#'   has been considered as a festive or anomalous day. }
#' 
#'   \item{epsilon}{: normalized error curve obtained as 
#'    (i_bias_free-i_restored)/i_restored^a.}
#' 
#'   \item{power_a}{: the power, a, which appears in the above expression of the 
#'   normalized error. }
#' 
#'   \item{si_distr}{: values of the distribution of the serial interval used in
#'    the EpiInvert estimation.}
#' 
#'   \item{shift_si_distr}{: shift of the above distribution of the serial 
#'   interval si_distr.}
#' 
#'   }
#' }
#'
#' @details
#' EpiInvert estimates the reproduction number Rt and a 
#' restored incidence curve by inverting the renewal equation : 
#' 
#' \deqn{{i(t) = \sum_k i(t-k)R(t-k)\Phi(k)}}
#' 
#' using a variational 
#' formulation. The theoretical foundations of the method can be found 
#' in references [1] and [2].
#' 
#' 
#' @author Luis Alvarez \email{lalvarez@ulpgc.es}
#' @references {
#' 
#' [1] Alvarez, L.; Colom, M.; Morel, J.D.; Morel, J.M. Computing the daily 
#' reproduction number of COVID-19 by inverting the renewal
#' equation using a variational technique. Proc. Natl. Acad. Sci. USA, 2021.
#' 
#' [2] Alvarez, Luis, Jean-David Morel, and Jean-Michel Morel. "Modeling 
#' COVID-19 Incidence by the Renewal Equation after Removal of Administrative 
#' Bias and Noise" Biology 11, no. 4: 540. 2022. 
#' 
#' [3] Ritchie, H. et al. Coronavirus Pandemic (COVID-19), OurWorldInData.org. 
#' Available online: 
#' https://ourworldindata.org/coronavirus-source-data 
#' (accessed on 5 May 2022).
#' 
#' }
#' @examples
#' ## load data on COVID-19 daily incidence up to 2022-05-05 for France, 
#' ## and Germany (taken from the official government data) and for UK and 
#' ## the USA taken from reference [3]
#' data(incidence)
#'
#' ## EpiInvert execution for USA with no festive days specification
#' ## using the incidence 70 days in the past
#' res <- EpiInvert(incidence$USA,
#' "2022-05-05",
#' "1000-01-01",
#' select_params(list(max_time_interval = 70))
#' )
#' 
#' ## Plot of the results
#' EpiInvert_plot(res)
#'
#' ## load data of festive days for France, Germany, UK and the USA
#' data(festives)
#' 
#' ## EpiInvert execution for France with festive days specification using
#' ## 70 days in the past 
#' res <- EpiInvert(incidence$FRA,"2022-05-05",festives$FRA,
#'                  select_params(list(max_time_interval = 70)))
#' 
#' ## Plot of the incidence between "2022-04-01" and "2022-05-01"
#' EpiInvert_plot(res,"incid","2022-04-01","2022-05-01")
#' 
#' ## load data of a serial interval
#' data("si_distr_data")
#' 
#' ## EpiInvert execution for Germany using the uploaded serial interval shifted
#' ## -2 days, using the incidence 70 days in the past
#' res <- EpiInvert(incidence$DEU,"2022-05-05",festives$DEU,
#'        EpiInvert::select_params(list(si_distr = si_distr_data,
#'        shift_si_distr=-2,max_time_interval = 70)))
#'        
#' ## Plot of the serial interval used (including the shift)
#' EpiInvert_plot(res,"SI")     
#'        
#' ## EpiInvert execution for UK changing the default values of the parametric
#' ## serial interval (using a shifted log-normal) using 70 days in the past
#' res <- EpiInvert(incidence$UK,"2022-05-05",festives$UK,
#'        EpiInvert::select_params(list(mean_si = 11,sd_si=6,shift_si=-1,max_time_interval = 70))
#'        )
#'        
#' ## Plot of the reproduction number Rt including an empiric 95\% confidence 
#' ## interval of the variation of Rt. To calculate Rt on each day t, EpiInvert 
#' ## uses the past days (t'<=t) and the future days (t'>t) when available. 
#' ## Therefore, the EpiInvert estimate of Rt varies when there are more days 
#' ## available. This confidence interval reflects the expected variation of Rt 
#' ## as a function of the number of days after t available.
#' EpiInvert_plot(res,"R")
#' 
#' @useDynLib EpiInvert, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @export
EpiInvert <- function(incid,
                      last_incidence_date,
                      festive_days = rep("1000-01-01",2),
                      config = EpiInvert::select_params()) {

  # CHECK IF incid IS A NUMERIC VECTOR
  if(is.numeric(incid)!=TRUE){
    stop("The first argument of EpiInvert() is not numeric")
  }
  
  last_incidence_date <- format(as.Date(last_incidence_date), "%Y-%m-%d")
  festive_days <- format(as.Date(festive_days), "%Y-%m-%d")
  
  # CHECK IF last_incidence_date IS A DATE IN THE FORMAT YYYY-MM-DD
  if(is.character(last_incidence_date)!=TRUE){
    stop("last_incidence_date must be in the format 'YYYY-MM-DD'")
  }
  #if( dat_time_check_fn(last_incidence_date)== FALSE){
  if(anyNA(as.Date(last_incidence_date, format= "%Y-%m-%d"))){
    stop("last_incidence_date must be in the format 'YYYY-MM-DD'")
  }

  # CHECK IF festive_days IS OF CHARACTER TYPE
  if(is.character(festive_days)!=TRUE){
    stop("festive_days must be in the format 'YYYY-MM-DD")
  }

  # WE KEEP ONLY THE FESTIVE DAYS IN THE FORMAT YYYY-MM-DD
  festive_days2  <- vector(mode='character',length=0)
  for (val in festive_days) {
    if(!anyNA(as.Date(val, format= "%Y-%m-%d"))){
    #if(dat_time_check_fn(val)==TRUE) 
      festive_days2 <- append(festive_days2,val)
    }
  }

  # WE STOP IF ANY OF THE FESTIVE DAYS HAS THE FORMAT YYYY-MM-DD
  if(length(festive_days2)==0){
    stop("festive_days must be in the format 'YYYY-MM-DD")
  }
  
  # COMPUTING THE LAST DAY
  # i <- length(incid)
  #  while (i >=0) {
  #    if(incid[i]!=0){
  #     break
  #   }
  #   last_incidence_date=as.Date(last_incidence_date)-1
  #  i <- i - 1
  #  }
  #  last_incidence_date=toString(last_incidence_date)
  
  # WE CALL THE MAIN Rcpp - C++ FUNCTION TO COMPUTE EpiInvert 
  results <- EpiInvertC(incid,
                        last_incidence_date,
                        festive_days2,
                        config$si_distr,
                        config$shift_si_distr,
                        config$max_time_interval,
                        config$mean_si,
                        config$sd_si,
                        config$shift_si,
                        config$Rt_regularization_weight,
                        config$seasonality_regularization_weight,
                        config$incidence_weekly_aggregated,
                        config$NweeksToKeepIncidenceSum
                      )
  
  class(results) <- "estimate_EpiInvert"
  return(results)


}



#' @title 
#' \code{EpiInvertForecast} computes a 28-day forecast of the restored incidence 
#' curve including a 95% confidence interval 
#' using the weekly seasonality, from the forecasted restored incidence curve we
#' also estimate a  28-day forecast of the original incidence curve. 
#'
#' @param EpiInvert_result output list of the EpiInvert execution including, in particular,  
#' the restored incidence curve and the seasonality. 
#'
#' @param restored_incidence_database a database including 27,418  samples of different
#' restored incidence curves computed by EpiInvert using real data. Each 
#' restored incidence curve includes the last 56 values of the sequence. That is
#' this database can be viewed as a matrix of size 27,418  X 56
#' 
#' @param type string with the forecast option. It can be "mean" or "median". 
#' 
#' @param NumberForecastAdditionalDays The number of forecast days is 28. With this
#' parameter you can add extra forecast days using linear extrapolation. 
#' 
#' @param trend_sentiment "a priori" knowledge about the future incidence evolution. 
#'    == 0 means that you are neutral about the future trend 
#'    > 0  means that you expect that the future trend is higher than the expected one 
#'         using all database curves. the value represents the percentage of database 
#'         curves removed before computing the forecast The curves removed are the ones 
#'         with lowest growth in the last 28 days.  
#'    < 0  means that you expect that the future trend is higher than the expected one 
#'         using all database curves. the meaning of the value is similar to the previous 
#'         case, but removing the curves with the highest growth in the last 28 days.
#'    
#'
#' 
#' @return {
#'   a list with components:
#'   \itemize{
#'   
#'      \item{dates}{: a vector of dates corresponding to the forecast days. }
#'   
#'      \item{i_restored_forecast}{: a numeric vector with the forecast of the 
#'      restored incidence curve for the next 28 days.}
#'   
#'      \item{i_original_forecast}{: a numeric vector with the forecast of the 
#'      original incidence curve for the next 28 days.}
#'
#'      \item{i_restored_forecast_CI50}{: radius of an empiric confidence interval,  
#'        with percentile 50, for the restored incidence forecast following the 
#'        number of days passed since the current day. 
#'      }
#'   
#'      \item{i_restored_forecast_CI75}{: radius of an empiric confidence interval,  
#'         with percentile 75, for the restored incidence forecast following the 
#'         number of days passed since the current day. 
#'      }
#'
#'     \item{i_restored_forecast_CI90}{: radius of an empiric confidence interval,  
#'        with percentile 90, for the restored incidence forecast following the 
#'        number of days passed since the current day. 
#'     }
#'
#'     \item{i_restored_forecast_CI95}{: radius of an empiric confidence interval,  
#'        with percentile 95, for the restored incidence forecast following the 
#'        number of days passed since the current day. 
#'     }
#' 
#'   }
#' }
#'
#' @details
#' EpiInvertForecast estimates a forecast of the restored incidence curve 
#' using a weighted average of 27,418  restored incidence curves previously 
#' estimated by EpiInvert and stored in the database "restored_incidence_database". 
#' The weight, in the average computation, of each restored incidence curve 
#' of the database depends on the similarity between the current curve in the last
#' 28 days and the first 28 days of the database curve. Each database curve 
#' contains 56 days. The first 28 days are used for comparison with the current 
#' curve and the last 28 days are used for forecasting. 
#' 
#' 
#' 
#' @author Luis Alvarez \email{lalvarez@ulpgc.es}
#' @references {
#' 
#' [1] Alvarez, L.; Colom, M.; Morel, J.D.; Morel, J.M. Computing the daily 
#' reproduction number of COVID-19 by inverting the renewal
#' equation using a variational technique. Proc. Natl. Acad. Sci. USA, 2021.
#' 
#' [2] Alvarez, Luis, Jean-David Morel, and Jean-Michel Morel. "Modeling 
#' COVID-19 Incidence by the Renewal Equation after Removal of Administrative 
#' Bias and Noise" Biology 11, no. 4: 540. 2022. 
#' 
#' [3] Ritchie, H. et al. Coronavirus Pandemic (COVID-19), OurWorldInData.org. 
#' Available online: 
#' https://ourworldindata.org/coronavirus-source-data 
#' (accessed on 5 May 2022).
#' 
#' [4] Alvarez, Luis, Jean-David Morel, and Jean-Michel Morel. EpiInvertForecast
#' Available online: 
#' https://ctim.ulpgc.es/covid19/EpiInvertForecastPaper.html 
#' 
#' }
#' @examples
#' ## load data on COVID-19 daily incidence up to 2022-05-05 for France, 
#' ## and Germany (taken from the official government data) and for UK and 
#' ## the USA taken from reference [3]
#' data(incidence)
#' 
#' ## load of the database of restored incidence curves. 
#' data("restored_incidence_database")
#'
#' ## EpiInvert execution for USA with no festive days specification
#' ## using the incidence 90 days in the past
#' res <- EpiInvert(incidence$USA,
#' "2022-05-05",
#' "1000-01-01",
#' select_params(list(max_time_interval = 90))
#' )
#' 
#' ## EpiInvertForecast execution using the EpiInvert results obtained by USA
#' forecast <-  EpiInvertForecast(res,restored_incidence_database)
#'
#'
#' 
#' @useDynLib EpiInvert, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @export
EpiInvertForecast <- function(EpiInvert_result,restored_incidence_database,type="median",trend_sentiment=0,NumberForecastAdditionalDays=0) {
  
  # CHECK IF restored_incidence_database IS NUMERIC 
  if(is.numeric(restored_incidence_database)!=TRUE){
    stop("The second argument of EpiInvertForecast() is not numeric")
  }
  
  # LAST DAY OF AVAILABLE DATA
  str <- EpiInvert_result$dates[length(EpiInvert_result$dates)]
  
  # WE CALL THE MAIN Rcpp - C++ FUNCTION TO COMPUTE EpiInvertForeCast 
  results <- EpiInvertForecastC(
    EpiInvert_result$i_original,
    EpiInvert_result$i_restored,
    str,
    EpiInvert_result$seasonality,
    restored_incidence_database,
    type,
    NumberForecastAdditionalDays,
    trend_sentiment
  )
  
  return(results)
}

