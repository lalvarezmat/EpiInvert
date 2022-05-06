#' @title 
#' \code{EpiInvert} estimates the reproduction number Rt and a 
#' restored incidence curve from the original 
#' daily incidence curve and the serial interval distribution.
#'
#' @param incid The original incidence curve (a numeric vector).
#'
#' @param last_incidence_date The date of the last value of the incidende curve
#' in the format YYYY-MM-DD.
#' 
#' @param festive_days The festive or anomalous dates in the format YYYY-MM-DD 
#' (a character vector). In these dates we "a priori" expect that the incidence
#' value is perturbed.
#'
#' @param config An object of class \code{estimate_R_config}, as returned by 
#' function \code{make_config}. The element of config are : 
#' 
#' \itemize{
#'  \item{si_distr}{: a numeric vector with the distribution of the serial 
#'  interval. If this vector is empty, the serial interval is estimated using 
#'  a parametric shifted log-normal (the default value is an empty vector)}
#'  \item{shift_si_distr}{: shift of the above user provided serial interval. This shift can be negative, which 
#'  means that secondary cases can show symptoms before the primary cases (the default value is 0) }
#'  \item{max_time_interval}{: Maximum number of days used by 
#'  EpiInvert (the default value is 150, which means that EpiInvert uses the last 150 days. The computational cost strongly depends on this value) }
#'  \item{mean_si}{: mean of the parametric shifted log-normal to approximate
#'   the serial interval (the default value is 12.267893) }
#'  \item{sd_si}{: standard deviation of the parametric shifted log-normal to 
#'  approximate the serial interval (the default value is  5.667547) }
#'  \item{shift_si=}{: shift of the parametric shifted log-normal to approximate
#'   the serial interval (the default value is -5) }
#'  \item{Rt_regularization_weight}{: regularization weight for Rt in the 
#'  variational model used by EpiInvert (the default value is 5) }
#'  \item{seasonality_regularization_weight}{: regularization weight for the 
#'  weekly bias correction factors in the variational model used by EpiInvert (the default value is 5)}
#'  
#' }
#' 
#'
#' @return {
#'   an object of class \code{estimate_EpiInvert}, given by a list with 
#'   components:
#'   \itemize{
#'
#'   \item{i_original}{: the original daily incidence time series}
#'
#'   \item{i_festive}{: the incidence after correction of the festive days bias}
#'
#'   \item{i_bias_free}{: the incidence after correction of the festive and 
#'   weekly biases}
#'
#'   \item{i_restored}{: the restored incidence obtained using the renewal 
#'   equation}
#'
#'   \item{Rt}{: the reproduction number Rt obtained by inverting the renewal 
#'   equation}
#' 
#'   \item{Rt_CI95}{: 95% confidence interval radius for the value of Rt taking 
#'   into account the variation of Rt when more days are added to the estimation. 
#'   }
#' 
#'   \item{seasonality}{: the weekly bias correction factors}
#' 
#'   \item{dates}{: a vector of dates corresponding to the incidence curve }
#' 
#'   \item{festive}{: boolean associated to each incidence value to check if it 
#'   has been considered as a festive or anomalous day }
#' 
#'   \item{epsilon}{: normalized error curve obtained as 
#'    (i_bias_free-i_restored)/i_restored^a}
#' 
#'   \item{powr_a}{: the power, a, which appears in the above expression of the 
#'   normalized error }
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
#' restored incidence curve from the original 
#' daily incidence time series and the serial interval distribution.
#' 
#' @author Luis Alvarez \email{lalvarez@ulpgc.es}
#' @references {
#' Cori, A. et al. A new framework and software to estimate time-varying
#' reproduction numbers during epidemics (AJE 2013).
#' Wallinga, J. and P. Teunis. Different epidemic curves for severe acute
#' respiratory syndrome reveal similar impacts of control measures (AJE 2004).
#' Reich, N.G. et al. Estimating incubation period distributions with coarse
#' data (Statis. Med. 2009)
#' }
#' @examples
#' ## load data on pandemic flu in a school in 2009
#' data("Flu2009")
#'
#' ## estimate the reproduction number (method "non_parametric_si")
#' ## when not specifying t_start and t_end in config, they are set to estimate
#' ## the reproduction number on sliding weekly windows                          
#' res <- estimate_R(incid = Flu2009$incidence, 
#'                   method = "non_parametric_si",
#'                   config = make_config(list(si_distr = Flu2009$si_distr)))
#' plot(res)
#'
#' ## the second plot produced shows, at each each day,
#' ## the estimate of the reproduction number over the 7-day window 
#' ## finishing on that day.
#' 
#' ## to specify t_start and t_end in config, e.g. to have biweekly sliding
#' ## windows      
#' t_start <- seq(2, nrow(Flu2009$incidence)-13)   
#' t_end <- t_start + 13                 
#' res <- estimate_R(incid = Flu2009$incidence, 
#'                   method = "non_parametric_si",
#'                   config = make_config(list(
#'                       si_distr = Flu2009$si_distr, 
#'                       t_start = t_start, 
#'                       t_end = t_end)))
#' plot(res)
#'
#' ## the second plot produced shows, at each each day,
#' ## the estimate of the reproduction number over the 14-day window 
#' ## finishing on that day.
#'
#' ## example with an incidence object
#'
#' ## create fake data
#' library(incidence)
#' data <- c(0,1,1,2,1,3,4,5,5,5,5,4,4,26,6,7,9)
#' location <- sample(c("local","imported"), length(data), replace=TRUE)
#' location[1] <- "imported" # forcing the first case to be imported
#'
#' ## get incidence per group (location)
#' incid <- incidence(data, groups = location)
#'
#' ## Estimate R with assumptions on serial interval
#' res <- estimate_R(incid, method = "parametric_si",
#'                   config = make_config(list(
#'                   mean_si = 2.6, std_si = 1.5)))
#' plot(res)
#' ## the second plot produced shows, at each each day,
#' ## the estimate of the reproduction number over the 7-day window
#' ## finishing on that day.
#'
#' ## estimate the reproduction number (method "parametric_si")
#' res <- estimate_R(Flu2009$incidence, method = "parametric_si",
#'                   config = make_config(list(mean_si = 2.6, std_si = 1.5)))
#' plot(res)
#' ## the second plot produced shows, at each each day,
#' ## the estimate of the reproduction number over the 7-day window
#' ## finishing on that day.
#'
#' ## estimate the reproduction number (method "uncertain_si")
#' res <- estimate_R(Flu2009$incidence, method = "uncertain_si",
#'                   config = make_config(list(
#'                   mean_si = 2.6, std_mean_si = 1,
#'                   min_mean_si = 1, max_mean_si = 4.2,
#'                   std_si = 1.5, std_std_si = 0.5,
#'                   min_std_si = 0.5, max_std_si = 2.5,
#'                   n1 = 100, n2 = 100)))
#' plot(res)
#' ## the bottom left plot produced shows, at each each day,
#' ## the estimate of the reproduction number over the 7-day window
#' ## finishing on that day.
#'
#'


EpiInvert <- function(incid,
                      last_incidence_date,
                      festive_days = rep("1000-01-01",2),
                      config = EpiInvert::make_config()) {

  # CHECK IF incid IS A NUMERIC VECTOR
  if(is.numeric(incid)!=TRUE){
    stop("The first argument of EpiInvert() is not numeric")
  }
  
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
                        config$seasonality_regularization_weight
                      )
  
  class(results) <- "estimate_EpiInvert"
  return(results)


}
