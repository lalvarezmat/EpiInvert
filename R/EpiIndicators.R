#' @title 
#' \code{EpiIndicators_params} function to select EpiIndicators parameters
#'
#' @param x a list with parameters of EpiIndicators() function
#'
#' @return the EpiIndicators() parameters
#' 
#' @seealso \code{\link{EpiIndicators}} 
#' 
#' @export
EpiIndicators_params <- function(x=""){
  ## SET DEFAULTS
  defaults <- list(s_min = -10,
                   s_max = 25,
                   wr = 1000,
                   ws = 10,
                   s_init = -1e6,
                   s_end=-1e6,
                   r_init = -1e6,
                   r_end=-1e6
  )
  class(defaults) <- "EpiIndicators_config"
  ## MODIFY CONFIG WITH ARGUMENTS ##
  if(is.list(x)== TRUE){
    config <- x
    extra <- setdiff(names(config), names(defaults))
    if ((length(extra) > 0L)) {
      stop("Additional invalid options: ", paste(extra, collapse=", "))
    }
    config <- utils::modifyList(defaults, config)
    return(config)
  }
  else{
    return(defaults)
  }
}


#' @title 
#' \code{EpiIndicators} 
#' 
#' @description EpiIndicators estimates the ratio, r(t), and shift(delay), s(t), between 2 
#' epidemiological indicators f(t) and g(t) following the relation 
#'         r(t)*f(t) = g(t+s(t))
#'
#' @param df a dataframe with 3 columns: the first column corresponds to the date 
#' of each indicator value, the second column is the value of the first indicator 
#' f(t) and the third column is the value of the second indicator g(t). A zero 
#' value is expected in the case that the real value of an indicator is not 
#' available. Indicators must be smooth functions. So, for instance, the raw 
#' registered number of cases or deaths are not adequate to run the function. 
#' These particular indicators should be smoothed before executing EpiIndicators(), 
#' for instance you can use the restored indicator values obtained  by EpiInvert()
#'
#' 
#' @param config a list of the following optional parameters obtained using the 
#' function EpiIndicators_params(): 
#' s_min = -10,
#' \itemize{
#'  \item{s_min}{: min value allowed for the shift s(t) (default value -10)}
#'  \item{s_max}{: max value allowed for the shift s(t) (default value 25)}
#'  \item{wr}{: energy regularization parameter for the ratio  r(t) (default value 1000)}
#'  \item{ws}{: energy regularization parameter for the shift  s(t) (default value 10)}
#'  \item{s_init}{: manually fixed initial value (at time t=0) for s(t) (default value -1e6)
#'  by default s_init is not fixed and it is automatically estimated }
#'  \item{s_end}{: manually fixed final value (at the last time) for s(t) 
#'  (default value -1e6) by default s_end is not fixed and it is automatically 
#'  estimated}
#'  \item{r_init}{: manually fixed initial value (at time t=0) for r(t) (default value -1e6)
#'  by default r_init is not fixed and it is automatically estimated }
#'  \item{r_end}{: manually fixed final value (at the last time) for r(t) 
#'  (default value -1e6) by default r_end is not fixed and it is automatically 
#'  estimated}
#' }
#' 
#'
#' @return {
#'   A dataframe with the following columns :
#'   \itemize{
#'
#'   \item{date}{: the date of the indicator values. }
#'
#'   \item{f}{: the first indicator f(t).}
#'   
#'   \item{g}{: the second indicator g(t). }
#'   
#'   \item{r}{: the estimated ratio r(t) }
#'   
#'   \item{s}{: the estimated shift (delay)  s(t) }
#'   
#'   \item{f.r}{: the result of r(t)*f(t) }
#'   
#'   \item{g.s}{: the result of g(t+s(t))}
#' 
#'   }
#' }
#'
#' @details
#' EpiIndicators estimates the ratio, r(t), and shift(delay), s(t) between 2 
#' epidemiological indicators f(t) and g(t) following the relation 
#'           r(t)*f(t) = g(t+s(t))
#' a variational method is proposed to add regularity constraints to the 
#' estimates of r(t) and s(t).           
#' 
#' 
#' 
#' @author Luis Alvarez \email{lalvarez@ulpgc.es}
#'
#' @examples
#' ## load data of epidemiological indicators obtained from the World in data
#' ## organization
#' data("owid")
#'
#' ## Filter the data to get France epidemiological indicators
#' library(dplyr)
#' sel <- filter(owid,iso_code=="FRA")
#' 
#' ## Generate a dataframe with the dates and the cases and deaths restored 
#' ## using EpiInvert()
#' df<-data.frame(sel$date,sel$new_cases_restored_EpiInvert,sel$new_deaths_restored_EpiInvert)
#' 
#' ## Run EpiIndicators
#' res <- EpiIndicators(df)
#' 
#' ## Plot the results 
#' EpiIndicators_plot(res)
#' 
#' @useDynLib EpiInvert, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @export
EpiIndicators <- function(df,
                      config = EpiIndicators_params()){
  
  # CHECK IF df IS A DATAFRAME WITH 3 COLUMNS
  if(ncol(df)!=3){
    stop("The first argument is not a dataframe with 3 columns")
  }
  
  # CHECK WHETHER df[1,1] IS A DATE
  if(anyNA(as.Date(df[1,1], format= "%Y-%m-%d"))){
    stop("the first column of the dataframe must be in the format 'YYYY-MM-DD'")
  }
  
  # CHECK WHETHER df[1,2] IS A NUMBER
  if(is.numeric(df[1,2])!=TRUE){
    stop("the second column of the dataframe must be a number")
  }
  # CHECK WHETHER df[1,3] IS A NUMBER
  if(is.numeric(df[1,3])!=TRUE){
    stop("the third column of the dataframe must be a number")
  }
  
  # WE CALL THE MAIN Rcpp - C++ FUNCTION TO COMPUTE EpiIndicators
  results <- EpiIndicatorsC(df[,1],
                        df[,2],
                        df[,3],
                        config$s_min,
                        config$s_max,
                        config$wr,
                        config$ws,
                        config$r_init,
                        config$r_end,
                        config$s_init,
                        config$s_end
  )
  results <- as.data.frame(results)
  #class(results) <- "estimate_EpiIndicators"
  return(results)
  
}

#' @title 
#' \code{joint_indicators_by_date} 
#' 
#' @description generates a dataframe joining the dates and values
#' of 2 indicators
#'
#' @param date0 the dates of the first indicator.
#' 
#' @param i0 the values of the first indicator.
#' 
#' @param date1 the dates of the second indicator.
#' 
#' @param i1 the values of the second indicator.
#' 
#' @return { 
#'   A dataframe with the following columns :
#'   \itemize{
#'
#'   \item{date}{: all dates presented in any of the indicators. }
#'
#'   \item{f}{: the values of the first indicator. We assign 0
#'   in the case the data is not available for a given day.}
#'   
#'   \item{g}{: the values of the second indicator. We assign 0
#'   in the case the data is not available for a given day }
#' 
#'   }
#' }
#' @useDynLib EpiInvert, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @export
joint_indicators_by_date <- function(date0,i0,date1,i1)
{
  
  # CHECK DATES
  if(anyNA(as.Date(date0[1], format= "%Y-%m-%d"))){
    stop("the first parameter must be a date in the format 'YYYY-MM-DD'")
  }
  
  if(anyNA(as.Date(date1[1], format= "%Y-%m-%d"))){
    stop("the third parameter must be a date in the format 'YYYY-MM-DD'")
  }
  
  # CHECK NUMBERs
  if(is.numeric(i0[1])!=TRUE){
    stop("the second parameter must be a number")
  }
  if(is.numeric(i1[1])!=TRUE){
    stop("the fourth parameter must be a number")
  }
  
  # WE CALL THE MAIN Rcpp - C++ FUNCTION TO COMPUTE EpiIndicators
  results <- joint_indicators_by_dateC(as.character(date0),i0,as.character(date1),i1)
  results <- as.data.frame(results)
  
  return(results)
  
}

#' @title 
#' \code{apply_delay} 
#' 
#' @description apply a delay vector, s(t), to an indicator, g(t)
#'
#' @param g indicator.
#' 
#' @param s delay to apply.
#' 
#' 
#' @return { 
#'   A numeric vector with the result of apply the cevtor delay s to g. 
#' }
#' @useDynLib EpiInvert, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @export
apply_delay <- function(g,s)
{
  
  g<-as.numeric(g)
  s<-as.numeric(s)
  
  # WE CALL the C++ function 
  return(apply_shiftC(g,s))
  
}

