#' @title 
#' \code{select_params} function to select EpiInvert parameters
#'
#' @param x a list with elements of the class  \code{estimate_R_config}
#'
#' @return an object of class \code{estimate_R_config}.
#' 
#' @seealso \code{\link{EpiInvert}} 
#' 
#' @export
select_params <- function(x=""){

  ## SET DEFAULTS
  defaults <- list(si_distr = vector(mode="numeric", length=0),
                   shift_si_distr = 0,
                   max_time_interval = 150,
                   mean_si = 12.267893,
                   sd_si = 5.667547,
                   shift_si=-5.,
                   Rt_regularization_weight=5.,
                   seasonality_regularization_weight=5.,
                   incidence_weekly_aggregated=FALSE,
                   NweeksToKeepIncidenceSum=2)
  
  class(defaults) <- "estimate_R_config"

  ## MODIFY CONFIG WITH ARGUMENTS ##
  if(is.list(x)== TRUE){
    config <- x
    #config <- EpiInvert::modify_defaults(defaults, config)
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


#modify_defaults <- function(defaults, x, strict = TRUE) {
#extra <- setdiff(names(x), names(defaults))
# if (strict && (length(extra) > 0L)) {
#   stop("Additional invalid options: ", paste(extra, collapse=", "))
# }
# utils::modifyList(defaults, x, keep.null = TRUE) # keep.null is needed here

#}
