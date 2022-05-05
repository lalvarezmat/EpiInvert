make_config <- function(L=""){



  #config <- list(...)


  ## SET DEFAULTS
  defaults <- list(si_distr = vector(mode="numeric", length=0),
                   shift_si_distr = 0,
                   max_time_interval = 150,
                   mean_si = 12.267893,
                   sd_si = 5.667547,
                   shift_si=-5.,
                   Rt_regularization_weight=5.,
                   seasonality_regularization_weight=5.)

  ## MODIFY CONFIG WITH ARGUMENTS ##
  if(is.list(L)== TRUE){
    config <- L
    config <- modify_defaults(defaults, config)
    return(config)
  }
  else{
    return(defaults)
  }

}

## This function was contributed by Rich Fitzjohn. It modifies default arguments
## using user-provided values. The argument 'strict' triggers and error
## behaviour: if strict==TRUE: all new values need to be part of the defaults.

modify_defaults <- function(defaults, x, strict = TRUE) {
  extra <- setdiff(names(x), names(defaults))
  if (strict && (length(extra) > 0L)) {
    stop("Additional invalid options: ", paste(extra, collapse=", "))
  }
  utils::modifyList(defaults, x, keep.null = TRUE) # keep.null is needed here

}
