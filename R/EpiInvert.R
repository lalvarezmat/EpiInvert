dat_time_check_fn <- function(dat_time) {
  if (!anyNA(as.Date(dat_time, format= "%Y-%m-%d"))) return(TRUE)
  else return(FALSE)

}

EpiInvert <- function(incid,
                      last_incidence_date,
                      festive_days = rep("1000-01-01",2),
                      config = make_config()) {

  if(is.numeric(incid)!=TRUE){
    stop("The first argument of EpiInvert() is not numeric")
  }
  if(is.character(last_incidence_date)!=TRUE){
    stop("last_incidence_date must be in the format 'YYYY-MM-DD'")
  }
  if( dat_time_check_fn(last_incidence_date)== FALSE){
    stop("last_incidence_date must be in the format 'YYYY-MM-DD'")
  }

  if(is.character(festive_days)!=TRUE){
    stop("festive_days must be in the format 'YYYY-MM-DD")
  }

  festive_days2  <- vector(mode='character',length=0)

  for (val in festive_days) {
    if(dat_time_check_fn(val)==TRUE) festive_days2 <- append(festive_days2,val)
  }

  if(length(festive_days2)==0){
    stop("festive_days must be in the format 'YYYY-MM-DD")
  }


  return(EpiInvertC(incid,
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
  )


}
