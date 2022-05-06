test_that("EpiInvert using USA data", {
  data("incidence")
  data("festives")
  data("si_distr_data")
  
  # test using parametric serial interval 
  res <- EpiInvert(incidence$USA,incidence$date[length(incidence$date)],festives$USA)
  x <- round(res$power_a,digits=3)
  expect_equal(x, 1.097)
  
  # test using non parametric serial interval 
  res <- EpiInvert(incidence$USA,incidence$date[length(incidence$date)],festives$USA,EpiInvert::make_config(list(si_distr = si_distr_data,shift_si_distr=1)))
  x <- round(res$power_a,digits=3)
  expect_equal(x, 1.086)
  
  # test without using festive days
  res <- EpiInvert(incidence$USA,incidence$date[length(incidence$date)])
  x <- round(res$power_a,digits=3)
  expect_equal(x, 1.063)
  
  # test using a different time interval 
  res <-  EpiInvert(incidence$USA,incidence$date[length(incidence$date)],festives$USA, EpiInvert::make_config(list(max_time_interval = 100)))
  x <- round(res$power_a,digits=3)
  expect_equal(x, 0.595)
  
  # test using a different parametric serial interval 
  res <-  EpiInvert(incidence$USA,incidence$date[length(incidence$date)],festives$USA, EpiInvert::make_config(list( mean_si = 11,sd_si=7,shift_si=-2)))
  x <- round(res$power_a,digits=3)
  expect_equal(x,1.033)
  
  # test using different regularization values
  res <-  EpiInvert(incidence$USA,incidence$date[length(incidence$date)],festives$USA, EpiInvert::make_config(list( Rt_regularization_weight=10,seasonality_regularization_weight=20)))
  x <- round(res$power_a,digits=3)
  expect_equal(x,1.053)
  
  

  
})
