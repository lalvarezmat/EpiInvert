test_that("EpiInvert using USA data", {
  data("incidence")
  data("festives")
  data("si_distr_data")
  
  # test using parametric serial interval and the last 70 incidence values
  res <- EpiInvert(incidence$USA,
                   incidence$date[length(incidence$date)],
                   festives$USA,
                   EpiInvert::select_params(list(max_time_interval = 70))
                   )
  x <- round(res$power_a+0.0018224,digits=2)
  expect_equal(x,-0.31)
  
  # test using non parametric serial interval and the last 70 incidence values
  res <- EpiInvert(incidence$USA,
                   incidence$date[length(incidence$date)],
                   festives$USA,
                   EpiInvert::select_params(list(si_distr = si_distr_data,shift_si_distr=1,max_time_interval = 70))
                   )
  x <- round(res$power_a+0.00809614,digits=2)
  expect_equal(x, -0.22)
  
  # test without using festive days and the last 70 incidence values
  res <- EpiInvert(incidence$USA,
                   incidence$date[length(incidence$date)],
                   "1000-01-01",
                   EpiInvert::select_params(list(max_time_interval = 70))
                   )
  x <- round(res$power_a+0.0018224,digits=2)
  expect_equal(x, -0.31)
  
  # test using a different time interval 
  res <-  EpiInvert(incidence$USA,
                    incidence$date[length(incidence$date)],
                    festives$USA, 
                    
                    EpiInvert::select_params(list(max_time_interval = 80)))
  x <- round(res$power_a+0.0060459,digits=2)
  expect_equal(x, -0.43)
  
  # test using a different parametric serial interval and the last 70 incidence values
  res <-  EpiInvert(incidence$USA,
                    incidence$date[length(incidence$date)],
                    festives$USA, 
                    EpiInvert::select_params(list( mean_si = 11,sd_si=7,shift_si=-2,max_time_interval = 70))
                    )
  x <- round(res$power_a+0.0022343,digits=2)
  expect_equal(x,-0.19)
  
  # test using different regularization values and the last 70 incidence values
  res <-  EpiInvert(incidence$USA,
                    incidence$date[length(incidence$date)],
                    festives$USA, 
                    EpiInvert::select_params(list( Rt_regularization_weight=10,seasonality_regularization_weight=20,max_time_interval = 70))
                    )
  x <- round(res$power_a+0.0022683,digits=2)
  expect_equal(x,-0.32)
  
  
})
