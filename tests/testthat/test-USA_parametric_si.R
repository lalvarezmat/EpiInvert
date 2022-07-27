test_that("EpiInvert using USA data", {
  data("incidence")
  data("festives")
  data("si_distr_data")
  
  # test using parametric serial interval 
  res <- EpiInvert(incidence$USA,
                   incidence$date[length(incidence$date)],
                   festives$USA
                   )
  x <- round(res$power_a+0.003,digits=2)
  expect_equal(x, 1.1)
  
  # test using non parametric serial interval 
  res <- EpiInvert(incidence$USA,
                   incidence$date[length(incidence$date)],
                   festives$USA,
                   EpiInvert::select_params(list(si_distr = si_distr_data,shift_si_distr=1)))
  x <- round(res$power_a+0.004,digits=2)
  expect_equal(x, 1.09)
  
  # test without using festive days
  res <- EpiInvert(
    incidence$USA,
    incidence$date[length(incidence$date)]
    )
  x <- round(res$power_a+0.007,digits=2)
  expect_equal(x, 1.07)
  
  # test using a different time interval 
  res <-  EpiInvert(incidence$USA,
                    incidence$date[length(incidence$date)],
                    festives$USA, 
                    EpiInvert::select_params(list(max_time_interval = 100))
                    )
  x <- round(res$power_a-0.005,digits=2)
  expect_equal(x, 0.59)
  
  # test using a different parametric serial interval 
  res <-  EpiInvert(incidence$USA,
                    incidence$date[length(incidence$date)],
                    festives$USA, 
                    
                    EpiInvert::select_params(list( mean_si = 11,sd_si=7,shift_si=-2))
                    )
  x <- round(res$power_a-0.003,digits=2)
  expect_equal(x,1.03)
  
  # test using different regularization values
  res <-  EpiInvert(incidence$USA,
                    incidence$date[length(incidence$date)],
                    festives$USA, 
                    EpiInvert::select_params(list( Rt_regularization_weight=10,seasonality_regularization_weight=20)))
  x <- round(res$power_a-0.003,digits=2)
  expect_equal(x,1.05)
  
  # test using weekly aggregated incidence 
  data(incidence_weekly_aggregated)
  res <- EpiInvert(incidence_weekly_aggregated$FRA,"2022-05-05",festives$FRA,
                    select_params(list(incidence_weekly_aggregated = TRUE)))
  x <- round(res$power_a-0.001,digits=2)
  expect_equal(x,1.25)
  
  # EpiInvertForecast function test
  data("restored_incidence_database")
  forecast <-  EpiInvertForecast(res,restored_incidence_database)
  x <- round(forecast$i_restored_forecast[1]-0.19,digits=0)
  expect_equal(x,36024)
  
})
