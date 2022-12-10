# EpiInvert 0.3.1 (December 2022)

* I have received a message from the CRAN staff pointed out  an error running 
‘testthat.R’ in "r-release-macos-arm64" and "r-oldrel-macos-arm64"  due to  
differences in the quality of floating point arithmetic operations. I have 
amplified the  margin of error in the estimated number to avoid this problem. 

# EpiInvert 0.3.0 (December 2022)

* We include EpiIndicators, a procedure for the estimation of the delay and ratio between epidemiological indicators.

# EpiInvert 0.2.1 (October 2022)

* In the EpiInvertForecast we include the option of using the median, instead of
a weighted mean of the incidence curve dataset. 

* In the EpiInvertForecast we add a "trend sentiment" parameter to 
include "a priori" expectations about the incidence trend. The default value is 
"neutral" which means that you do not have any reason to believe that the prediction 
will be different to the initially predicted by the method. A positive "trend sentiment" 
means that you believe that the true value will be higher than the predicted one, 
and a negative value means the opposite (for instance, a lockdown has been implemented, 
and the prediction does not reflect this fact)

# EpiInvert 0.2.0 (July 2022)

* We include the management of weekly aggregated incidence data where every 
week a single data is communicated with the accumulated incidence in the last 7 days. 

* We include EpiInvertForecast, a procedure for the short time forecast of the 
restored incidence curve.

# EpiInvert 0.1.1 (May 2022)

* README.md has been simplified to avoid problems with the render of mathematical 
formulas.

* The EpiInvert vignette has been moved outside CRAN repository to avoid problems 
with the render of mathematical formulas.

* Some comments of the function descriptions have been modified.

* None of EpiInvert's functions, examples or tests have been changed in this version.

# EpiInvert 0.1.0 (May 2022)

* First EpiInvert version