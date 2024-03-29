## submission EpiInvert package version 0.3.1

* I have received a message from the CRAN staff pointed out  an error running 
‘testthat.R’ in "r-release-macos-arm64" and "r-oldrel-macos-arm64" in version 
0.3.0  due to  differences in the quality of floating point arithmetic operations. 
I have amplified the  margin of error in the estimated number to avoid this problem. 

* I had not detected this problem when using rhub::check() with macOS 10.13.6
compilers. 

## update 1 of EpiInvert package version 0.3.0 submission

Following the suggestions of the CRAN staff the Description field starts now 
with a capital letter and we have included in the Description some references
about our work. 

## submission EpiInvert package version 0.3.0

This is a submission of a new version of the EpiInvert package. The modifications in this new version are: 

* We include EpiIndicators, a procedure for the estimation of the delay and ratio between epidemiological indicators.

### R CMD check results
* 0 errors ✔ | 0 warnings ✔ | 0 notes ✔

### devtools::check_win_devel() results: 
Status: success (1 NOTE) : 

* Possibly misspelled words in DESCRIPTION:
  Variational (4:8)
  variational (11:20)

### rhub::check()

Option 1 (debian-clang-devel) : OK

Option 9 (macos-highsierra-release) : OK

Option 10 (macos-highsierra-release-cran)) : OK

Option 19 (windows-x86_64-release) : OK





## update 3 of EpiInvert package version 0.2.1 submission

Following the suggestions of the CRAN staff we have fixed the test to work with the 
Mac Compilers. 

### rhub::check()

Option 9 (macos-highsierra-release) : OK 
Option 10 (macos-highsierra-release-cran)) : OK


## update 2 of EpiInvert package version 0.2.1 submission

Following the suggestions of the CRAN staff we have corrected the title to take 
into account the "title case" format. 

## update 1 of EpiInvert package version 0.2.1 submission

Following the suggestions of the CRAN staff we have reduced the 
length of the title to less than 65 characters and we wrote package names, 
software names and API (application programming interface) names in single quotes 
in title and description



## submission EpiInvert package version 0.2.1

This is a submission of a new version of the EpiInvert package. The modifications in this new version are: 

* In the EpiInvertForecast we include the option of using the median, instead of
a weighted mean of the incidence curve dataset. 

* In the EpiInvertForecast we add a "trend sentiment" parameter to 
include "a priori" expectations about the incidence trend. The default value is 
"neutral" which means that you do not have any reason to believe that the prediction 
will be different to the initially predicted by the method. A positive "trend sentiment" 
means that you believe that the true value will be higher than the predicted one, 
and a negative value means the opposite (for instance, a lockdown has been implemented, 
and the prediction does not reflect this fact)

### R CMD check results
* 0 errors ✔ | 0 warnings ✔ | 0 notes ✔

### devtools::check_win_devel() results: 
Status: success (1 NOTE) : 

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Luis Alvarez <lalvarez@ulpgc.es>'

### rhub::check()

Option 1 (debian-clang-devel) : OK

Option 10 (macos-highsierra-release-cran)) : OK

Option 19 (windows-x86_64-release) : OK





## submission EpiInvert package version 0.2.0

This is a submission of a new version of the EpiInvert package. The modifications
in this new version are: 

* We include the management of weekly aggregated incidence data where a single 
accumulated incidence value is stored each week. 

* We include EpiInvertForecast, a procedure for the short time forecast of the 
restored incidence curve.

### R CMD check results
* 0 errors ✔ | 0 warnings ✔ | 0 notes ✔

### devtools::check_win_devel() results: 
Status: success (1 NOTE) : 

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Luis Alvarez <lalvarez@ulpgc.es>'


### rhub::check()

Option 1 (debian-clang-devel) : OK

Option 10 (macos-highsierra-release-cran)) : OK

Option 19 (windows-x86_64-release) : OK






## submission EpiInvert package version 0.1.1

This is a submission of a new version of the EpiInvert package. The modifications
in this new version are: 

* README.md has been simplified to avoid problems with the render of mathematical 
formulas.

* The EpiInvert vignette has been moved outside CRAN repository to avoid problems 
with the render of mathematical formulas.

* Some comments of the function descriptions have been modified.

* None of EpiInvert's functions, examples or tests have been changed in this version.

### R CMD check results
* 0 errors ✔ | 0 warnings ✔ | 0 notes ✔

### devtools::check_win_devel() results: 
Status: 1 NOTE : 

Found the following (possibly) invalid URLs:
  URL: https://www.pnas.org/doi/10.1073/pnas.2105112118
    From: README.md
    Status: 503
    Message: Service Unavailable

Found the following (possibly) invalid DOIs:
  DOI: 10.1073/pnas.2105112118
    From: DESCRIPTION
    Status: Service Unavailable
    Message: 503
    
I checked the URL and DOI number and both are correct. 

### rhub::check() results:
* (7) fedora-gcc-devel : OK
* (11) macos-highsierra-release-cran : OK
* (18) windows-x86_64-devel : OK