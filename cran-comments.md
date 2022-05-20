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