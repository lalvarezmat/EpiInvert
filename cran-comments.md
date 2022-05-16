## First submission of the EpiInvert package

### R CMD check results
* 0 errors ✔ | 0 warnings ✔ | 0 notes ✔

### devtools::check_win_devel() results: 
Status: 2 NOTEs : 

* Possibly misspelled words in DESCRIPTION:
  Variational (4:32)
  al (15:14, 16:14)
  et (15:11, 16:11)

* Found the following (possibly) invalid URLs: 
  URL: https://www.pnas.org/doi/10.1073/pnas.2105112118
    From: inst/doc/Demo.html
          README.md
    Status: 503
    Message: Service Unavailable

* Found the following (possibly) invalid DOIs:
  DOI: 10.1073/pnas.2105112118
    From: DESCRIPTION
    Status: Service Unavailable
    Message: 503
    
The first note is due to the fact that the string "et al." used to specify the 
name of the authors was considered a misspelled word

The second note concerns the url/doi address : 
https://www.pnas.org/doi/10.1073/pnas.2105112118

I checked the address and it is correct. 

### rhub::check() results:

* (1) debian-clang-devel : OK
* (7) fedora-gcc-devel : OK
* (11) macos-highsierra-release-cran : OK
* (18) windows-x86_64-devel : OK
