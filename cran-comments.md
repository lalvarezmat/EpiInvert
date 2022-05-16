## Resubmission
This is a resubmission. In this version I have, following the suggestions of 
Uwe Ligges (CRAN Teams):

* Reduced each example to less than 5 sec. 

Now, when executing devtools::check_win_devel() the EpiInvert-Ex.timings is : 

name	user	system	elapsed
EpiInvert	0.35	0.00	0.36

### R CMD check results
* 0 errors ✔ | 0 warnings ✔ | 0 notes ✔

### devtools::check_win_devel() results: 
Status: 1 NOTE: 

Possibly misspelled words in DESCRIPTION:
  Variational (5:32)
  al (15:14, 16:14)
  et (15:11, 16:11)

Found the following (possibly) invalid URLs:
  URL: https://www.pnas.org/doi/10.1073/pnas.2105112118
    From: inst/doc/Demo.html
          README.md
    Status: 503
    Message: Service Unavailable

Found the following (possibly) invalid DOIs:
  DOI: 10.1073/pnas.2105112118
    From: DESCRIPTION
    Status: Service Unavailable
    Message: 503
    
Concerning the misspellings: The word variational is correct. It refers to the 
mathematical variational methods. The string "et al." is commonly used to reduce 
the number of authors in a cite.

The URL https://www.pnas.org/doi/10.1073/pnas.2105112118 is a valid URL address.

The DOI number <10.1073/pnas.2105112118>  is correct. 
    


### rhub::check() results:

* (1) debian-clang-devel : OK
* (7) fedora-gcc-devel : OK
* (11) macos-highsierra-release-cran : OK
* (18) windows-x86_64-devel : OK
