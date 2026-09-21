setwd("C:/proj/repos/clinDRpack/clinDR/vignettes")


##
## to avoid a note in R cmd check, mv graphic files in figure
## to a new graphic-named folder in the vignette directory 
## and and change 'figure' references in the
## compiled .Rmd file
knitr::knit("DRmodeling.Rmd.orig", output = "DRmodeling.Rmd")
knitr::knit("gettingStarted.Rmd.orig", output = "gettingStarted.Rmd")
knitr::knit("mixturePrior.Rmd.orig", output = "mixturePrior.Rmd")

##Note that installation with Rstudio required creation of tar
## package that is subsequently installed
browseVignettes(package='clinDR')
  
