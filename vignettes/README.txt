To reduce the time needed to build vignettes during package build
the vignette code is saved in a .orig file.  This code is 
executed by a knit function in the script.R file.  It outputs
a .Rmd file with all of the R computing embedded. The .Rmd
file can be opened and knit applied a second time to create
a html/pdf output.

The .orig will store any graphical output in a directory named
'figure'.  When editing is complete, this folder should be re-named 
to a vignette-specific name (still in the vignette folder) and
references to 'figure' in the .Rmd file should be changed to the
new name.
