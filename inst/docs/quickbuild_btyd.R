# Modeled after newbuild.R
# Run this from the command line with Rscript:
# Rscript quickbuild_btyd.R
library("devtools")

setwd(file.path(Sys.getenv()['HOME'], 'BTYD'))

document()
build()
install(build_vignettes = TRUE)
check()
