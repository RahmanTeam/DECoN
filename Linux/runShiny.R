packrat::on()
print("BEGIN runShiny.R")

library(methods)
library(labeling)
library(shiny)
library(R.utils)

args=commandArgs(asValues=TRUE)

Data<-args$Rdata
cnvs <- args$cnvs

file.copy(Data, "shinyGUI/Data.RData",overwrite=TRUE)
file.copy(cnvs, "shinyGUI/cnvs.csv",overwrite=TRUE)
runApp("shinyGUI",launch.browser=F, host='0.0.0.0', port=5888)

print("END runShiny.R")
