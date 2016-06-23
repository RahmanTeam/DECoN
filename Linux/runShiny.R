packrat::on()

print("BEGIN runShiny.R")

library(methods)
library(labeling)
library(shiny)
library(R.utils)

args=commandArgs(asValues=TRUE)

Data<-args$Rdata

print(Data)

file.copy(Data, "shinyGUI/Data.RData",overwrite=TRUE)

runApp("shinyGUI",launch.browser=T)


print("END runShiny.R")
