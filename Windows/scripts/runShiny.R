renv::restore()

print("BEGIN runShiny.R")

library(methods)
library(labeling)
library(shiny)
library(R.utils)

args=commandArgs(TRUE)

Data<-args[1]

file.copy(Data, "scripts/shinyGUI/Data.RData",overwrite=TRUE)

runApp("scripts/shinyGUI",launch.browser=T)


print("END runShiny.R")
