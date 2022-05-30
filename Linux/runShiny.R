renv::restore()

print("BEGIN runShiny.R")

library(methods)
library(labeling)
library(shiny)
library(R.utils)
library(optparse)

option_list<-list(make_option('--RData',help='Summary RData file (required)',dest='RData'))

opt<-parse_args(OptionParser(option_list=option_list))

Data<-opt$RData

file.copy(Data, "scripts/shinyGUI/Data.RData",overwrite=TRUE)

runApp("scripts/shinyGUI",launch.browser=T)


print("END runShiny.R")
