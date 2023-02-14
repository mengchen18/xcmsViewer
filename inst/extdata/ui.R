library(xcmsViewer)
library(omicsViewer)
source("/home/shiny/app/landingPage.R")
# source("landingPage.R")

ui <- fluidPage(
  uiOutput("uis"),
  uiOutput("aout")
)