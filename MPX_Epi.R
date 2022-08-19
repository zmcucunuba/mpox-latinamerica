
#--------------------------------------------- 
#
# Script Analisis Descriptivo Tiempo Real
#
#--------------------------------------------- 
rm(list=ls())


library(tidyverse)
library(readr)

#------ Leer datos
# Fuente: PAHO dashboard https://shiny.pahobra.org/monkeypox/

# Agregated
acases <-"https://shiny.pahobra.org/monkeypox/session/f8b64fab92128c4787a68b467d97468f/download/data1?w="
acases <- read_csv(acases) #de la libreria readr

# Individual
icases  <- "https://shiny.pahobra.org/monkeypox/session/f8b64fab92128c4787a68b467d97468f/download/data2?w="
icases <-read_csv(icases) #de la libreria readr
