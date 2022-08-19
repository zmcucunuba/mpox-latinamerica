
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
# acases <-"https://shiny.pahobra.org/monkeypox/session/f8b64fab92128c4787a68b467d97468f/download/data1?w="
# acases <- read_csv(acases) #de la libreria readr

# Individual
# icases  <- "https://shiny.pahobra.org/monkeypox/session/f8b64fab92128c4787a68b467d97468f/download/data2?w="
# icases <-read_csv(icases) #de la libreria readr


acases <- read_csv("data/cases/mpx_data.csv")
icases <- read_csv("data/cases/mpx_linelist.csv")


# Cleaning data
names(acases) <- epitrix::clean_labels(names(acases))
names(icases) <- epitrix::clean_labels(names(icases))

acases <- acases %>% filter(classification == "Confirmed")


selected_countries <- c("COL", "BRA", "PER", "BOL", "MEX", "ARG")

ggplot(data = acases %>% filter (iso3 %in% selected_countries)) +
  geom_line(aes(x= date, y = cases)) +
  facet_wrap(~ iso3, scales = "free_y") +
  ggtitle(label = "Cummulative cases")

ggplot(data = acases %>% filter (iso3 %in% selected_countries)) +
  geom_line(aes(x= date, y = new_cases)) +
  facet_wrap(~ iso3, scales = "free_y") +
  ggtitle(label = "New cases")


ggplot(data = acases %>% filter (iso3 %in% "COL")) +
  geom_line(aes(x= date, y = new_cases)) +
  facet_wrap(~ iso3, scales = "free_y") +
  ggtitle(label = "New cases COL")



col_not_cases <- icases %>% 
  filter (iso3 %in% "COL") %>%
  group_by(report_date) %>%
  summarise (cases = n())


col_fis_cases <- icases %>% 
  filter (iso3 %in% "COL") %>%
  group_by(date_onset) %>%
  summarise (cases = n())

ggplot() +
  geom_col(data = col_not_cases, aes(x= report_date, y = cases, fill = "report_date"), 
           alpha = .5, colour = "black") +
  geom_col(data = col_fis_cases, aes(x= date_onset, y = cases, fill = "date_onset"),
           alpha = .5, colour = "black") +
  ggtitle(label = "COL cases") +
  guides(fill = guide_legend(title="")) 



scl_not_cases <- icases %>% 
  filter (iso3 %in% selected_countries) %>%
  group_by(report_date, iso3) %>%
  summarise (cases = n())

sc_fis_cases <- icases %>% 
  filter (iso3 %in% selected_countries) %>%
  group_by(date_onset, iso3) %>%
  summarise (cases = n())

ggplot() +
  geom_col(data = scl_not_cases, aes(x= report_date, y = cases, fill = "report_date"), 
           alpha = .5, colour = "black") +
  geom_col(data = sc_fis_cases, aes(x= date_onset, y = cases, fill = "date_onset"),
           alpha = .5, colour = "black") +
  ggtitle(label = "selected countries ") +
  guides(fill = guide_legend(title=""))  +
  facet_wrap(~ iso3, scales = "free_y") 
  

  





