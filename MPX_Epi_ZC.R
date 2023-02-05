#--------------------------------------------------------------------------------------------------------------------------------------------------- 
#
# Epidemiological findings and estimates of the instantaneous reproduction number and control strategies of the first Mpox outbreak in Latin America
#
#--------------------------------------------------------------------------------------------------------------------------------------------------- 

rm(list=ls())

library(tidyverse)
library(readr)
library(readxl)
library(incidence)
library(distcrete)
library(epitrix)
library(EpiEstim)
library(projections)
library(ggplot2)
library(magrittr)
library(binom)
library(ape)
library(outbreaker2)
library(here)
library(knitr)

#------ Reading data
# Source: PAHO dashboard https://shiny.pahobra.org/monkeypox/

# Agregated cases
# acases <-"https://shiny.pahobra.org/monkeypox/session/f8b64fab92128c4787a68b467d97468f/download/data1?w="
# acases <- read_csv(acases) #de la libreria readr

# Individual cases
# icases  <- "https://shiny.pahobra.org/monkeypox/session/f8b64fab92128c4787a68b467d97468f/download/data2?w="
# icases <-read_csv(icases) #de la libreria readr


acases <- read_csv("data/cases/mpx_data_1811.csv")
icases <- read_csv("data/cases/mpx_linelist_1811.csv")

#------ Data cleaning 
names(acases) <- epitrix::clean_labels(names(acases))
names(icases) <- epitrix::clean_labels(names(icases))


#Creating function for fixing dates
fix_date <- function(dtset, wrong_col){ 
  mistakes <- which(as.Date(wrong_col)<=as.Date("2022-05-01") ) 
  icases_clean <- dtset[-mistakes, ] 
  
  return(icases_clean)
  
}

#Running the function fix_date for each date
icases_clean <- fix_date(icases,icases$date_onset)

icases_clean <- fix_date(icases_clean,icases_clean$date_diagnosis)

icases_clean <- fix_date(icases_clean,icases_clean$report_date)

icases_clean <- fix_date(icases_clean,icases_clean$date_rash)

#Removing the cases in which date_onset >= date of diagnosis
mistakes <- which(icases_clean$date_onset >= icases_clean$date_diagnosis) 
icases_clean[mistakes, ] 
icases_clean <- icases_clean[-mistakes, ] 

#Filtering the aggregated cases by infection status (Confirmed)
acases <- acases %>% filter(classification == "Confirmed")

#Selection of the countries
selected_countries <- c("COL", "BRA", "PER", "MEX", "ARG",  "CHL")

#------ Cummulative cases for selected countries 

c_cum <- 
  ggplot(data = acases %>% filter (iso3 %in% selected_countries)) +
  geom_line(aes(x= date, y = cases), colour = "#4B0082", size = 1) +
  geom_area(aes(x= date, y = cases), fill = "#E6E6FA") +
  facet_wrap(~ iso3, scales = "free_y") +
  ylab(label = "Cumulative confirmed cases\nof Mpox by country") +
  xlab ("Date of report") +
  theme_clean()  

c_cum

library(ggthemes)
#------ New cases for selected countries  (COL, BRA, PER, MEX, ARG)

c_new <-
  ggplot(data = acases %>% filter (iso3 %in% selected_countries)) +
  geom_line(aes(x= date, y = new_cases), colour = "#4B0082", size = 1) +
  geom_area(aes(x= date, y = new_cases), fill = "#E6E6FA") +
  facet_wrap(~ iso3, scales = "free_y") +
  xlab ("Date of report 2022") +
  theme_clean()  +
  ylab(label = "New confirmed cases\nof Mpox by country") 
  

cowplot::plot_grid(c_cum, c_new, nrow = 2, labels = "AUTO")

#------ Cases by country according to date_onset and reporte_date

# lscountries <- c("PER","BRA","MEX","ARG","COL","CHL")

f_plot_epicurve <- function(country) {
  
  report_date_cases <- icases_clean %>% 
    filter (iso3 %in% country) %>%
    group_by(report_date) %>%
    summarise (cases = n())
  
  date_onset_cases <- icases_clean %>% 
    filter (iso3%in% country) %>%
    group_by(date_onset) %>%
    summarise (cases = n())
  
  plot <- ggplot() +
    geom_col(data = report_date_cases, aes(x= report_date, y = cases, fill = "date_report"), 
             alpha = .5, colour = "black", size = 0.1) +
    geom_col(data = date_onset_cases, aes(x= date_onset, y = cases, fill = "date_onset"),
             alpha = .5, colour = "black",  size = 0.1) +
    ggtitle(label = paste(country)) +
    theme_light() +
    guides(fill = guide_legend(title="")) +
    theme(legend.position="none") 
  
  # 
  return(plot)
  
}

cowplot::plot_grid(f_plot_epicurve ("BRA"),
                   f_plot_epicurve ("PER"),
                   f_plot_epicurve ("COL"),
                   f_plot_epicurve ("MEX"),
                   f_plot_epicurve ("ARG"),
                   f_plot_epicurve ("CHL")
                   
                   + 
                     theme(legend.position=c(0.3,.95), 
                           legend.background = element_blank()),
                   nrow = 2)

#------------Estimation of the instantaneous reproduction number R(t)

#To truncate the daily incidence

fRunRt <- function(icases_clean, country_iso, time_start_delay = 2) {
  
  
  i_weekly <- incidence((icases_clean %>% filter (iso3 %in% country_iso))$date_onset,interval=7)
  i_daily <- incidence((icases_clean %>% filter (iso3 %in% country_iso))$date_onset)
  n_weeks_to_discard <- 2
  min_date <- min(i_daily$dates)
  max_date <- max(i_daily$dates) - n_weeks_to_discard * 7
  window <- 1 #in weeks
  
# To truncate the weekly incidence 
  i_weekly_trunc <- subset(i_weekly, 
                           from = min_date, 
                           to = max_date) 
  
# time_start_delay <- 2  
  config <- make_config(mean_prior = 1.5/7,
                        mean_si = 8.5/7, # media de la distribución
                        std_si = 5.0/7,  # desviación estándar de la distribución
                        t_start = time_start_delay:(length(i_weekly_trunc$counts)-window),         # día de inicio de la ventana de tiempo
                        t_end = (time_start_delay + 1):length(i_weekly_trunc$counts)
  ) # último día de la ventana de tiempo
  
# use estimate_R using method = "parametric_si"
  
  data_incidence <- data.frame (I = i_weekly_trunc$counts, 
                                date = as.Date(i_weekly_trunc$weeks), 
                                week = i_weekly_trunc$weeks) 
  data_incidence$t_start <- 1: nrow(data_incidence)
  
  Rt <- estimate_R(data_incidence, 
                   method = "parametric_si", 
                   config = config)
  
  Rt_out <- data.frame(t_start = Rt$R$t_start,
                       t_end = Rt$R$t_end,
                       Rt_mean = Rt$R$`Mean(R)`, 
                       Rt_lower = Rt$R$`Quantile.0.025(R)`,
                       Rt_upper = Rt$R$`Quantile.0.975(R)`)
  
  data_incidence <- merge(data_incidence, Rt_out, by = "t_start", all = TRUE) %>%
    mutate(country = country_iso)
  
  return(data_incidence)
  
}

rt <- rbind(
  fRunRt(icases_clean, "BRA", time_start_delay = 2),
  fRunRt(icases_clean, "ARG", time_start_delay = 2),
  fRunRt(icases_clean, "CHL", time_start_delay = 2),
  fRunRt(icases_clean, "MEX", time_start_delay = 2),
  fRunRt(icases_clean, "PER", time_start_delay = 2),
  fRunRt(icases_clean, "COL", time_start_delay = 2)
)

plot_f <- function(rt_data) {
  my_breaks <- as.Date(c("2022-06-01", "2022-07-01", "2022-08-01",
                 "2022-09-01", "2022-10-01", "2022-11-01", 
                 "2022-12-01"))
  
  rt_plot <- ggplot(rt_data) +
    geom_ribbon(aes(x = date, ymin= Rt_lower, ymax=Rt_upper), fill = "grey") +
    geom_line( aes(x = date, y = Rt_mean), colour = "#4B0082", size = 1) +
    facet_wrap(~ country, scales = "free_y") +
    theme_clean() + ylab ("") + xlab("") +
    coord_cartesian(ylim = c(0, 5)) +
    scale_x_date(limits = c(as.Date("2022-06-01"), as.Date("2022-12-01")),
                 breaks = my_breaks) 
    # geom_area(aes(x= date, y = new_cases), fill = "#E6E6FA") +
  
  inc_plot  <- ggplot(rt_data) +
    geom_col( aes(x = date, y = I), color = "#4B0082") +
    facet_wrap(~ country, scales = "free_y") +
    theme_clean()  + xlab ("") + ylab("") +
    scale_x_date(limits = c(as.Date("2022-06-01"), as.Date("2022-12-01")),
                 breaks = my_breaks)
  
  plot_final <- cowplot::plot_grid(inc_plot, rt_plot, nrow = 2, align = TRUE)
  
  return(plot_final)
  
}

BRA <- plot_f (rt_data = fRunRt(icases_clean, "BRA", time_start_delay = 2))
COL <- plot_f (rt_data = fRunRt(icases_clean, "COL", time_start_delay = 2))
ARG <- plot_f (rt_data = fRunRt(icases_clean, "ARG", time_start_delay = 2))
CHL <- plot_f (rt_data = fRunRt(icases_clean, "CHL", time_start_delay = 2))
PER <- plot_f (rt_data = fRunRt(icases_clean, "PER", time_start_delay = 2))
MEX <- plot_f (rt_data = fRunRt(icases_clean, "MEX", time_start_delay = 2))

cowplot::plot_grid(ARG, BRA, CHL, COL, MEX, PER, nrow = 2)

rt %>% filter(country == "CHL") 
rt %>% filter(country == "COL") 
rt %>% filter(country == "BRA") 
rt %>% filter(country == "MEX") 
rt %>% filter(country == "PER") 
rt %>% filter(country == "ARG")

rt %>% filter(country == "CHL") %>% summarise (max(Rt_mean, na.rm = TRUE))
rt %>% filter(country == "COL") %>% summarise (max(Rt_mean, na.rm = TRUE))
rt %>% filter(country == "BRA") %>% summarise (max(Rt_mean, na.rm = TRUE))
rt %>% filter(country == "MEX") %>% summarise (max(Rt_mean, na.rm = TRUE))
rt %>% filter(country == "PER") %>% summarise (max(Rt_mean, na.rm = TRUE))
rt %>% filter(country == "ARG") %>% summarise (max(Rt_mean, na.rm = TRUE))

rt %>% filter(country == "CHL") %>% summarise (min(Rt_mean, na.rm = TRUE))
rt %>% filter(country == "COL") %>% summarise (min(Rt_mean, na.rm = TRUE))
rt %>% filter(country == "BRA") %>% summarise (min(Rt_mean, na.rm = TRUE))
rt %>% filter(country == "MEX") %>% summarise (min(Rt_mean, na.rm = TRUE))
rt %>% filter(country == "PER") %>% summarise (min(Rt_mean, na.rm = TRUE))
rt %>% filter(country == "ARG") %>% summarise (min(Rt_mean, na.rm = TRUE))

#------------Sensitivity analysis

#Using Prior mean R(t) = 3.8

# cfRunRt <- function(icases_clean, country_iso, time_start_delay = 2) {
#   
#   i_weekly <- incidence((icases_clean %>% filter (iso3 %in% country_iso))$date_onset,interval=7)
#   i_daily <- incidence((icases_clean %>% filter (iso3 %in% country_iso))$date_onset)
#   n_weeks_to_discard <- 2
#   min_date <- min(i_daily$dates)
#   max_date <- max(i_daily$dates) - n_weeks_to_discard * 7
#   window <- 1 #in weeks
#   
#   # To truncate the weekly incidence 
#   i_weekly_trunc <- subset(i_weekly, 
#                            from = min_date, 
#                            to = max_date) 
#   
#   # time_start_delay <- 2  
#   config_2 <- make_config(mean_prior = 3.8/7,
#                         mean_si = 8.5/7, # media de la distribución
#                         std_si = 5.0/7,  # desviación estándar de la distribución
#                         t_start = time_start_delay:(length(i_weekly_trunc$counts)-window),         # día de inicio de la ventana de tiempo
#                         t_end = (time_start_delay + 1):length(i_weekly_trunc$counts)
#   ) # último día de la ventana de tiempo
#   
#   # use estimate_R using method = "parametric_si"
#   
#   data_incidence <- data.frame (I = i_weekly_trunc$counts, 
#                                 date = as.Date(i_weekly_trunc$weeks), 
#                                 week = i_weekly_trunc$weeks) 
#   data_incidence$t_start <- 1: nrow(data_incidence)
#   
#   Rt <- estimate_R(data_incidence, 
#                    method = "parametric_si", 
#                    config = config_2)
#   
#   Rt_out <- data.frame(t_start = Rt$R$t_start,
#                        t_end = Rt$R$t_end,
#                        Rt_mean = Rt$R$`Mean(R)`, 
#                        Rt_lower = Rt$R$`Quantile.0.025(R)`,
#                        Rt_upper = Rt$R$`Quantile.0.975(R)`)
#   
#   data_incidence <- merge(data_incidence, Rt_out, by = "t_start", all = TRUE) %>%
#     mutate(country = country_iso)
#   
#   return(data_incidence)
#   
# }
# 
# rt <- rbind(
#   fRunRt(icases_clean, "BRA", time_start_delay = 2),
#   fRunRt(icases_clean, "ARG", time_start_delay = 2),
#   fRunRt(icases_clean, "CHL", time_start_delay = 2),
#   fRunRt(icases_clean, "MEX", time_start_delay = 2),
#   fRunRt(icases_clean, "PER", time_start_delay = 2),
#   fRunRt(icases_clean, "COL", time_start_delay = 2)
# )
# 
# plot_f <- function(rt_data) {
#   
#   rt_plot <- ggplot(rt_data) +
#     geom_ribbon(aes(x = date, ymin= Rt_lower, ymax=Rt_upper), fill = "grey") +
#     geom_line( aes(x = date, y = Rt_mean), color = "black") +
#     facet_wrap(~ country, scales = "free_y") +
#     theme_bw() + ylab ("") + xlab("") +
#     coord_cartesian(ylim = c(0, 5))
#   
#   inc_plot  <- ggplot(rt_data) +
#     geom_col( aes(x = date, y = I), color = "black") +
#     facet_wrap(~ country, scales = "free_y") +
#     theme_bw()  + xlab ("") + ylab("")
#   
#   plot_final <- cowplot::plot_grid(inc_plot, rt_plot, nrow = 2, align = TRUE)
#   
#   return(plot_final)
#   
# }
# 
# BRA <- plot_f (rt_data = fRunRt(icases_clean, "BRA", time_start_delay = 2))
# COL <- plot_f (rt_data = fRunRt(icases_clean, "COL", time_start_delay = 2))
# ARG <- plot_f (rt_data = fRunRt(icases_clean, "ARG", time_start_delay = 2))
# CHL <- plot_f (rt_data = fRunRt(icases_clean, "CHL", time_start_delay = 2))
# PER <- plot_f (rt_data = fRunRt(icases_clean, "PER", time_start_delay = 2))
# MEX <- plot_f (rt_data = fRunRt(icases_clean, "MEX", time_start_delay = 2))
# 
# cowplot::plot_grid(ARG, BRA, CHL, COL, MEX, PER, nrow = 2)
# 
# rt %>% filter(country == "CHL") 
# rt %>% filter(country == "COL") 
# rt %>% filter(country == "BRA") 
# rt %>% filter(country == "MEX") 
# rt %>% filter(country == "PER") 
# rt %>% filter(country == "ARG")
# 
# rt %>% filter(country == "CHL") %>% summarise (max(Rt_mean, na.rm = TRUE))
# rt %>% filter(country == "COL") %>% summarise (max(Rt_mean, na.rm = TRUE))
# rt %>% filter(country == "BRA") %>% summarise (max(Rt_mean, na.rm = TRUE))
# rt %>% filter(country == "MEX") %>% summarise (max(Rt_mean, na.rm = TRUE))
# rt %>% filter(country == "PER") %>% summarise (max(Rt_mean, na.rm = TRUE))
# rt %>% filter(country == "ARG") %>% summarise (max(Rt_mean, na.rm = TRUE))
# 
# rt %>% filter(country == "CHL") %>% summarise (min(Rt_mean, na.rm = TRUE))
# rt %>% filter(country == "COL") %>% summarise (min(Rt_mean, na.rm = TRUE))
# rt %>% filter(country == "BRA") %>% summarise (min(Rt_mean, na.rm = TRUE))
# rt %>% filter(country == "MEX") %>% summarise (min(Rt_mean, na.rm = TRUE))
# rt %>% filter(country == "PER") %>% summarise (min(Rt_mean, na.rm = TRUE))
# rt %>% filter(country == "ARG") %>% summarise (min(Rt_mean, na.rm = TRUE))


#------------Descriptive analysis 

#proportion of male
icases_clean %>% filter(country_en=="Brazil") %>% count(gender=="MALE") #91.8 
icases_clean %>% filter(country_en=="Argentina") %>% count(gender=="MALE") # 98.5 
icases_clean %>% filter(country_en=="Chile") %>% count(gender=="MALE") # 98.6
icases_clean %>% filter(country_en=="Colombia") %>% count(gender=="MALE") # 97.1
icases_clean %>% filter(country_en=="Mexico") %>% count(gender=="MALE") # 97.5
icases_clean %>% filter(country_en=="Peru") %>% count(gender=="MALE") #98.2

#median age
Bra <- icases_clean %>% filter(country_en=="Brazil") 
median(Bra$age_years, na.rm=TRUE)

Arg <- icases_clean %>% filter(country_en=="Argentina") 
median(Arg$age_years, na.rm=TRUE)

Chl <- icases_clean %>% filter(country_en=="Chile") 
median(Chl$age_years, na.rm=TRUE)

Col <- icases_clean %>% filter(country_en=="Colombia")
median(Col$age_years, na.rm=TRUE)

Mexico <- icases_clean %>% filter(country_en=="Mexico")
median(Mexico$age_years, na.rm=TRUE)

Per <- icases_clean %>% filter(country_en=="Peru") 
median(Per$age_years, na.rm=TRUE)


#------------Map of mpox confirmed cases by country

library(ggplot2)             
library(tidyverse)   
library(dplyr)

mapdata <- map_data("world") ##ggplot2
View(mapdata)

icases_clean_map <- icases_clean %>% rename(region=country_en) 
icases_clean_map %>% count(region) -> casos

mapdata <- left_join(mapdata, casos, by="region")

# mapdata1<-mapdata %>% filter(region %in% c("Colombia","Peru","Chile","Argentina","Brazil","Mexico"))
# View(mapdata1)

map1<-ggplot(mapdata, aes( x = long, y = lat, group=group)) +
  geom_polygon(aes(fill = n), color = "white")
map1

latlimits <- c(-55, 30) 
longlimits <- c(-115, -35) 

map2 <- map1 + 
  # scale_fill_gradient(name = "Cumulative Mpox\nconfirmed cases by country", low = "light blue", high =  "dark blue", na.value = "grey50")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        rect = element_blank()) +
  coord_cartesian(xlim = longlimits, ylim = latlimits) +
  scale_fill_viridis_c( begin = 0, end = 0.5, direction = -1,
                        name = "Cumulative Mpox\nconfirmed cases by country")+
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))


map2




