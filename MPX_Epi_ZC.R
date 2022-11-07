
#--------------------------------------------- 
#
# Descriptive analysis of MPX outbreak in Latin America 
#
#--------------------------------------------- 

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


acases <- read_csv("data/cases/mpx_data_06_Nov.csv")
icases <- read_csv("data/cases/mpx_linelist_06_ Nov.csv")

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
  geom_line(aes(x= date, y = cases)) +
  facet_wrap(~ iso3, scales = "free_y") +
  ylab(label = "Cumulative confirmed cases\nof monkeypox by country") +
  xlab ("Date")

c_cum


#------ New cases for selected countries  (COL, BRA, PER, MEX, ARG)

c_new <-
  ggplot(data = acases %>% filter (iso3 %in% selected_countries)) +
  geom_line(aes(x= date, y = new_cases)) +
  facet_wrap(~ iso3, scales = "free_y") +
  ylab(label = "New confirmed cases\nof monkeypox by country") +
  xlab ("Date")


cowplot::plot_grid(c_cum, c_new, nrow = 2, labels = "AUTO")




#------ Cases by country according to date_onset and reporte_date

# lscountries <- c("PER","BRA","MEX","ARG","COL", "CHL", "ECU")

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



#------------Estimacion de la transmisibilidad variable en el tiempo (Rt) 

#Para truncar la incidencia diaria

fRunRt <- function(icases_clean, country_iso, time_start_delay = 2) {
  
  
  i_weekly <- incidence((icases_clean %>% filter (iso3 %in% country_iso))$report_date,interval=7)
  i_daily <- incidence((icases_clean %>% filter (iso3 %in% country_iso))$report_date)
  n_weeks_to_discard <- 2
  min_date <- min(i_daily$dates)
  max_date <- max(i_daily$dates) - n_weeks_to_discard * 7
  
  # To truncate the weekly incidence 
  i_weekly_trunc <- subset(i_weekly, 
                           from = min_date, 
                           to = max_date) 
  
  # time_start_delay <- 2  
  config <- make_config(mean_si = 9.8/7, # media de la distribución
                        std_si = 9.1/7,  # desviación estándar de la distribución
                        t_start = time_start_delay:(length(i_weekly_trunc$counts)-1),         # día de inicio de la ventana de tiempo
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
  fRunRt(icases_clean, "ARG", time_start_delay = 3),
  fRunRt(icases_clean, "CHL", time_start_delay = 2),
  fRunRt(icases_clean, "MEX", time_start_delay = 2),
  fRunRt(icases_clean, "PER", time_start_delay = 2),
  fRunRt(icases_clean, "COL", time_start_delay = 2)
)


plot_f <- function(rt_data) {
  
  
  rt_plot <- ggplot(rt_data) +
    geom_ribbon(aes(x = date, ymin= Rt_lower, ymax=Rt_upper), fill = "grey") +
    geom_line( aes(x = date, y = Rt_mean), color = "black") +
    facet_wrap(~ country, scales = "free_y") +
    theme_bw() + ylab ("") + xlab("") +
    coord_cartesian(ylim = c(0, 5))
  
  
  inc_plot  <- ggplot(rt_data) +
    geom_col( aes(x = date, y = I), color = "black") +
    facet_wrap(~ country, scales = "free_y") +
    theme_bw()  + xlab ("") + ylab("")
  
  plot_final <- cowplot::plot_grid(inc_plot, rt_plot, nrow = 2, align = TRUE)
  
  return(plot_final)
  
  
}


BRA <- plot_f (rt_data = fRunRt(icases_clean, "BRA", time_start_delay = 2))
COL <- plot_f (rt_data = fRunRt(icases_clean, "COL", time_start_delay = 2))
ARG <- plot_f (rt_data = fRunRt(icases_clean, "ARG", time_start_delay = 3))
CHL <- plot_f (rt_data = fRunRt(icases_clean, "CHL", time_start_delay = 2))
PER <- plot_f (rt_data = fRunRt(icases_clean, "PER", time_start_delay = 2))
MEX <- plot_f (rt_data = fRunRt(icases_clean, "MEX", time_start_delay = 2))

cowplot::plot_grid(ARG, BRA, CHL, COL, MEX, PER, nrow = 2)

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



rt %>% filter(country == "CHL") %>% summarise (mean(Rt_mean, na.rm = TRUE))
rt %>% filter(country == "COL") %>% summarise (mean(Rt_mean, na.rm = TRUE))
rt %>% filter(country == "BRA") %>% summarise (mean(Rt_mean, na.rm = TRUE))
rt %>% filter(country == "MEX") %>% summarise (mean(Rt_mean, na.rm = TRUE))
rt %>% filter(country == "PER") %>% summarise (mean(Rt_mean, na.rm = TRUE))
rt %>% filter(country == "ARG") %>% summarise (mean(Rt_mean, na.rm = TRUE))



rt %>% filter(country == "CHL") 
rt %>% filter(country == "COL") 
rt %>% filter(country == "BRA") 
rt %>% filter(country == "MEX") 
rt %>% filter(country == "PER") 
rt %>% filter(country == "ARG")


#numero de casos por millon de habitantes (25-oct-2022)
COL_pm <- 3298/51600000 *1e6
BRA_pm <- 8978/212600000 *1e6
ARG_pm <- 627/45380000 * 1e6
PER_pm <- 2981/32970000 * 1e6
CHL_pm <- 1127/19120000  * 1e6
MEX_pm <- 2468/128900000 * 1e6



#cuantos son hombres
icases_clean %>% filter(country_en=="Brazil") %>% count(gender=="MALE")
icases_clean %>% filter(country_en=="Argentina") %>% count(gender=="MALE")
icases_clean %>% filter(country_en=="Chile") %>% count(gender=="MALE")
icases_clean %>% filter(country_en=="Colombia") %>% count(gender=="MALE")
icases_clean %>% filter(country_en=="Mexico") %>% count(gender=="MALE")
icases_clean %>% filter(country_en=="Peru") %>% count(gender=="MALE")

#cuantos son msm
icases_clean %>% filter(country_en=="Brazil") %>% count(sexual_orientation=="MSM")
icases_clean %>% filter(country_en=="Argentina") %>% count(sexual_orientation=="MSM")
icases_clean %>% filter(country_en=="Chile") %>% count(sexual_orientation=="MSM")
icases_clean %>% filter(country_en=="Colombia") %>% count(sexual_orientation=="MSM")
icases_clean %>% filter(country_en=="Mexico") %>% count(sexual_orientation=="MSM")
icases_clean %>% filter(country_en=="Peru") %>% count(sexual_orientation=="MSM")










