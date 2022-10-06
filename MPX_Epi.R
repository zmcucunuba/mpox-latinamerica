
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
library(epicontacts)
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


acases <- read_csv("data/cases/mpx_data_23Sep.csv")
icases <- read_csv("data/cases/mpx_linelist_23Sep.csv")

#------ Data cleaning 
names(acases) <- epitrix::clean_labels(names(acases))
names(icases) <- epitrix::clean_labels(names(icases))


#Creating function for fixing dates
fix_date <- function(base, wrong_col){ 
  mistakes <- which(as.Date(wrong_col)<=as.Date("2022-05-01") ) 
  icases_clean <- base[-mistakes, ] 
  
  return(icases_clean)
  
}

#Running the function fix_date for each date
icases_clean <- fix_date(icases,icases$date_onset)

icases_clean<- fix_date(icases_clean,icases_clean$date_diagnosis)

icases_clean<- fix_date(icases_clean,icases_clean$date_rash)

#Removing the cases in which date_onset >= date of diagnosis
mistakes <- which(icases_clean$date_onset >= icases_clean$date_diagnosis) 
icases_clean[mistakes, ] 
icases_clean <- icases_clean[-mistakes, ] 

#Filtering the aggregated cases by infection status (Confirmed)
acases <- acases %>% filter(classification == "Confirmed")

#Selection of the countries
selected_countries <- c("COL", "BRA", "PER", "MEX", "ARG")

#------ Cummulative cases for selected countries (COL, BRA, PER, MEX, ARG)



c_cum <- 
  ggplot(data = acases %>% filter (iso3 %in% selected_countries)) +
  geom_line(aes(x= date, y = cases)) +
  facet_wrap(~ iso3, scales = "free_y") +
  ylab(label = "Casos acumulados\nde viruela símica") +
  xlab ("Fecha")


#------ New cases for selected countries  (COL, BRA, PER, MEX, ARG)

c_new <-
  ggplot(data = acases %>% filter (iso3 %in% selected_countries)) +
  geom_line(aes(x= date, y = new_cases)) +
  facet_wrap(~ iso3, scales = "free_y") +
  ylab(label = "Casos nuevos\nde viruela símica") +
  xlab ("Fecha")


cowplot::plot_grid(c_cum, c_new, nrow = 2, labels = "AUTO")




#------ Cases by country according to date_onset and reporte_date

lscountries <- c("PER","BRA","MEX","ARG","COL")

for(country in lscountries){
  
  report_date_cases <- icases_clean%>% 
    filter (iso3 %in% country) %>%
    group_by(report_date) %>%
    summarise (cases = n())
  
  date_onset_cases <- icases_clean%>% 
    filter (iso3%in% country) %>%
    group_by(date_onset) %>%
    summarise (cases = n())
  
  plot <- ggplot() +
    geom_col(data = report_date_cases, aes(x= report_date, y = cases, fill = "date_report"), 
             alpha = .5, colour = "black") +
    geom_col(data = date_onset_cases, aes(x= date_onset, y = cases, fill = "date_onset"),
             alpha = .5, colour = "black") +
    ggtitle(label = paste(country," cases")) +
    guides(fill = guide_legend(title="")) 
  print(plot)
}


#--------Curves of daily and weekly incidence by country

#Second loop
for(country in lscountries){
  
  #By day
  i_daily <- incidence((icases_clean %>% filter (iso3 %in% country))$date_onset)
  print(plot(i_daily, border = "black"))
  
  #By week
  i_weekly <- incidence((icases_clean %>% filter (iso3 %in% country))$date_onset,interval=7)
  print(plot(i_weekly,border = "black"))
}

#-------Estimation of the growth rate and doubling time by country using a log-linear model 

#Third loop
for(country in lscountries){
  
  #Log-transformed incidence graph
  plot<-  ggplot(as.data.frame(i_weekly)) + 
    geom_point(aes(x = dates, y = log(counts))) + 
    scale_x_incidence(i_weekly) +
    xlab("date") +
    ylab(paste("log weekly incidence ",country)) + 
    theme_minimal()
  print(plot)
}

#Four loop
for(country in lscountries){
  

  #Fitting a log-linear model to the weekly incidence data
  f <- incidence::fit(i_weekly)
  print(plot(i_weekly, fit = f))
  
  #Weeks to discard at the end of the epicurve 
  n_weeks_to_discard<- 2
  min_date<- min(i_daily$dates)
  max_date <- max(i_daily$dates) - n_weeks_to_discard * 7
  
  # To truncate the weekly incidence 
  i_weekly_trunc <- subset(i_weekly, 
                           from = min_date, 
                           to = max_date) 

  
  #Re-assemble and plot the log-linear model, but using the truncated i_weekly_trunc data. 
  #The results should look like the following:
  
  f<- incidence::fit(i_weekly_trunc)
  print(plot(i_weekly_trunc, fit = f))
  
  #Note the summary statistics of your adjustment:
  summary(f$model)
  
  #Does the model fit the data?
  adjRsq_model_fit <- summary(f$model)$adj.r.squared
  
  #What is the estimated growth rate of the epidemic?
  daily_growth_rate<- f$model$coefficients['dates.x']
  #Confidence interval:
  daily_growth_rate_CI <- confint(f$model, 'dates.x', level=0.95)
  
 #What is the epidemic doubling time?
  doubling_time_days <- log(2) / daily_growth_rate
  #Confidence interval:
  doubling_time_days_CI <- log(2) / rev(daily_growth_rate_CI)
  
}


#Fifth loop

for(country in lscountries){

#Para truncar la incidencia diaria

i_daily_trunc <- subset(i_daily, 
                            from = min_date, 
                            to = max_date)

i_daily_trunc$counts

config <- make_config(mean_si = 9.8, # media de la distribución 
                          std_si = 9.1,  # desviación estándar de la distribución 
                          t_start = 2,         # día de inicio de la ventana de tiempo
                          t_end = length(i_daily_trunc$counts)) # último día de la ventana de tiempo


R <- estimate_R(incid = i_daily_trunc,
                    method = c("parametric_si"),
                    si_data = NULL,
                    si_sample = NULL,
                    config = config)

print(plot(R, legend = FALSE))

#Extraiga la mediana y los intervalos de credibilidad del ( CrI) para el número de reproducción de la siguiente manera:

R_median <- R$R$`Median(R)`


R_CrI <- c(R$R$`Quantile.0.025(R)`, R$R$`Quantile.0.975(R)`)
R_CrI

#------------Estimacion de la transmisibilidad variable en el tiempo (Rt) 

config= make_config(list(mean_si = 9.8, std_si = 9.1))  

# use estimate_R using method = "parametric_si"
Rt <- estimate_R(i_daily_trunc, method = "parametric_si", 
                     si_data = NULL,
                     config = config)

# mire las estimaciones de Rt más recientes:
tail(Rt$R[, c("t_start", "t_end", "Median(R)", 
                  "Quantile.0.025(R)", "Quantile.0.975(R)")])

print(plot(Rt, legend = FALSE))

}


