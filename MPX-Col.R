#--------------------------------------------- 
#
# Script Analisis Descriptivo Tiempo Real Colombia
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

#Leer datos

url <- "https://www.datos.gov.co/api/views/tmet-yeek/rows.csv?accessType=DOWNLOAD&bom=true&format=true"
col <- read_csv(url) #de la libreria readr


#------ Limpiar los datos

names(col) <- epitrix::clean_labels(names(col))

#------- Ajustando  formato de fechas a Y/M/D

library(lubridate)
col$f_i_s<- dmy(col$fecha_de_inicio_de_sintomas)
tail(table(col$f_i_s))

col$f_d<- dmy(col$fecha_diagnostico)
tail(table(col$f_d))

col$f_n<- dmy(col$fecha_notificacion)
tail(table(col$f_n))

col$f_e<- dmy(col$fecha_de_exantema)
tail(table(col$f_e))

col$f_t_s<- dmy(col$fecha_de_terminacion_del_seguimiento)
tail(table(col$f_t_s))

#-----Corregir fechas (eliminar las fechas de inicio de sintomas que sean mayores o iguales a la fecha de diagnostico)

mistakes <- which(col$f_i_s >= col$f_d)
col[mistakes, ]

col_clean <- col[-mistakes, ]

#--------Curvas de incidencia diaria de Colombia desde enero 2022

i_daily_col <- incidence(col_clean$f_i_s)

plot(i_daily_col, border = "black")


#--------Curvas de incidencia semanal de Colombia desde enero 2022

i_weekly_col <- incidence(col_clean$f_i_s,interval=7)

plot(i_weekly_col,border = "black")

#-------Estimación de la tasa de crecimiento de Colombia mediante un modelo log-lineal 

#Grafico de la incidencia transformada logaritmicamente

ggplot(as.data.frame(i_weekly_col)) + 
  geom_point(aes(x = dates, y = log(counts))) + 
  scale_x_incidence(i_weekly_col) +
  xlab("date") +
  ylab("log weekly incidence Colombia") + 
  theme_minimal()

#Ajuste un modelo log-lineal a los datos de incidencia semanal

f_col <- incidence::fit(i_weekly_col)

plot(i_weekly_col, fit = f_col)

#-------Encontrando una fecha límite adecuada para el modelo log-lineal de Colombia, en función de los retrasos observados

summary(as.numeric(col_clean$f_d - col_clean$f_i_s))

# Semanas a descartar al final de la epicurva
n_weeks_to_discard_col <- 3

min_date_col<- min(i_daily_col$dates)

max_date_col <- max(i_daily_col$dates) - n_weeks_to_discard_col * 7

# Para truncar la incidencia semanal 
i_weekly_trunc_col <- subset(i_weekly_col, 
                             from = min_date_col, 
                             to = max_date_col) # descarte las últimas semanas de datos

#Vuelva a montar y a graficar el modelo logarítmico lineal, pero utilizando los datos truncados i_weekly_trunc. 
#Los resultados deben ser como los siquientes:

f_col<- incidence::fit(i_weekly_trunc_col)
f_col

plot(i_weekly_trunc_col, fit = f_col)

#Observe las estadísticas resumidas de su ajuste:
summary(f_col$model)

#¿El modelo se ajusta bien a los datos?
adjRsq_model_fit_col <- summary(f_col$model)$adj.r.squared

#¿Cuál es la tasa de crecimiento estimada de la epidemia?
daily_growth_rate_col <- f_col$model$coefficients['dates.x']
daily_growth_rate_CI_col <- confint(f_col$model, 'dates.x', level=0.95)


#¿Cuál es el tiempo de duplicación de la epidemia?
doubling_time_days_col <- log(2) / daily_growth_rate_col
doubling_time_days_CI_col <- log(2) / rev(daily_growth_rate_CI_col)

#------------Estimacion de la transmisibilidad (R) en Colombia

#Para truncar la incidencia diaria
i_daily_trunc_col <- subset(i_daily_col, 
                            from = min_date_col, 
                            to = max_date_col)

i_daily_trunc_col$counts

config_col <- make_config(mean_si = 9.8, # media de la distribución 
                          std_si = 9.1,  # desviación estándar de la distribución 
                          t_start = 2,         # día de inicio de la ventana de tiempo
                          t_end = length(i_daily_trunc_col$counts)) # último día de la ventana de tiempo


R_col <- estimate_R(incid = i_daily_trunc_col,
                    method = c("parametric_si"),
                    si_data = NULL,
                    si_sample = NULL,
                    config = config_col)

plot(R_col, legend = FALSE) 

#Extraiga la mediana y los intervalos de credibilidad del ( CrI) para el número de reproducción de la siguiente manera:

R_median_col <- R_col$R$`Median(R)`

R_CrI_col <- c(R_col$R$`Quantile.0.025(R)`, R_col$R$`Quantile.0.975(R)`)

#------------Estimacion de la transmisibilidad variable en el tiempo (Rt) en Colombia

config_col_1 = make_config(list(mean_si = 9.8, std_si = 9.1))  

Rt_col <- estimate_R(i_daily_trunc_col, method = "parametric_si", 
                     si_data = NULL,
                     config = config_col_1)

tail(Rt_col$R[, c("t_start", "t_end", "Median(R)", 
                  "Quantile.0.025(R)", "Quantile.0.975(R)")])

plot(Rt_col, legend = FALSE)

