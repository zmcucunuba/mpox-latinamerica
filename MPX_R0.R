
#--------------------------------------------- 
#
# Script Analisis Descriptivo Tiempo Real
#
#--------------------------------------------- 


#En este archivo se encuentra la estimacion  de los numeros reproductivos de la epidemia en cinco paises de America Latina


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

#Fuente de datos: https://pubmed.ncbi.nlm.nih.gov/35817231/

#------------Estimacion del RO en Colombia

#Para truncar la incidencia diaria

i_daily_trunc_col <- subset(i_daily_col, 
                        from = min_date_col, 
                        to = max_date_col)


config_col <- make_config(mean_si = 9.8, # media de la distribución 
                      std_si = 9.1,  # desviación estándar de la distribución 
                      t_start = 2,         # día de inicio de la ventana de tiempo
                      t_end = length(i_daily_trunc_col$counts)) # último día de la ventana de tiempo


R_col <- estimate_R(incid = i_daily_trunc_col,
                method = c("parametric_si"),
                si_data = NULL,
                si_sample = NULL,
                config = config_col)


R_col

plot(R_col, legend = FALSE) 

#Extraiga la mediana y los intervalos de credibilidad del ( CrI) para el número de reproducción de la siguiente manera:

R_median_col <- R_col$R_col$`Median(R_col)`
R_median_col

R_CrI_col <- c(R_col$R_col$`Quantile.0.025(R_col)`, R_col$R_col$`Quantile.0.975(R_col)`)
R_CrI_col

#me sale null


#------------Estimacion del RO en Brasil

#Para truncar la incidencia diaria

i_daily_trunc_bra <- subset(i_daily_bra, 
                            from = min_date_bra, 
                            to = max_date_bra)


config_bra <- make_config(mean_si = 9.8, # media de la distribución 
                          std_si = 9.1,  # desviación estándar de la distribución 
                          t_start = 2,         # día de inicio de la ventana de tiempo
                          t_end = length(i_daily_trunc_bra$counts)) # último día de la ventana de tiempo


R_bra <- estimate_R(incid = i_daily_trunc_bra,
                    method = c("parametric_si"),
                    si_data = NULL,
                    si_sample = NULL,
                    config = config_bra)


R_bra

plot(R_bra, legend = FALSE) 



#------------Estimacion del RO en Argentina

#Para truncar la incidencia diaria

i_daily_trunc_arg <- subset(i_daily_arg, 
                            from = min_date_arg, 
                            to = max_date_arg)


config_arg <- make_config(mean_si = 9.8, # media de la distribución 
                          std_si = 9.1,  # desviación estándar de la distribución 
                          t_start = 2,         # día de inicio de la ventana de tiempo
                          t_end = length(i_daily_trunc_arg$counts)) # último día de la ventana de tiempo


R_arg <- estimate_R(incid = i_daily_trunc_arg,
                    method = c("parametric_si"),
                    si_data = NULL,
                    si_sample = NULL,
                    config = config_arg)

R_arg

plot(R_arg, legend = FALSE) 



#------------Estimacion del RO en Peru

#Para truncar la incidencia diaria

i_daily_trunc_per <- subset(i_daily_per, 
                            from = min_date_per, 
                            to = max_date_per)


config_per <- make_config(mean_si = 9.8, # media de la distribución 
                          std_si = 9.1,  # desviación estándar de la distribución 
                          t_start = 2,         # día de inicio de la ventana de tiempo
                          t_end = length(i_daily_trunc_per$counts)) # último día de la ventana de tiempo


R_per <- estimate_R(incid = i_daily_trunc_per,
                    method = c("parametric_si"),
                    si_data = NULL,
                    si_sample = NULL,
                    config = config_per)


R_per

plot(R_per, legend = FALSE) 



#------------Estimacion del RO en Mexico

#Para truncar la incidencia diaria

i_daily_trunc_mex <- subset(i_daily_mex, 
                            from = min_date_mex, 
                            to = max_date_mex)


config_mex <- make_config(mean_si = 9.8, # media de la distribución 
                          std_si = 9.1,  # desviación estándar de la distribución 
                          t_start = 2,         # día de inicio de la ventana de tiempo
                          t_end = length(i_daily_trunc_mex$counts)) # último día de la ventana de tiempo


R_mex <- estimate_R(incid = i_daily_trunc_mex,
                    method = c("parametric_si"),
                    si_data = NULL,
                    si_sample = NULL,
                    config = config_mex)


R_mex

plot(R_mex, legend = FALSE) 









