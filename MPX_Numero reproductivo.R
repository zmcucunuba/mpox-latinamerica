
#--------------------------------------------- 
#
# Script Analisis Descriptivo Tiempo Real
#
#--------------------------------------------- 


#En este archivo se encuentra la estimacion de la transmisibilidad (R) y transmisibilidad variable en el tiempo (Rt) en cinco paises de
#America Latina


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

info <- readRDS("data/cases/info.RDS")



#Fuente de datos LA: https://pubmed.ncbi.nlm.nih.gov/35817231/

#------------Estimacion de la transmisibilidad (R) en Colombia

#Para truncar la incidencia diaria

i_daily_trunc_col <- subset(info$i_daily_col, 
                        from = info$min_date_col, 
                        to = info$max_date_col)

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


R_col

plot(R_col, legend = FALSE) 

#Extraiga la mediana y los intervalos de credibilidad del ( CrI) para el número de reproducción de la siguiente manera:

R_median_col <- R_col$R$`Median(R)`
R_median_col

R_CrI_col <- c(R_col$R$`Quantile.0.025(R)`, R_col$R$`Quantile.0.975(R)`)
R_CrI_col

#------------Estimacion de la transmisibilidad variable en el tiempo (Rt) en Colombia

config_col_1 = make_config(list(mean_si = 9.8, std_si = 9.1))  

# t_start y t_end se configuran automáticamente para estimar R en ventanas deslizantes para 1 semana de forma predeterminada.

# use estimate_R using method = "parametric_si"
Rt_col <- estimate_R(i_daily_trunc_col, method = "parametric_si", 
                 si_data = NULL,
                 config = config_col_1)
Rt_col

# mire las estimaciones de Rt más recientes:
tail(Rt_col$R[, c("t_start", "t_end", "Median(R)", 
              "Quantile.0.025(R)", "Quantile.0.975(R)")])



plot(Rt_col, legend = FALSE)

#------------Estimacion de la transmisibilidad (R)  en Brasil

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

R_median_bra <- R_col$R$`Median(R)`
R_median_bra

R_CrI_bra <- c(R_bra$R$`Quantile.0.025(R)`, R_bra$R$`Quantile.0.975(R)`)
R_CrI_bra



#------------Estimacion de la transmisibilidad variable en el tiempo (Rt) en Brasil

config_bra_1 = make_config(list(mean_si = 9.8, std_si = 9.1))  

# t_start y t_end se configuran automáticamente para estimar R en ventanas deslizantes para 1 semana de forma predeterminada.

# use estimate_R using method = "parametric_si"
Rt_bra <- estimate_R(i_daily_trunc_bra, method = "parametric_si", 
                     si_data = NULL,
                     config = config_bra_1)
Rt_bra

# mire las estimaciones de Rt más recientes:
tail(Rt_bra$R[, c("t_start", "t_end", "Median(R)", 
                  "Quantile.0.025(R)", "Quantile.0.975(R)")])


plot(Rt_bra, legend = FALSE)


#------------Estimacion de la transmisibilidad (R) en Argentina

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

R_median_arg <- R_arg$R$`Median(R)`
R_median_arg

R_CrI_arg <- c(R_arg$R$`Quantile.0.025(R)`, R_bra$R$`Quantile.0.975(R)`)
R_CrI_arg



#------------Estimacion de la transmisibilidad variable en el tiempo (Rt) en Argentina

config_arg_1 = make_config(list(mean_si = 9.8, std_si = 9.1))  

# t_start y t_end se configuran automáticamente para estimar R en ventanas deslizantes para 1 semana de forma predeterminada.

# use estimate_R using method = "parametric_si"
Rt_arg <- estimate_R(i_daily_trunc_arg, method = "parametric_si", 
                     si_data = NULL,
                     config = config_arg_1)
Rt_arg

# mire las estimaciones de Rt más recientes:
tail(Rt_arg$R[, c("t_start", "t_end", "Median(R)", 
                  "Quantile.0.025(R)", "Quantile.0.975(R)")])


plot(Rt_arg, legend = FALSE)


#------------Estimacion de la transmisibilidad (R)  en Peru

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

R_median_per <- R_per$R$`Median(R)`
R_median_per

R_CrI_per <- c(R_per$R$`Quantile.0.025(R)`, R_bra$R$`Quantile.0.975(R)`)
R_CrI_per


#------------Estimacion de la transmisibilidad variable en el tiempo (Rt) en Peru

config_per_1 = make_config(list(mean_si = 9.8, std_si = 9.1))  

# t_start y t_end se configuran automáticamente para estimar R en ventanas deslizantes para 1 semana de forma predeterminada.

# use estimate_R using method = "parametric_si"
Rt_per <- estimate_R(i_daily_trunc_per, method = "parametric_si", 
                     si_data = NULL,
                     config = config_per_1)
Rt_per

# mire las estimaciones de Rt más recientes:
tail(Rt_per$R[, c("t_start", "t_end", "Median(R)", 
                  "Quantile.0.025(R)", "Quantile.0.975(R)")])


plot(Rt_per, legend = FALSE)



#------------Estimacion de la transmisibilidad (R)  en Mexico

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


R_median_mex <- R_mex$R$`Median(R)`
R_median_mex

R_CrI_mex <- c(R_mex$R$`Quantile.0.025(R)`, R_bra$R$`Quantile.0.975(R)`)
R_CrI_mex




#------------Estimacion de la transmisibilidad variable en el tiempo (Rt) en Mexico

config_mex_1 = make_config(list(mean_si = 9.8, std_si = 9.1))  

# t_start y t_end se configuran automáticamente para estimar R en ventanas deslizantes para 1 semana de forma predeterminada.

# use estimate_R using method = "parametric_si"
Rt_mex <- estimate_R(i_daily_trunc_mex, method = "parametric_si", 
                     si_data = NULL,
                     config = config_mex_1)
Rt_mex

# mire las estimaciones de Rt más recientes:
tail(Rt_mex$R[, c("t_start", "t_end", "Median(R)", 
                  "Quantile.0.025(R)", "Quantile.0.975(R)")])


plot(Rt_mex, legend = FALSE)





