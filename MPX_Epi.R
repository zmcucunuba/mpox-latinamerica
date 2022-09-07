
#--------------------------------------------- 
#
# Script Analisis Descriptivo Tiempo Real
#
#--------------------------------------------- 


#En este archivo se encuentra el calculo y graficas de la incidencia de MPXV en paises de America Latina, asi como la estimacion e interpretacion de
#la tasa de crecimiento y el tiempo de duplicación de la epidemia discriminado por territorio


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

#------ Leer datos
# Fuente: PAHO dashboard https://shiny.pahobra.org/monkeypox/

# Agregated
# acases <-"https://shiny.pahobra.org/monkeypox/session/f8b64fab92128c4787a68b467d97468f/download/data1?w="
# acases <- read_csv(acases) #de la libreria readr

# Individual
# icases  <- "https://shiny.pahobra.org/monkeypox/session/f8b64fab92128c4787a68b467d97468f/download/data2?w="
# icases <-read_csv(icases) #de la libreria readr


acases <- read_csv("data/cases/mpx_data(1sep).csv")
icases <- read_csv("data/cases/mpx_linelist (1sep).csv")

#------ Limpiar los datos


names(acases) <- epitrix::clean_labels(names(acases))
names(icases) <- epitrix::clean_labels(names(icases))

mistakes_onset <- which(icases$date_onset<=as.Date("2021-12-31") ) 
icases_clean_onset <- icases[-mistakes_onset, ] 

#########

corregir_fechas <- function(bd, col_sucia){ 
  mistakes <- which(as.Date(col_sucia)<=as.Date("2021-12-31") ) 
  icases_clean <- bd[-mistakes, ] 
 
  return(icases_clean)
  
}

icases_clean <- corregir_fechas(icases,icases$date_onset)

icases_clean<- corregir_fechas(icases_clean,icases_clean$date_diagnosis)

icases_clean<- corregir_fechas(icases_clean,icases_clean$date_rash)


#######


mistakes <- which(icases_clean$date_onset >= icases_clean$date_diagnosis) 
icases_clean[mistakes, ] #Se remueven observaciones con fechas de inicio de sintomas >= a fecha de diagnostico
icases_clean <- icases_clean[-mistakes, ] 

acases <- acases %>% filter(classification == "Confirmed")
selected_countries <- c("COL", "BRA", "PER", "BOL", "MEX", "ARG")


#------ Grafica de casos acumulados de paises seleccionados (COL, BRA, PER, BOL, MEX, ARG)

ggplot(data = acases %>% filter (iso3 %in% selected_countries)) +
  geom_line(aes(x= date, y = cases)) +
  facet_wrap(~ iso3, scales = "free_y") +
  ggtitle(label = "Cummulative cases")

#------ Grafica de casos nuevos de paises seleccionados  (COL, BRA, PER, BOL, MEX, ARG)

ggplot(data = acases %>% filter (iso3 %in% selected_countries)) +
  geom_line(aes(x= date, y = new_cases)) +
  facet_wrap(~ iso3, scales = "free_y") +
  ggtitle(label = "New cases")

#------ Grafica de casos nuevos de Colombia

ggplot(data = acases %>% filter (iso3 %in% "COL")) +
  geom_line(aes(x= date, y = new_cases)) +
  facet_wrap(~ iso3, scales = "free_y") +
  ggtitle(label = "New cases COL")

#------ Grafica de casos de Colombia discriminados por fecha de inicio de sintomas y fecha de reporte


col_not_cases <- icases_clean%>% 
  filter (iso3 %in% "COL") %>%
  group_by(report_date) %>%
  summarise (cases = n())


col_fis_cases <- icases_clean%>% 
  filter (iso3%in% "COL") %>%
  group_by(date_onset) %>%
  summarise (cases = n())

ggplot() +
  geom_col(data = col_not_cases, aes(x= report_date, y = cases, fill = "date_report"), 
           alpha = .5, colour = "black") +
  geom_col(data = col_fis_cases, aes(x= date_onset, y = cases, fill = "date_onset"),
           alpha = .5, colour = "black") +
  ggtitle(label = "COL cases") +
  guides(fill = guide_legend(title="")) 

#------ Grafica de casos de Peru discriminado por inicio de sintomas y fecha de reporte


per_not_cases <- icases_clean %>% 
  filter (iso3 %in% "PER") %>%
  group_by(report_date) %>%
  summarise (cases = n())


per_fis_cases <- icases_clean %>% 
  filter (iso3 %in% "PER") %>%
  group_by(date_onset) %>%
  summarise (cases = n())

ggplot() +
  geom_col(data = per_not_cases, aes(x= report_date, y = cases, fill = "date_report"), 
           alpha = .5, colour = "black") +
  geom_col(data = per_fis_cases, aes(x= date_onset, y = cases, fill = "date_onset"),
           alpha = .5, colour = "black") +
  ggtitle(label = "PER cases") +
  guides(fill = guide_legend(title="")) 

#------ Grafica de casos de Argentina discriminado por inicio de sintomas y fecha de reporte


arg_not_cases <- icases_clean %>% 
  filter (iso3 %in% "ARG") %>%
  group_by(report_date) %>%
  summarise (cases = n())


arg_fis_cases <- icases_clean %>% 
  filter (iso3 %in% "ARG") %>%
  group_by(date_onset) %>%
  summarise (cases = n())

ggplot() +

  geom_col(data = arg_not_cases, aes(x= report_date, y = cases, fill = "date_report"), 
           alpha = .5, colour = "black") +
  geom_col(data = arg_fis_cases, aes(x= date_onset, y = cases, fill = "date_onset"),
           alpha = .5, colour = "black") +
  ggtitle(label = "ARG cases") +
  guides(fill = guide_legend(title="")) 

  geom_col(data = scl_not_cases, aes(x= report_date, y = cases, fill = "report_date"), 
           alpha = .5, colour = "black", size = 0.1) +
  geom_col(data = sc_fis_cases, aes(x= date_onset, y = cases, fill = "date_onset"),
           alpha = .5, colour = "black", size = 0.1) +
  ggtitle(label = "selected countries ") +
  guides(fill = guide_legend(title=""))  +
  facet_wrap(~ iso3, scales = "free_y") 


  
#------ Grafica de casos de Mexico discriminado por inicio de sintomas y fecha de reporte


mex_not_cases <- icases_clean %>% 
  filter (iso3 %in% "MEX") %>%
  group_by(report_date) %>%
  summarise (cases = n())


mex_fis_cases <- icases_clean %>% 
  filter (iso3 %in% "MEX") %>%
  group_by(date_onset) %>%
  summarise (cases = n())

ggplot() +
  geom_col(data = mex_not_cases, aes(x= report_date, y = cases, fill = "date_report"), 
           alpha = .5, colour = "black") +
  geom_col(data = mex_fis_cases, aes(x= date_onset, y = cases, fill = "date_onset"),
           alpha = .5, colour = "black") +
  ggtitle(label = "MEX cases") +
  guides(fill = guide_legend(title="")) 


#------ Grafica de casos de Brasil discriminado por inicio de sintomas y fecha de reporte 


bra_not_cases <- icases_clean %>% 
  filter (iso3 %in% "BRA") %>%
  group_by(report_date) %>%
  summarise (cases = n())


bra_fis_cases <- icases_clean %>% 
  filter (iso3 %in% "BRA") %>%
  group_by(date_onset) %>%
  summarise (cases = n())

ggplot() +
  geom_col(data = bra_not_cases, aes(x= report_date, y = cases, fill = "date_report"), 
           alpha = .5, colour = "black") +
  geom_col(data = bra_fis_cases, aes(x= date_onset, y = cases, fill = "date_onset"),
           alpha = .5, colour = "black") +
  ggtitle(label = "BRA cases") +
  guides(fill = guide_legend(title="")) 



#--------Curvas de incidencia diaria de Colombia desde enero 2022

i_daily_col <- incidence((icases_clean %>% filter (iso3 %in% "COL"))$date_onset)

i_daily_col



plot(i_daily_col, border = "black")


#--------Curvas de incidencia semanal de Colombia desde enero 2022

i_weekly_col <- incidence((icases_clean %>% filter (iso3 %in% "COL"))$date_onset,interval=7)

i_weekly_col



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
f_col


plot(i_weekly_col, fit = f_col)

#-------Encontrando una fecha límite adecuada para el modelo log-lineal de Colombia, en función de los retrasos observados

summary(as.numeric(icases_clean%>% filter (iso3 %in% "COL")$date_diagnosis - icases_clean%>% filter (iso3 %in% "COL")$date_onset))
#no me corre
#VERIFICAR


#¿cuántas semanas debe descartar al final de la epicurva? 

# Semanas a descartar al final de la epicurva
n_weeks_to_discard_col <- 2

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


#Puede observar la bondad del ajuste (Rsquared), la pendiente estimada (tasa de crecimiento/growth rate) y el tiempo de duplicación correspondiente 
#como se muestra a continuación:

#¿El modelo se ajusta bien a los datos?
adjRsq_model_fit_col <- summary(f_col$model)$adj.r.squared

#¿Cuál es la tasa de crecimiento estimada de la epidemia?
daily_growth_rate_col <- f_col$model$coefficients['dates.x']


# intervalo de confianza:
daily_growth_rate_CI_col <- confint(f_col$model, 'dates.x', level=0.95)


#¿Cuál es el tiempo de duplicación de la epidemia?
doubling_time_days_col <- log(2) / daily_growth_rate_col

# intervalo de confianza:
doubling_time_days_CI_col <- log(2) / rev(daily_growth_rate_CI_col)



#--------Curvas de incidencia diaria de Brasil desde enero 2022 (corre raro)

i_daily_bra <- incidence((icases_clean %>% filter (iso3 %in% "BRA"))$date_onset)

i_daily_bra


plot(i_daily_bra, border = "black")


#--------Curvas de incidencia semanal de Brasil desde enero 2022 (corre raro)

i_weekly_bra<- incidence((icases_clean %>% filter (iso3 %in% "BRA"))$date_onset,interval=7)

i_weekly_bra



plot(i_weekly_bra,border = "black")

#-------Estimación de la tasa de crecimiento de Brasil mediante un modelo log-lineal 

#Grafico de la incidencia transformada logaritmicamente

ggplot(as.data.frame(i_weekly_bra)) + 
  geom_point(aes(x = dates, y = log(counts))) + 
  scale_x_incidence(i_weekly_bra) +
  xlab("date") +
  ylab("log weekly incidence Brasil") + 
  theme_minimal()

#Ajuste un modelo log-lineal a los datos de incidencia semanal

f_bra <- incidence::fit(i_weekly_bra)
f_bra



plot(i_weekly_bra, fit = f_bra)

#-------Encontrando una fecha límite adecuada para el modelo log-lineal de Brasil, en función de los retrasos observados

summary(as.numeric(icases_clean%>% filter (iso3 %in% "BRA")$date_diagnosis - icases_clean%>% filter (iso3 %in% "BRA")$date_onset))
#no me corre
#VERIFICAR


#¿cuántas semanas debe descartar al final de la epicurva? 

# Semanas a descartar al final de la epicurva
n_weeks_to_discard_bra <- 2

min_date_bra<- min(i_daily_bra$dates)

max_date_bra <- max(i_daily_bra$dates) - n_weeks_to_discard_bra * 7


# Para truncar la incidencia semanal 
i_weekly_trunc_bra <- subset(i_weekly_bra, 
                             from = min_date_bra, 
                             to = max_date_bra) # descarte las últimas semanas de datos


#Vuelva a montar y a graficar el modelo logarítmico lineal, pero utilizando los datos truncados i_weekly_trunc. 
#Los resultados deben ser como los siquientes:

f_bra<- incidence::fit(i_weekly_trunc_bra)
f_bra



plot(i_weekly_trunc_bra, fit = f_bra)

#Observe las estadísticas resumidas de su ajuste:

summary(f_bra$model)


#Puede observar la bondad del ajuste (Rsquared), la pendiente estimada (tasa de crecimiento/growth rate) y el tiempo de duplicación correspondiente 
#como se muestra a continuación:

#¿El modelo se ajusta bien a los datos?
adjRsq_model_fit_bra <- summary(f_bra$model)$adj.r.squared


#¿Cuál es la tasa de crecimiento estimada de la epidemia?
daily_growth_rate_bra <- f_bra$model$coefficients['dates.x']


# intervalo de confianza:
daily_growth_rate_CI_bra <- confint(f_bra$model, 'dates.x', level=0.95)


#¿Cuál es el tiempo de duplicación de la epidemia?
doubling_time_days_bra <- log(2) / daily_growth_rate_bra

# intervalo de confianza:
doubling_time_days_CI_bra <- log(2) / rev(daily_growth_rate_CI_bra)





#--------Curvas de incidencia diaria de Argentina desde enero 2022

i_daily_arg <- incidence((icases_clean %>% filter (iso3 %in% "ARG"))$date_onset)

i_daily_arg




plot(i_daily_arg, border = "black")


#--------Curvas de incidencia semanal de Argentina desde enero 2022

i_weekly_arg<- incidence((icases_clean %>% filter (iso3 %in% "ARG"))$date_onset,interval=7)

i_weekly_arg




plot(i_weekly_arg,border = "black")

#-------Estimación de la tasa de crecimiento de Argentina mediante un modelo log-lineal 

#Grafico de la incidencia transformada logaritmicamente

ggplot(as.data.frame(i_weekly_arg)) + 
  geom_point(aes(x = dates, y = log(counts))) + 
  scale_x_incidence(i_weekly_arg) +
  xlab("date") +
  ylab("log weekly incidence Argentina") + 
  theme_minimal()

#Ajuste un modelo log-lineal a los datos de incidencia semanal

f_arg <- incidence::fit(i_weekly_arg)
f_arg




plot(i_weekly_arg, fit = f_arg)

#-------Encontrando una fecha límite adecuada para el modelo log-lineal de Argentina, en función de los retrasos observados

summary(as.numeric(icases_clean%>% filter (iso3 %in% "ARG")$date_diagnosis - icases_clean%>% filter (iso3 %in% "ARG")$date_onset))
#no me corre
#VERIFICAR


#¿cuántas semanas debe descartar al final de la epicurva? 

# Semanas a descartar al final de la epicurva
n_weeks_to_discard_arg <- 2

min_date_arg<- min(i_daily_arg$dates)

max_date_arg<- max(i_daily_arg$dates) - n_weeks_to_discard_arg * 7


# Para truncar la incidencia semanal 
i_weekly_trunc_arg <- subset(i_weekly_arg, 
                             from = min_date_arg, 
                             to = max_date_arg) # descarte las últimas semanas de datos


#Vuelva a montar y a graficar el modelo logarítmico lineal, pero utilizando los datos truncados i_weekly_trunc. 
#Los resultados deben ser como los siquientes:

f_arg<- incidence::fit(i_weekly_trunc_arg)
f_arg



plot(i_weekly_trunc_arg, fit = f_arg)

#Observe las estadísticas resumidas de su ajuste:

summary(f_arg$model)




#Puede observar la bondad del ajuste (Rsquared), la pendiente estimada (tasa de crecimiento/growth rate) y el tiempo de duplicación correspondiente 
#como se muestra a continuación:

#¿El modelo se ajusta bien a los datos?
adjRsq_model_fit_arg <- summary(f_arg$model)$adj.r.squared


#¿Cuál es la tasa de crecimiento estimada de la epidemia?
daily_growth_rate_arg <- f_arg$model$coefficients['dates.x']


# intervalo de confianza:
daily_growth_rate_CI_arg <- confint(f_arg$model, 'dates.x', level=0.95)

#¿Cuál es el tiempo de duplicación de la epidemia?
doubling_time_days_arg <- log(2) / daily_growth_rate_arg

# intervalo de confianza:
doubling_time_days_CI_arg <- log(2) / rev(daily_growth_rate_CI_arg)





#--------Curvas de incidencia diaria de Peru desde enero 2022

i_daily_per <- incidence((icases_clean %>% filter (iso3 %in% "PER"))$date_onset)

i_daily_per



plot(i_daily_per, border = "black")


#--------Curvas de incidencia semanal de Peru desde enero 2022

i_weekly_per<- incidence((icases_clean %>% filter (iso3 %in% "PER"))$date_onset,interval=7)

i_weekly_per



plot(i_weekly_per,border = "black")

#-------Estimación de la tasa de crecimiento de Peru mediante un modelo log-lineal 

#Grafico de la incidencia transformada logaritmicamente

ggplot(as.data.frame(i_weekly_per)) + 
  geom_point(aes(x = dates, y = log(counts))) + 
  scale_x_incidence(i_weekly_per) +
  xlab("date") +
  ylab("log weekly incidence Peru") + 
  theme_minimal()

#Ajuste un modelo log-lineal a los datos de incidencia semanal

f_per <- incidence::fit(i_weekly_per)
f_per


plot(i_weekly_per, fit = f_per)

#-------Encontrando una fecha límite adecuada para el modelo log-lineal de Peru, en función de los retrasos observados

summary(as.numeric(icases_clean%>% filter (iso3 %in% "PER")$date_diagnosis - icases_clean%>% filter (iso3 %in% "PER")$date_onset))
#no me corre
#VERIFICAR


#¿cuántas semanas debe descartar al final de la epicurva? 

# Semanas a descartar al final de la epicurva
n_weeks_to_discard_per <- 2

min_date_per<- min(i_daily_per$dates)

max_date_per<- max(i_daily_per$dates) - n_weeks_to_discard_per * 7


# Para truncar la incidencia semanal 
i_weekly_trunc_per <- subset(i_weekly_per, 
                             from = min_date_per, 
                             to = max_date_per) # descarte las últimas semanas de datos


#Vuelva a montar y a graficar el modelo logarítmico lineal, pero utilizando los datos truncados i_weekly_trunc. 
#Los resultados deben ser como los siquientes:

f_per<- incidence::fit(i_weekly_trunc_per)
f_per



plot(i_weekly_trunc_per, fit = f_per)

#Observe las estadísticas resumidas de su ajuste:

summary(f_per$model)


#Puede observar la bondad del ajuste (Rsquared), la pendiente estimada (tasa de crecimiento/growth rate) y el tiempo de duplicación correspondiente 
#como se muestra a continuación:

#¿El modelo se ajusta bien a los datos?
adjRsq_model_fit_per <- summary(f_per$model)$adj.r.squared


#¿Cuál es la tasa de crecimiento estimada de la epidemia?
daily_growth_rate_per <- f_per$model$coefficients['dates.x']


# intervalo de confianza:
daily_growth_rate_CI_per <- confint(f_per$model, 'dates.x', level=0.95)


#¿Cuál es el tiempo de duplicación de la epidemia?
doubling_time_days_per <- log(2) / daily_growth_rate_per


# intervalo de confianza:
doubling_time_days_CI_per <- log(2) / rev(daily_growth_rate_CI_per)





#--------Curvas de incidencia diaria de Mexico desde enero 2022

i_daily_mex <- incidence((icases_clean %>% filter (iso3 %in% "MEX"))$date_onset)

i_daily_mex



plot(i_daily_mex, border = "black")


#--------Curvas de incidencia semanal de Mexico desde enero 2022

i_weekly_mex<- incidence((icases_clean %>% filter (iso3 %in% "MEX"))$date_onset,interval=7)

i_weekly_mex



plot(i_weekly_mex,border = "black")

#-------Estimación de la tasa de crecimiento de Mexico mediante un modelo log-lineal 

#Grafico de la incidencia transformada logaritmicamente

ggplot(as.data.frame(i_weekly_mex)) + 
  geom_point(aes(x = dates, y = log(counts))) + 
  scale_x_incidence(i_weekly_mex) +
  xlab("date") +
  ylab("log weekly incidence Mexico") + 
  theme_minimal()

#Ajuste un modelo log-lineal a los datos de incidencia semanal

f_mex <- incidence::fit(i_weekly_mex)
f_mex



plot(i_weekly_mex, fit = f_mex)

#-------Encontrando una fecha límite adecuada para el modelo log-lineal de Mexico, en función de los retrasos observados

summary(as.numeric(icases_clean%>% filter (iso3 %in% "MEX")$date_diagnosis - icases_clean%>% filter (iso3 %in% "MEX")$date_onset))
#no me corre
#VERIFICAR


#¿cuántas semanas debe descartar al final de la epicurva? 

# Semanas a descartar al final de la epicurva
n_weeks_to_discard_mex <- 2

min_date_mex<- min(i_daily_mex$dates)

max_date_mex<- max(i_daily_mex$dates) - n_weeks_to_discard_mex* 7


# Para truncar la incidencia semanal 
i_weekly_trunc_mex <- subset(i_weekly_mex, 
                             from = min_date_mex, 
                             to = max_date_mex) # descarte las últimas semanas de datos


#Vuelva a montar y a graficar el modelo logarítmico lineal, pero utilizando los datos truncados i_weekly_trunc. 
#Los resultados deben ser como los siquientes:

f_mex<- incidence::fit(i_weekly_trunc_mex)
f_mex

plot(i_weekly_trunc_mex, fit = f_mex)

#Observe las estadísticas resumidas de su ajuste:

summary(f_mex$model)


#Puede observar la bondad del ajuste (Rsquared), la pendiente estimada (tasa de crecimiento/growth rate) y el tiempo de duplicación correspondiente 
#como se muestra a continuación:

#¿El modelo se ajusta bien a los datos?
adjRsq_model_fit_mex <- summary(f_mex$model)$adj.r.squared


#¿Cuál es la tasa de crecimiento estimada de la epidemia?
daily_growth_rate_mex <- f_mex$model$coefficients['dates.x']


# intervalo de confianza:
daily_growth_rate_CI_mex <- confint(f_mex$model, 'dates.x', level=0.95)


#¿Cuál es el tiempo de duplicación de la epidemia?
doubling_time_days_mex <- log(2) / daily_growth_rate_mex


# intervalo de confianza:
doubling_time_days_CI_mex <- log(2) / rev(daily_growth_rate_CI_mex)



#Guardar elementos importantes...
info <- list(i_daily_col=i_daily_col,
             min_date_col=min_date_col,
             max_date_col=max_date_col,
             i_daily_bra=i_daily_bra,
             min_date_bra=min_date_bra,
             max_date_bra=max_date_bra,
             i_daily_per=i_daily_per,
             min_date_per=min_date_per,
             max_date_per=max_date_per,
             i_daily_arg=i_daily_arg,
             min_date_arg=min_date_arg,
             max_date_arg=max_date_arg,
             i_daily_mex=i_daily_mex,
             min_date_mex=min_date_mex,
             max_date_mex=max_date_mex)

saveRDS(info, "data/cases/info.RDS")





