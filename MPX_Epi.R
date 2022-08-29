
#--------------------------------------------- 
#
# Script Analisis Descriptivo Tiempo Real
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

#------ Limpiar los datos

icases[-c(6198, 11892, 1218, 6098), ]  #se remueven observaciones con fechas antes de 2022

icases_clean <- icases[-c(6198, 11892, 1218, 6098), ]  

names(acases) <- epitrix::clean_labels(names(acases))
names(icases_clean) <- epitrix::clean_labels(names(icases_clean))


mistakes <- which(icases_clean$date_onset >= icases_clean$date_diagnosis) 
icases_clean[mistakes, ] #se remueven observaciones con fechas de inicio de sintomas >= a fecha de diagnostico

icases_clean <- icases_clean[-mistakes, ] 

acases <- acases %>% filter(classification == "Confirmed")

selected_countries <- c("COL", "BRA", "PER", "BOL", "MEX", "ARG")


#------ Grafica de casos acumulados de paises seleccionados

ggplot(data = acases %>% filter (iso3 %in% selected_countries)) +
  geom_line(aes(x= date, y = cases)) +
  facet_wrap(~ iso3, scales = "free_y") +
  ggtitle(label = "Cummulative cases")

#------ Grafica de casos nuevos de paises seleccionados

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


#------ Grafica de casos de Bolivia discriminado por inicio de sintomas y fecha de reporte (ERROR, NO CORRE, VERIFICAR)


bol_not_cases <- icases_clean %>% 
  filter (iso3 %in% "BOL") %>%
  group_by(report_date) %>%
  summarise (cases = n())


bol_fis_cases <- icases_clean %>% 
  filter (iso3 %in% "BOL") %>%
  group_by(date_onset) %>%
  summarise (cases = n())

ggplot() +
  geom_col(data = bol_not_cases, aes(x= report_date, y = cases, fill = "date_report"), 
           alpha = .5, colour = "black") +
  geom_col(data = bol_fis_cases, aes(x= date_onset, y = cases, fill = "date_onset"),
           alpha = .5, colour = "black") +
  ggtitle(label = "BOL cases") +
  guides(fill = guide_legend(title="")) 


#--------Curvas de incidencia diaria de paises seleccionados desde enero 2022

i_daily_selected_countries <- incidence((icases_clean %>% filter (iso3 %in% selected_countries))$date_onset)
                     
i_daily_selected_countries

# <incidence object>
#   [3365 cases from days 2022-02-27 to 2022-08-11]
# 
# $counts: matrix with 166 rows and 1 columns
# $n: 3365 cases in total
# $dates: 166 dates marking the left-side of bins
# $interval: 1 day
# $timespan: 166 days
# $cumulative: FALSE


plot(i_daily_selected_countries, border = "black")


#--------Curvas de incidencia semanal de paises seleccionados desde enero 2022

i_weekly_selected_countries <- incidence((icases_clean %>% filter (iso3 %in% selected_countries))$date_onset,interval=7)
                    
i_weekly_selected_countries

# 
# > i_weekly
# <incidence object>
#   [3365 cases from days 2022-02-21 to 2022-08-08]
# [3365 cases from ISO weeks 2022-W08 to 2022-W32]
# 
# $counts: matrix with 25 rows and 1 columns
# $n: 3365 cases in total
# $dates: 25 dates marking the left-side of bins
# $interval: 7 days
# $timespan: 169 days
# $cumulative: FALSE

plot(i_weekly_selected_countries,border = "black")

#-------Estimación de la tasa de crecimientode paises seleccionados mediante un modelo log-lineal 

#Grafico de la incidencia transformada logaritmicamente

ggplot(as.data.frame(i_weekly_selected_countries)) + 
  geom_point(aes(x = dates, y = log(counts))) + 
  scale_x_incidence(i_weekly_selected_countries) +
  xlab("date") +
  ylab("log weekly incidence selected countries") + 
  theme_minimal()

#Ajuste un modelo log-lineal a los datos de incidencia semanal

f_selected_countries <- incidence::fit(i_weekly_selected_countries)
f_selected_countries

# <incidence_fit object>
#   
#   $model: regression of log-incidence over time
# 
# $info: list containing the following items:
#   $r (daily growth rate):
#   [1] 0.04720902
# 
# $r.conf (confidence interval):
#   2.5 %     97.5 %
#   [1,] 0.02868264 0.06573541
# 
# $doubling (doubling time in days):
#   [1] 14.68252
# 
# $doubling.conf (confidence interval):
#   2.5 %   97.5 %
#   [1,] 10.5445 24.16609
# 
# $pred: data.frame of incidence predictions (16 rows, 5 columns)
# 
plot(i_weekly_selected_countries, fit = f_selected_countries)

#-------Encontrando una fecha límite adecuada para el modelo log-lineal de paises seleccionados, en función de los retrasos observados

#Utilizando la gráfica del logaritmo (incidencia) que graficó anteriormente, y pensando en por qué el crecimiento exponencial no puede observarse 
#en las últimas semanas, elija una fecha límite y ajuste el modelo logarítmico lineal a una sección adecuada de la epicurva donde crea que puede 
#estimar de manera más confiable la tasa de crecimiento r, y el tiempo de duplicación.

#Es posible que desee examinar cuánto tiempo después de la aparición de los síntomas hasta el diagnostico; para obtener un reporte de una fecha
#especifica, siga estos comandos. 

summary(as.numeric(icases_clean$date_diagnosis -  icases_clean$date_onset)) #no me corre con paises seleccionados sino con el total de las americas

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   1.000   5.000   7.000   8.083  10.000 158.000   11798 
#   
  
#¿cuántas semanas debe descartar al final de la epicurva? 

# Semanas a descartar al final de la epicurva
n_weeks_to_discard_selected_countries <- 1

min_date_selected_countries <- min(i_daily_selected_countries$dates)
#2022-02-27
max_date_selected_countries <- max(i_daily_selected_countries$dates) - n_weeks_to_discard_selected_countries * 7
#2022-08-04

# Para truncar la incidencia semanal 
i_weekly_trunc_selected_countries <- subset(i_weekly_selected_countries, 
                         from = min_date_selected_countries, 
                         to = max_date_selected_countries) # descarte las últimas semanas de datos

# incidencia diaria truncada.No la usamos para la regresión lineal pero se puede usar más adelante

#Vuelva a montar y a graficar el modelo logarítmico lineal, pero utilizando los datos truncados i_weekly_trunc. 
#Los resultados deben ser como los siquientes:

f_selected_countries <- incidence::fit(i_weekly_trunc_selected_countries)
f_selected_countries
# 
# <incidence_fit object>
#   
#   $model: regression of log-incidence over time
# 
# $info: list containing the following items:
#   $r (daily growth rate):
#   [1] 0.06973557
# 
# $r.conf (confidence interval):
#   2.5 %     97.5 %
#   [1,] 0.0510655 0.08840564
# 
# $doubling (doubling time in days):
#   [1] 9.939651
# 
# $doubling.conf (confidence interval):
#   2.5 %   97.5 %
#   [1,] 7.840532 13.57369
# 
# $pred: data.frame of incidence predictions (14 rows, 5 columns)

plot(i_weekly_trunc_selected_countries, fit = f_selected_countries)

#Observe las estadísticas resumidas de su ajuste:

# summary(f_selected_countries$model)
# 
# Call:
#   stats::lm(formula = log(counts) ~ dates.x, data = df)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.6784 -0.9454  0.1188  0.7971  1.7386 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -1.982706   0.754223  -2.629    0.022 *  
#   dates.x      0.069736   0.008569   8.138 3.15e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.066 on 12 degrees of freedom
# Multiple R-squared:  0.8466,	Adjusted R-squared:  0.8338 
# F-statistic: 66.23 on 1 and 12 DF,  p-value: 3.154e-06

#Puede observar la bondad del ajuste (Rsquared), la pendiente estimada (tasa de crecimiento/growth rate) y el tiempo de duplicación correspondiente 
#como se muestra a continuación:

#¿El modelo se ajusta bien a los datos?
adjRsq_model_fit_selected_countries <- summary(f_selected_countries$model)$adj.r.squared
#0.83

#¿Cuál es la tasa de crecimiento estimada de la epidemia?
daily_growth_rate_selected_countries<- f_selected_countries$model$coefficients['dates.x']
#0.069

# intervalo de confianza:
daily_growth_rate_CI_selected_countries <- confint(f_selected_countries$model, 'dates.x', level=0.95)
#0.0511-0.0884

#¿Cuál es el tiempo de duplicación de la epidemia?
doubling_time_days_selected_countries <- log(2) / daily_growth_rate_selected_countries
#9.94
# intervalo de confianza:
doubling_time_days_CI_selected_countries <- log(2) / rev(daily_growth_rate_CI_selected_countries)
#7.84-13.57




#--------Curvas de incidencia diaria de Colombia desde enero 2022

i_daily_col <- incidence((icases_clean %>% filter (iso3 %in% "COL"))$date_onset)

i_daily_col

# <incidence object>
#   [111 cases from days 2022-05-29 to 2022-08-11]
# 
# $counts: matrix with 75 rows and 1 columns
# $n: 111 cases in total
# $dates: 75 dates marking the left-side of bins
# $interval: 1 day
# $timespan: 75 days
# $cumulative: FALSE


plot(i_daily_col, border = "black")


#--------Curvas de incidencia semanal de Colombia desde enero 2022

i_weekly_col <- incidence((icases_clean %>% filter (iso3 %in% "COL"))$date_onset,interval=7)

i_weekly_col


# <incidence object>
#   [111 cases from days 2022-05-23 to 2022-08-08]
# [111 cases from ISO weeks 2022-W21 to 2022-W32]
# 
# $counts: matrix with 12 rows and 1 columns
# $n: 111 cases in total
# $dates: 12 dates marking the left-side of bins
# $interval: 7 days
# $timespan: 78 days
# $cumulative: FALSE


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
# 
# <incidence_fit object>
#   
#   $model: regression of log-incidence over time
# 
# $info: list containing the following items:
#   $r (daily growth rate):
#   [1] 0.04295371
# 
# $r.conf (confidence interval):
#   2.5 %     97.5 %
#   [1,] 0.01727046 0.06863696
# 
# $doubling (doubling time in days):
#   [1] 16.13707
# 
# $doubling.conf (confidence interval):
#   2.5 %   97.5 %
#   [1,] 10.09874 40.13485
# 
# $pred: data.frame of incidence predictions (11 rows, 5 columns)

plot(i_weekly_col, fit = f_col)

#-------Encontrando una fecha límite adecuada para el modelo log-lineal de Colombia, en función de los retrasos observados

summary(as.numeric(icases_clean%>% filter (iso3 %in% "COL")$date_diagnosis - icases_clean%>% filter (iso3 %in% "COL")$date_onset))
#no me corre
#VERIFICAR


#¿cuántas semanas debe descartar al final de la epicurva? 

# Semanas a descartar al final de la epicurva
n_weeks_to_discard_col <- 1

min_date_col<- min(i_daily_col$dates)
#2022-05-29
max_date_col <- max(i_daily_col$dates) - n_weeks_to_discard_col * 7
#2022-08-04

# Para truncar la incidencia semanal 
i_weekly_trunc_col <- subset(i_weekly_col, 
                         from = min_date_col, 
                         to = max_date_col) # descarte las últimas semanas de datos


#Vuelva a montar y a graficar el modelo logarítmico lineal, pero utilizando los datos truncados i_weekly_trunc. 
#Los resultados deben ser como los siquientes:

f_col<- incidence::fit(i_weekly_trunc_col)
f_col

# <incidence_fit object>
#   
#   $model: regression of log-incidence over time
# 
# $info: list containing the following items:
#   $r (daily growth rate):
#   [1] 0.06147808
# 
# $r.conf (confidence interval):
#   2.5 %     97.5 %
#   [1,] 0.03031934 0.09263682
# 
# $doubling (doubling time in days):
#   [1] 11.2747
# 
# $doubling.conf (confidence interval):
#   2.5 %   97.5 %
#   [1,] 7.482415 22.86155
# 
# $pred: data.frame of incidence predictions (9 rows, 5 columns)

plot(i_weekly_trunc_col, fit = f_col)

#Observe las estadísticas resumidas de su ajuste:

summary(f_col$model)

# Call:
#   stats::lm(formula = log(counts) ~ dates.x, data = df)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.9394 -0.4842 -0.2214  0.5353  1.2372 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept) -0.75925    0.55875  -1.359   0.2164   
# dates.x      0.06148    0.01318   4.666   0.0023 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.7656 on 7 degrees of freedom
# Multiple R-squared:  0.7567,	Adjusted R-squared:  0.7219 
# F-statistic: 21.77 on 1 and 7 DF,  p-value: 0.002299

#Puede observar la bondad del ajuste (Rsquared), la pendiente estimada (tasa de crecimiento/growth rate) y el tiempo de duplicación correspondiente 
#como se muestra a continuación:

#¿El modelo se ajusta bien a los datos?
adjRsq_model_fit_col <- summary(f_col$model)$adj.r.squared
#0.7219

#¿Cuál es la tasa de crecimiento estimada de la epidemia?
daily_growth_rate_col <- f_col$model$coefficients['dates.x']
#0.0615

# intervalo de confianza:
daily_growth_rate_CI_col <- confint(f_col$model, 'dates.x', level=0.95)
#0.0303-0.0926

#¿Cuál es el tiempo de duplicación de la epidemia?
doubling_time_days_col <- log(2) / daily_growth_rate_col
#11.3
# intervalo de confianza:
doubling_time_days_CI_col <- log(2) / rev(daily_growth_rate_CI_col)
#7.48-22.86




#--------Curvas de incidencia diaria de Brasil desde enero 2022

i_daily_bra <- incidence((icases_clean %>% filter (iso3 %in% "BRA"))$date_onset)

i_daily_bra

# <incidence object>
#   [2165 cases from days 2022-02-27 to 2022-08-04]
# 
# $counts: matrix with 159 rows and 1 columns
# $n: 2165 cases in total
# $dates: 159 dates marking the left-side of bins
# $interval: 1 day
# $timespan: 159 days
# $cumulative: FALSE


plot(i_daily_bra, border = "black")


#--------Curvas de incidencia semanal de Brasil desde enero 2022

i_weekly_bra<- incidence((icases_clean %>% filter (iso3 %in% "BRA"))$date_onset,interval=7)

i_weekly_bra


# <incidence object>
#   [2165 cases from days 2022-02-21 to 2022-08-01]
# [2165 cases from ISO weeks 2022-W08 to 2022-W31]
# 
# $counts: matrix with 24 rows and 1 columns
# $n: 2165 cases in total
# $dates: 24 dates marking the left-side of bins
# $interval: 7 days
# $timespan: 162 days
# $cumulative: FALSE

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

# <incidence_fit object>
#   
#   $model: regression of log-incidence over time
# 
# $info: list containing the following items:
#   $r (daily growth rate):
#   [1] 0.06507167
# 
# $r.conf (confidence interval):
#   2.5 %     97.5 %
#   [1,] 0.04362775 0.08651558
# 
# $doubling (doubling time in days):
#   [1] 10.65206
# 
# $doubling.conf (confidence interval):
#   2.5 %   97.5 %
#   [1,] 8.011819 15.88776
# 
# $pred: data.frame of incidence predictions (11 rows, 5 columns)


plot(i_weekly_bra, fit = f_bra)

#-------Encontrando una fecha límite adecuada para el modelo log-lineal de Brasil, en función de los retrasos observados

summary(as.numeric(icases_clean%>% filter (iso3 %in% "BRA")$date_diagnosis - icases_clean%>% filter (iso3 %in% "BRA")$date_onset))
#no me corre
#VERIFICAR


#¿cuántas semanas debe descartar al final de la epicurva? 

# Semanas a descartar al final de la epicurva
n_weeks_to_discard_bra <- 1

min_date_bra<- min(i_daily_bra$dates)
#2022-02-27
max_date_bra <- max(i_daily_bra$dates) - n_weeks_to_discard_bra * 7
#2022-07-28

# Para truncar la incidencia semanal 
i_weekly_trunc_bra <- subset(i_weekly_bra, 
                             from = min_date_bra, 
                             to = max_date_bra) # descarte las últimas semanas de datos


#Vuelva a montar y a graficar el modelo logarítmico lineal, pero utilizando los datos truncados i_weekly_trunc. 
#Los resultados deben ser como los siquientes:

f_bra<- incidence::fit(i_weekly_trunc_bra)
f_bra

# <incidence_fit object>
#   
#   $model: regression of log-incidence over time
# 
# $info: list containing the following items:
#   $r (daily growth rate):
#   [1] 0.06507167
# 
# $r.conf (confidence interval):
#   2.5 %     97.5 %
#   [1,] 0.04362775 0.08651558
# 
# $doubling (doubling time in days):
#   [1] 10.65206
# 
# $doubling.conf (confidence interval):
#   2.5 %   97.5 %
#   [1,] 8.011819 15.88776

plot(i_weekly_trunc_bra, fit = f_bra)

#Observe las estadísticas resumidas de su ajuste:

summary(f_bra$model)

# Call:
#   stats::lm(formula = log(counts) ~ dates.x, data = df)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.5781 -0.7398  0.4066  0.6839  1.3728 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -1.600547   0.842377  -1.900   0.0899 .  
# dates.x      0.065072   0.009479   6.865 7.35e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.9943 on 9 degrees of freedom
# Multiple R-squared:  0.8396,	Adjusted R-squared:  0.8218 
# F-statistic: 47.12 on 1 and 9 DF,  p-value: 7.353e-05


#Puede observar la bondad del ajuste (Rsquared), la pendiente estimada (tasa de crecimiento/growth rate) y el tiempo de duplicación correspondiente 
#como se muestra a continuación:

#¿El modelo se ajusta bien a los datos?
adjRsq_model_fit_bra <- summary(f_bra$model)$adj.r.squared
#0.82

#¿Cuál es la tasa de crecimiento estimada de la epidemia?
daily_growth_rate_bra <- f_bra$model$coefficients['dates.x']
#0.0651

# intervalo de confianza:
daily_growth_rate_CI_bra <- confint(f_bra$model, 'dates.x', level=0.95)
#0.0436-0.0865

#¿Cuál es el tiempo de duplicación de la epidemia?
doubling_time_days_bra <- log(2) / daily_growth_rate_bra
#10.7
# intervalo de confianza:
doubling_time_days_CI_bra <- log(2) / rev(daily_growth_rate_CI_bra)
#8.01-15.89




#--------Curvas de incidencia diaria de Argentina desde enero 2022

i_daily_arg <- incidence((icases_clean %>% filter (iso3 %in% "ARG"))$date_onset)

i_daily_arg

# <incidence object>
#   [70 cases from days 2022-05-15 to 2022-08-09]
# 
# $counts: matrix with 87 rows and 1 columns
# $n: 70 cases in total
# $dates: 87 dates marking the left-side of bins
# $interval: 1 day
# $timespan: 87 days
# $cumulative: FALSE



plot(i_daily_arg, border = "black")


#--------Curvas de incidencia semanal de Argentina desde enero 2022

i_weekly_arg<- incidence((icases_clean %>% filter (iso3 %in% "ARG"))$date_onset,interval=7)

i_weekly_arg

# <incidence object>
#   [70 cases from days 2022-05-09 to 2022-08-08]
# [70 cases from ISO weeks 2022-W19 to 2022-W32]
# 
# $counts: matrix with 14 rows and 1 columns
# $n: 70 cases in total
# $dates: 14 dates marking the left-side of bins
# $interval: 7 days
# $timespan: 92 days
# $cumulative: FALSE


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

# <incidence_fit object>
#   
#   $model: regression of log-incidence over time
# 
# $info: list containing the following items:
#   $r (daily growth rate):
#   [1] 0.03091585
# 
# $r.conf (confidence interval):
#   2.5 %     97.5 %
#   [1,] 0.01669146 0.04514025
# 
# $doubling (doubling time in days):
#   [1] 22.42045
# 
# $doubling.conf (confidence interval):
#   2.5 %   97.5 %
#   [1,] 15.35541 41.52706
# 
# $pred: data.frame of incidence predictions (11 rows, 5 columns)


plot(i_weekly_arg, fit = f_arg)

#-------Encontrando una fecha límite adecuada para el modelo log-lineal de Argentina, en función de los retrasos observados

summary(as.numeric(icases_clean%>% filter (iso3 %in% "ARG")$date_diagnosis - icases_clean%>% filter (iso3 %in% "ARG")$date_onset))
#no me corre
#VERIFICAR


#¿cuántas semanas debe descartar al final de la epicurva? 

# Semanas a descartar al final de la epicurva
n_weeks_to_discard_arg <- 1

min_date_arg<- min(i_daily_arg$dates)
#2022-05-15
max_date_arg<- max(i_daily_arg$dates) - n_weeks_to_discard_arg * 7
#2022-08-02

# Para truncar la incidencia semanal 
i_weekly_trunc_arg <- subset(i_weekly_arg, 
                             from = min_date_arg, 
                             to = max_date_arg) # descarte las últimas semanas de datos


#Vuelva a montar y a graficar el modelo logarítmico lineal, pero utilizando los datos truncados i_weekly_trunc. 
#Los resultados deben ser como los siquientes:

f_arg<- incidence::fit(i_weekly_trunc_arg)
f_arg

# <incidence_fit object>
#   
#   $model: regression of log-incidence over time
# 
# $info: list containing the following items:
#   $r (daily growth rate):
#   [1] 0.04731912
# 
# $r.conf (confidence interval):
#   2.5 %     97.5 %
#   [1,] 0.04023337 0.05440486
# 
# $doubling (doubling time in days):
#   [1] 14.64835
# 
# $doubling.conf (confidence interval):
#   2.5 %   97.5 %
#   [1,] 12.74054 17.22816
# 
# $pred: data.frame of incidence predictions (9 rows, 5 columns)

plot(i_weekly_trunc_arg, fit = f_arg)

#Observe las estadísticas resumidas de su ajuste:

summary(f_arg$model)

# Call:
#   stats::lm(formula = log(counts) ~ dates.x, data = df)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.275856 -0.022033  0.008647  0.055378  0.309201 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.190315   0.127064  -1.498    0.178    
# dates.x      0.047319   0.002997  15.791  9.9e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1741 on 7 degrees of freedom
# Multiple R-squared:  0.9727,	Adjusted R-squared:  0.9688 
# F-statistic: 249.4 on 1 and 7 DF,  p-value: 9.896e-07


#Puede observar la bondad del ajuste (Rsquared), la pendiente estimada (tasa de crecimiento/growth rate) y el tiempo de duplicación correspondiente 
#como se muestra a continuación:

#¿El modelo se ajusta bien a los datos?
adjRsq_model_fit_arg <- summary(f_arg$model)$adj.r.squared
#0.968

#¿Cuál es la tasa de crecimiento estimada de la epidemia?
daily_growth_rate_arg <- f_arg$model$coefficients['dates.x']
#0.0473

# intervalo de confianza:
daily_growth_rate_CI_arg <- confint(f_arg$model, 'dates.x', level=0.95)
#0.0402-0.0544

#¿Cuál es el tiempo de duplicación de la epidemia?
doubling_time_days_arg <- log(2) / daily_growth_rate_arg
#14.6
# intervalo de confianza:
doubling_time_days_CI_arg <- log(2) / rev(daily_growth_rate_CI_arg)
#12.7-17.2




#--------Curvas de incidencia diaria de Peru desde enero 2022

i_daily_per <- incidence((icases_clean %>% filter (iso3 %in% "PER"))$date_onset)

i_daily_per

# <incidence object>
#   [773 cases from days 2022-06-15 to 2022-08-11]
# 
# $counts: matrix with 58 rows and 1 columns
# $n: 773 cases in total
# $dates: 58 dates marking the left-side of bins
# $interval: 1 day
# $timespan: 58 days
# $cumulative: FALSE


plot(i_daily_per, border = "black")


#--------Curvas de incidencia semanal de Peru desde enero 2022

i_weekly_per<- incidence((icases_clean %>% filter (iso3 %in% "PER"))$date_onset,interval=7)

i_weekly_per

# <incidence object>
#   [773 cases from days 2022-06-13 to 2022-08-08]
# [773 cases from ISO weeks 2022-W24 to 2022-W32]
# 
# $counts: matrix with 9 rows and 1 columns
# $n: 773 cases in total
# $dates: 9 dates marking the left-side of bins
# $interval: 7 days
# $timespan: 57 days
# $cumulative: FALSE


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

# <incidence_fit object>
#   
#   $model: regression of log-incidence over time
# 
# $info: list containing the following items:
#   $r (daily growth rate):
#   [1] 0.03780395
# 
# $r.conf (confidence interval):
#   2.5 %     97.5 %
#   [1,] -0.01774982 0.09335771
# 
# $doubling (doubling time in days):
#   [1] 18.33531
# 
# $doubling.conf (confidence interval):
#   2.5 %    97.5 %
#   [1,] 7.424638 -39.05095
# 
# $pred: data.frame of incidence predictions (9 rows, 5 columns)


plot(i_weekly_per, fit = f_per)

#-------Encontrando una fecha límite adecuada para el modelo log-lineal de Peru, en función de los retrasos observados

summary(as.numeric(icases_clean%>% filter (iso3 %in% "PER")$date_diagnosis - icases_clean%>% filter (iso3 %in% "PER")$date_onset))
#no me corre
#VERIFICAR


#¿cuántas semanas debe descartar al final de la epicurva? 

# Semanas a descartar al final de la epicurva
n_weeks_to_discard_per <- 1

min_date_per<- min(i_daily_per$dates)
#2022-06-15
max_date_per<- max(i_daily_per$dates) - n_weeks_to_discard_per * 7
#2022-08-04

# Para truncar la incidencia semanal 
i_weekly_trunc_per <- subset(i_weekly_per, 
                             from = min_date_per, 
                             to = max_date_per) # descarte las últimas semanas de datos


#Vuelva a montar y a graficar el modelo logarítmico lineal, pero utilizando los datos truncados i_weekly_trunc. 
#Los resultados deben ser como los siquientes:

f_per<- incidence::fit(i_weekly_trunc_per)
f_per

# <incidence_fit object>
#   
#   $model: regression of log-incidence over time
# 
# $info: list containing the following items:
#   $r (daily growth rate):
#   [1] 0.06491749
# 
# $r.conf (confidence interval):
#   2.5 %     97.5 %
#   [1,] 0.03289493 0.09694005
# 
# $doubling (doubling time in days):
#   [1] 10.67736
# 
# $doubling.conf (confidence interval):
#   2.5 %   97.5 %
#   [1,] 7.150266 21.07155
# 
# $pred: data.frame of incidence predictions (7 rows, 5 columns)

plot(i_weekly_trunc_per, fit = f_per)

#Observe las estadísticas resumidas de su ajuste:

summary(f_per$model)

# Call:
#   stats::lm(formula = log(counts) ~ dates.x, data = df)
# 
# Residuals:
#   1        2        3        4        5        6        7 
# -0.56678  0.01489  0.41694  0.39056  0.14489  0.19694 -0.59746 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  2.73746    0.35152   7.788 0.000559 ***
#   dates.x      0.06492    0.01246   5.211 0.003435 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.4614 on 5 degrees of freedom
# Multiple R-squared:  0.8445,	Adjusted R-squared:  0.8134 
# F-statistic: 27.16 on 1 and 5 DF,  p-value: 0.003435

#Puede observar la bondad del ajuste (Rsquared), la pendiente estimada (tasa de crecimiento/growth rate) y el tiempo de duplicación correspondiente 
#como se muestra a continuación:

#¿El modelo se ajusta bien a los datos?
adjRsq_model_fit_per <- summary(f_per$model)$adj.r.squared
#0.81

#¿Cuál es la tasa de crecimiento estimada de la epidemia?
daily_growth_rate_per <- f_per$model$coefficients['dates.x']
#0.0649

# intervalo de confianza:
daily_growth_rate_CI_per <- confint(f_per$model, 'dates.x', level=0.95)
#0.0329-0.0969

#¿Cuál es el tiempo de duplicación de la epidemia?
doubling_time_days_per <- log(2) / daily_growth_rate_per
#10.7
# intervalo de confianza:
doubling_time_days_CI_per <- log(2) / rev(daily_growth_rate_CI_per)
#7.15-21.07




#--------Curvas de incidencia diaria de Mexico desde enero 2022

i_daily_mex <- incidence((icases_clean %>% filter (iso3 %in% "MEX"))$date_onset)

i_daily_mex

# <incidence object>
#   [246 cases from days 2022-05-19 to 2022-08-11]
# 
# $counts: matrix with 85 rows and 1 columns
# $n: 246 cases in total
# $dates: 85 dates marking the left-side of bins
# $interval: 1 day
# $timespan: 85 days
# $cumulative: FALSE

plot(i_daily_mex, border = "black")


#--------Curvas de incidencia semanal de Mexico desde enero 2022

i_weekly_mex<- incidence((icases_clean %>% filter (iso3 %in% "MEX"))$date_onset,interval=7)

i_weekly_mex

# <incidence object>
#   [246 cases from days 2022-05-16 to 2022-08-08]
# [246 cases from ISO weeks 2022-W20 to 2022-W32]
# 
# $counts: matrix with 13 rows and 1 columns
# $n: 246 cases in total
# $dates: 13 dates marking the left-side of bins
# $interval: 7 days
# $timespan: 85 days
# $cumulative: FALSE

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

# <incidence_fit object>
#   
#   $model: regression of log-incidence over time
# 
# $info: list containing the following items:
#   $r (daily growth rate):
#   [1] 0.03600945
# 
# $r.conf (confidence interval):
#   2.5 %     97.5 %
#   [1,] 0.01804768 0.05397123
# 
# $doubling (doubling time in days):
#   [1] 19.24903
# 
# $doubling.conf (confidence interval):
#   2.5 %   97.5 %
#   [1,] 12.8429 38.40644
# 
# $pred: data.frame of incidence predictions (12 rows, 5 columns)

plot(i_weekly_mex, fit = f_mex)

#-------Encontrando una fecha límite adecuada para el modelo log-lineal de Mexico, en función de los retrasos observados

summary(as.numeric(icases_clean%>% filter (iso3 %in% "MEX")$date_diagnosis - icases_clean%>% filter (iso3 %in% "MEX")$date_onset))
#no me corre
#VERIFICAR


#¿cuántas semanas debe descartar al final de la epicurva? 

# Semanas a descartar al final de la epicurva
n_weeks_to_discard_mex <- 1

min_date_mex<- min(i_daily_mex$dates)
#2022-06-19
max_date_mex<- max(i_daily_mex$dates) - n_weeks_to_discard_mex* 7
#2022-08-04

# Para truncar la incidencia semanal 
i_weekly_trunc_mex <- subset(i_weekly_mex, 
                             from = min_date_mex, 
                             to = max_date_mex) # descarte las últimas semanas de datos


#Vuelva a montar y a graficar el modelo logarítmico lineal, pero utilizando los datos truncados i_weekly_trunc. 
#Los resultados deben ser como los siquientes:

f_mex<- incidence::fit(i_weekly_trunc_mex)
f_mex
# 
# <incidence_fit object>
#   
#   $model: regression of log-incidence over time
# 
# $info: list containing the following items:
#   $r (daily growth rate):
#   [1] 0.04132656
# 
# $r.conf (confidence interval):
#   2.5 %     97.5 %
#   [1,] 0.03035736 0.05229577
# 
# $doubling (doubling time in days):
#   [1] 16.77244
# 
# $doubling.conf (confidence interval):
#   2.5 %   97.5 %
#   [1,] 13.25436 22.83292
# 
# $pred: data.frame of incidence predictions (10 rows, 5 columns)

plot(i_weekly_trunc_mex, fit = f_mex)

#Observe las estadísticas resumidas de su ajuste:

summary(f_mex$model)

# Call:
#   stats::lm(formula = log(counts) ~ dates.x, data = df)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.45061 -0.15874 -0.00753  0.23601  0.41614 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 1.346044   0.192003   7.011 0.000111 ***
#   dates.x     0.041327   0.004757   8.688  2.4e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3024 on 8 degrees of freedom
# Multiple R-squared:  0.9042,	Adjusted R-squared:  0.8922 
# F-statistic: 75.48 on 1 and 8 DF,  p-value: 2.4e-05

#Puede observar la bondad del ajuste (Rsquared), la pendiente estimada (tasa de crecimiento/growth rate) y el tiempo de duplicación correspondiente 
#como se muestra a continuación:

#¿El modelo se ajusta bien a los datos?
adjRsq_model_fit_mex <- summary(f_mex$model)$adj.r.squared
#0.89

#¿Cuál es la tasa de crecimiento estimada de la epidemia?
daily_growth_rate_mex <- f_mex$model$coefficients['dates.x']
#0.0413

# intervalo de confianza:
daily_growth_rate_CI_mex <- confint(f_mex$model, 'dates.x', level=0.95)
#0.0304-0.0523

#¿Cuál es el tiempo de duplicación de la epidemia?
doubling_time_days_mex <- log(2) / daily_growth_rate_mex
#16.8

# intervalo de confianza:
doubling_time_days_CI_mex <- log(2) / rev(daily_growth_rate_CI_mex)
#13.3-22.8















