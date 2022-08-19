
# Script para entender los datos de UK
rm(list=ls())

library(tidyverse)
library(haven)


nat2000 <- read_sav("data/contacts/natsal_2000_for_archive.sav")
nat2010 <- read_sav("data/contacts/eul_level_natsal_2010_web_data_for_archive.sav")
wave1 <- read_sav("data/contacts/natsal_covid_wave_1_archive.sav")
nat2010arch <- read_sav("data/contacts/eul_natsal_2010_for_archive.sav")

str(wave1)

# Num ocasiones 4 semanas sexo opuesto
wave1$num_opp_4k <- as.numeric(wave1$OppSex4wks_w1)
wave1$num_opp_4k[wave1$num_opp_4k == -9] <- NA   # Prefer not to say
wave1$num_opp_4k[wave1$num_opp_4k == 9999] <- NA # Prefer not to responde

wave1 %>% filter (num_opp_4k > 100) %>% summarise(n= n())

summary(wave1$num_opp_4k)
xx <- summary(wave1$num_opp_4k)
hist (wave1$num_opp_4k[wave1$num_opp_4k < 50], breaks = 30, main = "Num Encuentros Opp Sex. 4 sem")
abline(v = xx[4], col = "red")


# Num nuevas parejas sexuales sexo opuesto 4 semanas 

wave1$num_opp_new_4k <- as.numeric(wave1$D_Oppnonew4wksMerge_w1)
wave1$num_opp_new_4k[wave1$num_opp_new_4k == -9] <- NA   # Prefer not to say
wave1$num_opp_new_4k[wave1$num_opp_new_4k == 9999] <- NA # Prefer not to responde
wave1$num_opp_new_4k[wave1$num_opp_new_4k == 9995] <- 1 # at least 1 partner
wave1$num_opp_new_4k[wave1$num_opp_new_4k == 9996] <- 2 # at least 2 partners

wave1 %>% filter (num_opp_new_4k > 100) %>% summarise(n= n())

summary(wave1$num_opp_new_4k)
yy <- summary(wave1$num_opp_new_4k)
hist (wave1$num_opp_new_4k[wave1$num_opp_new_4k < 100], breaks = 300, 
      main = "Num Nuevas Parejas Sex Opp 4 sem")
abline(v = yy[4], col = "red")



wave1$num_opp_1y <- as.numeric(wave1$Opp1yr_w1)
wave1$num_opp_1y[wave1$num_opp_1y == -9] <- NA   # Prefer not to say
wave1$num_opp_1y[wave1$num_opp_1y == 9999] <- NA # Prefer not to responde
wave1$num_opp_1y[wave1$num_opp_1y == 9995] <- 1 # at least 1 partner
wave1$num_opp_1y[wave1$num_opp_1y == 9996] <- 2 # at least 2 partners


summary(wave1$num_opp_1y)
yy <- summary(wave1$num_opp_1y)
hist (wave1$num_opp_1y[wave1$num_opp_1y < 100], breaks = 300, 
      main = "Num Nuevas Parejas Sex Opp 4 sem")
abline(v = yy[4], col = "red")

