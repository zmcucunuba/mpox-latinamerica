######## Monkeypox model#######

rm(list =ls())

# Libraries
library(deSolve)
library(ggplot2)
library(cowplot)



p     <- 1   # Probabilidad p de volverse E dado que tiene un contacto con I. 

# b_gg  <- 1     # Num contactos general-general
# b_gl  <- 0.05  # Num contactos general-HSH low risk
# b_gh  <- 1e-6  # Num contactos general-HSH high risk
# b_lg  <- 0.05  # Num contactos HSH low risk-general
# b_ll  <- 0.25  # Num contactos HSH low risk-HSH low risk
# b_lh  <- 1e-2  # Num contactos HSH low risk-HSH high risk
# b_hg  <- 1e-6  # Num contactos HSH high risk-general risk
# b_hl  <- 1e-2  # Num contactos HSH high risk-HSH low risk
# b_hh  <- 1     # Num contactos HSH high risk-HSH high risk

# -------------------------------------------------------
# ----------- Number of sex partners during the last year 
# ------------------------------------------------------

# Experimento David Santiago

# Num parejas sexuales en últimas 4 semanas
npg <- 1.5/365
npl <- 5/365
nph <- 50/365

b_gg  <- npg * 0.9   # Num contactos general-general
b_gl  <- npg * 0.09  # Num contactos general-HSH low risk
b_gh  <- npg * 0.01  # Num contactos general-HSH high risk

b_ll  <- npl * 0.9   # Num contactos HSH low risk-HSH low risk
b_lg  <- npl * 0.05  # Num contactos HSH low risk-HSH high risk
b_lh  <- npl * 0.05  # Num contactos HSH high risk-HSH high risk

b_hh  <- nph * 0.9   # Num contactos HSH high risk-HSH low risk
b_hl  <- nph * 0.09  # Num contactos HSH high risk-general risk
b_hg  <- nph * 0.01  # Num contactos HSH low risk-general

# b_gg  <- 1.5/365 # Num contactos general-general
# b_gl  <- 0.5/365 # Num contactos general-HSH low risk
# b_gh  <- 0.2/365 # Num contactos general-HSH high risk
# b_ll  <- 3/365 # Num contactos HSH low risk-HSH low risk
# b_lh  <- 2/365 # Num contactos HSH low risk-HSH high risk
# b_hh  <- 50/365 # Num contactos HSH high risk-HSH high risk

# b_lg  <- b_gl  # Num contactos HSH low risk-general
# b_hg  <- b_gh  # Num contactos HSH high risk-general risk
# b_hl  <- b_lh  # Num contactos HSH high risk-HSH low risk


# ----- Tamaño de las poblaciones

p_HSH       <- 0.024 #(https://www.sdp.gov.co/sites/default/files/boletin_15.pdf)
p_HSH_l     <- 0.8 #supuesto.. buscar.. 
N_total     <- 5.425e6  # Num pob general 18-65a de Bogota 
                  #(SALUDATA-https://saludata.saludcapital.gov.co/osb/index.php/datos-de-salud/demografia/piramidepoblacional/)
N_HSH       <- N_total/2 * p_HSH  # N_total/2 porque corresponde sólo a hombres que son aprox el 50% de N_g
N_g         <- N_total - N_HSH
N_l         <- N_HSH * p_HSH_l  # Num pob HSH low risk
N_h         <- N_HSH * (1- p_HSH_l)  # Num pob HSH high risk

latency     <- 7.6 #Periodo de incubacion, Charniga el al. 2022 https://doi.org/10.1101/2022.06.22.22276713
inf_period  <- 21 # Dias # CDC https://www.cdc.gov/poxvirus/monkeypox/clinicians/monitoring.html
est_S       <- 0.9 #Susceptibilidad estimada de la infeccion en Colombia
#Taube et al. 2022 https://doi.org/10.1101/2022.07.29.22278217
#pendiente estimar nosotros dado que Colombia suspendió vacunación en 1984)


alpha <- 1/latency
gamma <- 1/inf_period

# ---------------------------------------
# ----------- Specifying parameters
# ---------------------------------------

parameters <- c(beta_gg = b_gg/N_g,#Tasa =  num promedio de contactos gg / Total poblacion general por unidad de tiempo (dia)
                beta_gl = b_gl/N_g,#Tasa =  num promedio de contactos gl / Total poblacion general por unidad de tiempo (dia)
                beta_gh = b_gh/N_g,#Tasa =  num promedio de contactos gh / Total poblacion general por unidad de tiempo (dia)
                beta_lg = b_lg/N_l,#Tasa =  num promedio de contactos lg / Total poblacion LowRisk por unidad de tiempo (dia)
                beta_ll = b_ll/N_l,#Tasa =  num promedio de contactos ll / Total poblacion LowRisk por unidad de tiempo (dia)
                beta_lh = b_lh/N_l,#Tasa =  num promedio de contactos lh / Total poblacion LowRisk por unidad de tiempo (dia)
                beta_hg = b_hg/N_h,#Tasa = num promedio de contactos hg / Total de poblacion HighRisk por unidad de tiempo (dia)
                beta_hl = b_hl/N_h,#Tasa = num promedio de contactos hl / Total de poblacion HighRisk por unidad de tiempo (dia)
                beta_hh = b_hh/N_h,#Tasa = num promedio de contactos hh / Total de poblacion HighRisk por unidad de tiempo (dia)
                alpha = alpha, 
                gamma = gamma)



# ---------------------------------------
# ----------- # Initial conditions (State)
# ---------------------------------------

xstart <- c(S_g = N_g * est_S, 
            E_g = 0, 
            I_g = 0,
            R_g = N_g * (1-est_S),
            
            S_l = N_l * est_S, 
            E_l = 0, 
            I_l = 0, 
            R_l = N_l * (1-est_S),
            
            S_h = N_h * est_S, 
            E_h = 0, 
            I_h = 1, 
            R_h = N_h * (1-est_S)
)


# --------------------------------------
#------------   Time 
#---------------------------------------

tdays  <- seq(1, 365 * 5 , by = 1)



# ---------------------------------------
# ------------    Model
# ---------------------------------------

#------ Compartimentos del modelo
# S_g <- x[1]    #  Susceptible general
# E_g <- x[2]    #  Expuesto general
# I_g <- x[3]    #  Infeccioso general
# R_g <- x[4]    #  Recuperado general
# S_l <- x[5]    #  Susceptible Low risk
# E_l <- x[6]    #  Expuesto Low risk
# I_l <- x[7]    #  Infeccioso low risk
# R_l <- x[8]    #  Recuperado low risk
# S_h <- x[9]    #  Susceptible high risk
# E_h <- x[10]   #  Expuesto high risk
# I_h <- x[11]   #  Infeccioso high risk
# R_h <- x[12]   #  Recuperado high risk


mpxmodel <- function(tdays, xstart, parameters) {
  
  with(as.list(c(xstart, parameters)),{
    
    #------ Force-of-Infection (Incidence per susceptible population)
    lambda_g <- p*(beta_gg*I_g + beta_gl*I_l + beta_gh*I_h)*S_g
    lambda_l <- p*(beta_lg*I_g + beta_ll*I_l + beta_lh*I_h)*S_l
    lambda_h <- p*(beta_hg*I_g + beta_hl*I_l + beta_hh*I_h)*S_h
    
    #------ Rate of Change
    dS_g <- -lambda_g
    dE_g <- lambda_g - alpha * E_g
    dI_g <- alpha * E_g - gamma * I_g
    dR_g <- gamma * I_g
    
    dS_l <- -lambda_l
    dE_l <- lambda_l - alpha * E_l
    dI_l <- alpha * E_l - gamma * I_l
    dR_l <- gamma * I_l
    
    dS_h <- -lambda_h
    dE_h <- lambda_h - alpha * E_h
    dI_h <- alpha * E_h - gamma * I_h
    dR_h <- gamma * I_h
    # return the rate of change
    return(list(c(dS_g ,dE_g, dI_g, dR_g, dS_l, dE_l, dI_l, dR_l, dS_h, dE_h, dI_h, dR_h)))
  }) 
}


# --------------------------------------
#------------   Solving the equation
#---------------------------------------

out <- ode(y = xstart,times=tdays,fun=mpxmodel, parms=parameters)   
out.df <-as.data.frame(out)
# plot(out, xlab = "tiempo", ylab = "Población(t)")


# --------------------------------------
#------------   Graphics
#---------------------------------------

mytheme4 <- theme(text=element_text(colour="black")) +
  theme(panel.grid = element_line(colour = "white")) +
  theme(panel.background = element_rect(fill = "#B2B2B2")) +
  theme_bw() +
  theme(legend.title=element_text(size=12,face="bold"),
        legend.background = element_rect(fill='white', size=0.5,linetype="solid"),
        legend.text=element_text(size=10),
        legend.key=element_rect(colour="white",
                                fill='white',
                                size= 0.25,
                                linetype="solid"))
theme_set(mytheme4)

x_limit_days <- 500

p1 <- 
  ggplot(out.df,aes(x=time)) +
  ggtitle(bquote("HSH high risk")) +
  geom_line(aes(y=S_h,colour="Susceptible")) +
  geom_line(aes(y=E_h, colour= "Exposed"))  +
  geom_line(aes(y=I_h,colour="Infectious")) +
  geom_line(aes(y=R_h,colour="Recovered")) +
  # coord_cartesian(xlim = c(0, x_limit_days)) +
  ylab(label="Number") + xlab(label="Time (days)") +
  scale_colour_manual("",
                      breaks=c("Susceptible", "Exposed", "Infectious","Recovered"),
                      values=c("blue","orange", "red","darkgreen"))


p2 <-
  ggplot(out.df,aes(x=time))+
  ggtitle(bquote("HSH Low risk")) +
  geom_line(aes(y=S_l,colour="Susceptible"))+
  geom_line(aes(y=E_l, colour= "Exposed")) +
  geom_line(aes(y=I_l,colour="Infectious"))+
  geom_line(aes(y=R_l,colour="Recovered"))+
  ylab(label="Number") + xlab(label="Time (days)") +
  # coord_cartesian(xlim = c(0, x_limit_days)) +
  scale_colour_manual("",
                      breaks=c("Susceptible", "Exposed", "Infectious","Recovered"),
                      values=c("blue","orange", "red","darkgreen"))




p3 <-
  ggplot(out.df,aes(x=time))+
  ggtitle(bquote("General pop")) +
  geom_line(aes(y=S_g,colour="Susceptible"))+
  geom_line(aes(y=E_g, colour= "Exposed")) +
  geom_line(aes(y=I_g,colour="Infectious"))+
  geom_line(aes(y=R_g,colour="Recovered"))+
  ylab(label="Number") + xlab(label="Time (days)") +
  # coord_cartesian(xlim = c(0, x_limit_days)) +
  scale_colour_manual("",
                      breaks=c("Susceptible", "Exposed", "Infectious","Recovered"),
                      values=c("blue","orange", "red","darkgreen"))




pp <- cowplot::plot_grid(p1, p2, p3, nrow = 3, align = "hv", labels = "AUTO")
pp
# save_plot("figs/pp.png", pp, base_width = 7, base_height = 7)






