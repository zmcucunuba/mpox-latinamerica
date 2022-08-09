######## Monkeypox model#######

rm(list =ls())

# Libraries
library(deSolve)


p     <- 0.8   # Probabilidad p de ser infeccioso cuando el encuentro se hace con alguien que está contagiado
b_gg  <- 1     # Num contactos general-general
b_gl  <- 0.05  # Num contactos general-HSH low risk
b_gh  <- 1e-6  # Num contactos general-HSH high risk
b_lg  <- 0.05  # Num contactos HSH low risk-general
b_ll  <- 0.25  # Num contactos HSH low risk-HSH low risk
b_lh  <- 1e-2  # Num contactos HSH low risk-HSH high risk
b_hg  <- 1e-6  # Num contactos HSH high risk-general risk
b_hl  <- 1e-2  # Num contactos HSH high risk-HSH low risk
b_hh  <- 1     # Num contactos HSH high risk-HSH high risk
N_g   <- 1.e6  # Num pob general
N_l   <- 2.e4  # Num pob HSH low risk
N_h   <- 1.5e2 # Num pob HSH high risk

latency  <- 7.6 #Periodo de incubacion, Charniga el al. 2022 https://doi.org/10.1101/2022.06.22.22276713
e_to_i_rate <- 1/latency # Tasa de expuestos a infecciosos 
inf_period  <- 21 # Dias # CDC https://www.cdc.gov/poxvirus/monkeypox/clinicians/monitoring.html
recovery_rate <- 1/inf_period  # Tasa de recueracion (de infeccioso a recuperado) 
est_S       <- 0.9 #Susceptibilidad estimada de la infeccion en Colombia
#Taube et al. 2022 https://doi.org/10.1101/2022.07.29.22278217
#pendiente estimar nosotros dado que Colombia suspendió vacunación en 1984)

alpha <-  0.01 
gamma <- 1/21

# ---------------------------------------
# ----------- Specifying parameters
# ---------------------------------------

parameters <- c(beta_gg = b_gg/N_g, 
                beta_gl = b_gl/N_g, 
                beta_gh = b_gh/N_g,
                beta_lg = b_lg/N_l, 
                beta_ll = b_ll/N_l, 
                beta_lh = b_lh/N_l,
                beta_hg = b_hg/N_h, 
                beta_hl = b_hl/N_h, 
                beta_hh = b_hh/N_h,
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

tdays  <- seq(1, 365 * 2 , by = 1)



# ---------------------------------------
# ------------    Model
# ---------------------------------------

  #------ Compartimentos del modelo

  Sg <- x[1]    #  Susceptible general
  Eg <- x[2]    #  Expuesto general
  Ig <- x[3]    #  Infeccioso general
  Rg <- x[4]    #  Recuperado general
  Sl <- x[5]    #  Susceptible Low risk
  El <- x[6]    #  Expuesto Low risk
  Il <- x[7]    #  Infeccioso low risk
  Rl <- x[8]    #  Recuperado low risk
  Sh <- x[9]    #  Susceptible high risk
  Eh <- x[10]   #  Expuesto high risk
  Ih <- x[11]   #  Infeccioso high risk
  Rh <- x[12]   #  Recuperado high risk

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
out.df<-as.data.frame(out)
plot(out, xlab = "tiempo", ylab = "Población(t)")


# --------------------------------------
#------------   Graphics
#---------------------------------------

require(ggplot2)

mytheme4 <- theme(text=element_text(colour="black")) +
  theme(panel.grid = element_line(colour = "white")) +
  theme(panel.background = element_rect(fill = "#B2B2B2")) +
  theme_bw() 
theme_set(mytheme4)

title <- bquote("Monkeypox basic SEIR model")
subtit <- bquote("Group HSH high risk")

res<-ggplot(out.df,aes(x=time))+
  ggtitle(bquote(atop(bold(.(title)))))+
  geom_line(aes(y=S_h,colour="Susceptible"))+
  geom_line(aes(y=E_h, colour= "Exposed")) +
  geom_line(aes(y=I_h,colour="Infectious"))+
  geom_line(aes(y=R_h,colour="Recovered"))+
  ylab(label="Number")+
  xlab(label="Time (days)")+
  theme(legend.justification=c(1,0), legend.position=c(1,0.5))+
  theme(legend.title=element_text(size=12,face="bold"),
        legend.background = element_rect(fill='#FFFFFF',
                                         size=0.5,linetype="solid"),
        legend.text=element_text(size=10),
        legend.key=element_rect(colour="#FFFFFF",
                                fill='#C2C2C2',
                                size=0.25,
                                linetype="solid"))+
  scale_colour_manual("Compartments",
                      breaks=c("Susceptible", "Exposed", "Infectious","Recovered"),
                      values=c("blue","yellow", "red","darkgreen"))

print(res)

subtit <- bquote("Group HSH Low risk")

res<-ggplot(out.df,aes(x=time))+
  ggtitle(bquote(atop(bold(.(title)))))+
  geom_line(aes(y=S_l,colour="Susceptible"))+
  geom_line(aes(y=E_l, colour= "Exposed")) +
  geom_line(aes(y=I_l,colour="Infectious"))+
  geom_line(aes(y=R_l,colour="Recovered"))+
  ylab(label="Number")+
  xlab(label="Time (days)")+
  theme(legend.justification=c(1,0), legend.position=c(1,0.5))+
  theme(legend.title=element_text(size=12,face="bold"),
        legend.background = element_rect(fill='#FFFFFF',
                                         size=0.5,linetype="solid"),
        legend.text=element_text(size=10),
        legend.key=element_rect(colour="#FFFFFF",
                                fill='#C2C2C2',
                                size=0.25,
                                linetype="solid"))+
  scale_colour_manual("Compartments",
                      breaks=c("Susceptible", "Exposed", "Infectious","Recovered"),
                      values=c("blue","yellow", "red","darkgreen"))

print(res)

subtit <- bquote("General population")

res<-ggplot(out.df,aes(x=time))+
  ggtitle(bquote(atop(bold(.(title)))))+
  geom_line(aes(y=S_g,colour="Susceptible"))+
  geom_line(aes(y=E_g, colour= "Exposed")) +
  geom_line(aes(y=I_g,colour="Infectious"))+
  geom_line(aes(y=R_g,colour="Recovered"))+
  ylab(label="Number")+
  xlab(label="Time (days)")+
  theme(legend.justification=c(1,0), legend.position=c(1,0.5))+
  theme(legend.title=element_text(size=12,face="bold"),
        legend.background = element_rect(fill='#FFFFFF',
                                         size=0.5,linetype="solid"),
        legend.text=element_text(size=10),
        legend.key=element_rect(colour="#FFFFFF",
                                fill='#C2C2C2',
                                size=0.25,
                                linetype="solid"))+
  scale_colour_manual("Compartments",
                      breaks=c("Susceptible", "Exposed", "Infectious","Recovered"),
                      values=c("blue","yellow", "red","darkgreen"))

print(res)








