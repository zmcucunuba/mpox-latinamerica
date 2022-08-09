

# Monkeypox model
rm(list =ls())

# Libraries
library(deSolve)



p     <- 0.8   # Probabilidad p de ser infeccioso cuando el encuentro se hace con alguien que está contagiado
b_gg  <- 1     # Num contactos general-general
b_gh  <- 0.05  # Num contactos general-general
b_gh  <- 1e-6  # Num contactos general-pr
b_lg  <- 0.05  # Num contactos HSH-general
b_ll  <- 0.25  # Num contactos 
b_lh  <- 1e-2  # Num contactos 
b_hg  <- 1e-6  # Num contactos 
b_hl  <- 1e-2  # Num contactos 
b_hh  <- 1     # 
N_g   <- 1.e6  # Num pob general
N_l   <- 2.e4  # Num pob HSH low risk
N_h   <- 1.5e2 # Num pob HSH high risk

latency  <- 7.6 # Charniga el al. 2022 https://doi.org/10.1101/2022.06.22.22276713
e_to_i_rate <- user(1) # 1/latency
recovery_rate <- user(1) # 1/infectious period

alpha       <- 0.01  # From Exposed to Infectious (¿tal vez muy pequeño?, 1/7.6 = 0.13 1/latency)
inf_period  <- 21 # Days
gamma       <- 1/inf_period  # From Infectious to Recovered

est_S       <- 0.9 # Ref: XXXX (pendiente estimar nosotros dado que COlombia suspendió vacunación en 198xx)




# ---------------------------------------
# ------------    Model
# ---------------------------------------


mpxmodel <- function(t, state, parameters) {
  
  #------ Compartimentos del modelo
  Sg <- x[1]    #  Susceptible 
  Eg <- x[2]    #  Exposed general
  Ig <- x[3]    #  Infecious general
  Rg <- x[4]    #  Recovered general
  
  #### Terminar compartimentos
  
  with(as.list(c(state, parameters)),{
    #------ Force-of-Infection (Incidence per susceptible population)
    lambda_g <- p*(beta_gg*I_g + beta_gh*I_h + beta_gp*I_p)*S_g
    lambda_h <- p*(beta_hg*I_g + beta_hh*I_h + beta_hp*I_p)*S_h
    lambda_p <- p*(beta_pg*I_g + beta_ph*I_h + beta_pp*I_p)*S_p
    
    #------ Rate of Change
    dS_g <- -lambda_g
    dE_g <- lambda_g - alpha * E_g
    dI_g <- alpha * E_g - gamma * I_g
    dR_g <- gamma * I_g
    
    dS_h <- -lambda_h
    dE_h <- lambda_h - alpha * E_h
    dI_h <- alpha * E_h - gamma * I_h
    dR_h <- gamma * I_h
    
    dS_p <- -lambda_p
    dE_p <- lambda_p - alpha * E_p
    dI_p <- alpha * E_p - gamma * I_p
    dR_p <- gamma * I_p
    # return the rate of change
    return(list(c(dS_g ,dE_g, dI_g, dR_g, dS_h, dE_h, dI_h, dR_h, dS_p, dE_p, dI_p, dR_p)))
  }) 
}



# Time 
tdays  <- seq(1, 365 * 2 , by = 1)

# ---------------------------------------
# ----------- Specifying parameters
# ---------------------------------------

parameters <- c(beta_gg = b_gg/N_g, 
                beta_gh = b_gh/N_g, 
                beta_gp = b_gp/N_g,
                beta_hg = b_hg/N_h, 
                beta_hh = b_hh/N_h, 
                beta_hp = b_hp/N_h,
                beta_pg = b_pg/N_p, 
                beta_ph = b_ph/N_p, 
                beta_pp = b_pp/N_p,
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



# Solving the equations
out <- as.data.frame(ode(y      = xstart,   # meaning??  
                         times  = tdays,   # meaning??     
                         fun    = mpxmodel,   # meaning??  
                         parms  = parameters)
)  # meaning??   



plot(out, xlab = "tiempo", ylab = "Población(t)")


