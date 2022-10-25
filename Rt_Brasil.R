#loads EpiEstim library
library(EpiEstim)

#imports RJ state file
MPV.RJ <- read.csv("~/RStudio/MPV RJ.csv", sep=";")

#create the dataframes with the Rt estimates for RJ state
Rt_RJ <- estimate_R(MPV.RJ,
                    method = "parametric_si",
                    config = make_config(list(
                      mean_si = 9.8,
                      std_si = 9.1
                    )))

#plots the epidemic curve graph for RJ state, from the 20th day onward
plot(Rt_RJ,
     what = "incid",
     options_R = list(col = palette(),
                      transp = 0.2,
                      xlim = c(20,45),
                      ylim = c(0,4),
                      xlab = "Time",
                      ylab = "Incidence"),
     legend = TRUE)

#plots the Rt estimate graph for RJ state, from the 20th day onward
plot(Rt_RJ,
     what = "R",
     options_R = list(col = palette(),
                      transp = 0.2,
                      xlim = c(20,45),
                      ylim = c(0,4),
                      xlab = "Time",
                      ylab = "Rt"),
     legend = TRUE)

#generates a CSV with RJ Rt estimates data
write.csv(Rt_RJ$R, "/Users/isaacschrarstzhaupt/Desktop/R.RJ.csv")

#imports SP state file
MPV.SP <- read.csv("~/RStudio/MPV SP.csv", sep=";")

#create the dataframes with the Rt estimates for SP state
Rt_SP <- estimate_R(MPV.SP,
                    method = "parametric_si",
                    config = make_config(list(
                      mean_si = 9.8,
                      std_si = 9.1
                    )))
#plots the epidemic curve graph for SP state, from the 20th day onward
plot(Rt_SP,
     what = "incid",
     options_R = list(col = palette(),
                      transp = 0.2,
                      xlim = c(20,45),
                      ylim = c(0,4),
                      xlab = "Time",
                      ylab = "Incidence"),
     legend = TRUE)

#plots the Rt estimate graph for SP state, from the 20th day onward
plot(Rt_SP,
     what = "R",
     options_R = list(col = palette(),
                      transp = 0.2,
                      xlim = c(20,45),
                      ylim = c(0,4),
                      xlab = "Time",
                      ylab = "Rt"),
     legend = TRUE)

#generates a CSV with SP Rt estimates data
write.csv(Rt_SP$R, "/Users/isaacschrarstzhaupt/Desktop/R.SP.csv")

#imports MG state file
MPV.MG <- read.csv("~/RStudio/MPV MG.csv", sep=";")

#create the dataframes with the Rt estimates for MG state
Rt_MG <- estimate_R(MPV.MG,
                    method = "parametric_si",
                    config = make_config(list(
                      mean_si = 9.8,
                      std_si = 9.1
                    )))

#plots the epidemic curve graph for MG state, from the 20th day onward
plot(Rt_MG,
     what = "incid",
     options_R = list(col = palette(),
                      transp = 0.2,
                      xlim = c(20,45),
                      ylim = c(0,4),
                      xlab = "Time",
                      ylab = "Incidence"),
     legend = TRUE)

#plots the Rt estimate graph for MG state, from the 20th day onward
plot(Rt_MG,
     what = "R",
     options_R = list(col = palette(),
                      transp = 0.2,
                      xlim = c(20,45),
                      ylim = c(0,4),
                      xlab = "Time",
                      ylab = "Rt"),
     legend = TRUE)

#generates a CSV with MG Rt estimates data
write.csv(Rt_MG$R, "/Users/isaacschrarstzhaupt/Desktop/R.MG.csv")

#imports GO+FD state file
MPV.DF.GO <- read.csv("~/RStudio/MPV DF-GO.csv", sep=";")

#create the dataframes with the Rt estimates for GO+FD state
Rt_DFGO <- estimate_R(MPV.DF.GO,
                      method = "parametric_si",
                      config = make_config(list(
                        mean_si = 9.8,
                        std_si = 9.1
                      )))

#plots the epidemic curve graph for GO+FD state, from the 20th day onward

plot(Rt_DFGO,
     what = "incid",
     options_R = list(col = palette(),
                      transp = 0.2,
                      xlim = c(20,45),
                      ylim = c(0,4),
                      xlab = "Time",
                      ylab = "Incidence"),
     legend = TRUE)

#plots the Rt estimate graph for GO+FD state, from the 20th day onward
plot(Rt_DFGO,
     what = "R",
     options_R = list(col = palette(),
                      transp = 0.2,
                      xlim = c(20,45),
                      ylim = c(0,4),
                      xlab = "Time",
                      ylab = "Rt"),
     legend = TRUE)

#generates a CSV with GO+FD Rt estimates data
write.csv(Rt_DFGO$R, "/Users/isaacschrarstzhaupt/Desktop/R.DF.GO.csv")