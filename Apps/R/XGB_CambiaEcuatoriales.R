rm(list = ls())
library(xgboost)
library(data.table)

# Cambia las coordenadas ecuatoriales
# Se giran sistematicamente las observaciones dentro de un  rango de ascensiones rectas y declinaciones
# para un rango de variaciones de r (distancia)


CaE <- function(xyz) {
    
  # Convierte Cartesianas a Ecuatoriales

  x <- xyz[, 1]
  y <- xyz[, 2]
  z <- xyz[, 3]
  distancia   <- sqrt(x * x + y * y + z * z)
  declinacion <- asin(z / distancia);
  ascension   <- declinacion
  for (i in 1:length(x)){
    if (x[i] == 0)
    {
      if (y[i] > 0)
      {
        ascension[i] <- pi / 2.0;
      } else {
        ascension[i] <- pi + pi / 2.0;
      }
    } else if (x[i] < 0) {
      v            <- y[i] / (distancia[i] * cos(declinacion[i]));
      ascension[i] <- pi - asin(v);
    } else {
      if (y[i] > 0)
      {
        v            <- y[i] / (distancia[i] * cos(declinacion[i]));
        ascension[i] <- asin(v);
      } else {
        v            <- y[i] / (distancia[i] * cos(declinacion[i]));
        ascension[i] <- 2.0 * pi + asin(v);
      }
    }
  }
  return (cbind(ascension, declinacion, distancia))
}
EaC <- function(adr){
    
  # Convierte Ecuatoriales a Cartesianas

  a <- adr[, 1]
  d <- adr[, 2]
  r <- adr[, 3]
  l <- r * cos(d);
  x <- l * cos(a);
  y <- l * sin(a);
  z <- r * sin(d);
  return (cbind(x, y, z))
}

desamb <- "_Lucas"

# horas

desde_ar <- +1.94
hasta_ar <- +1.95
inc_ar   <- 1/60/60

# grados

desde_de <- -4.5
hasta_de <- -3.5
inc_de   <- 1/60

# celdas

desde_r  <- +2
hasta_r  <- +5
inc_r    <- 0.1

# escala

u_min  <- 0.9
u_max  <- 1.1
u_inc  <- 0.01

ponderar <- TRUE
nr  <- 620
eta <- 0.084
md  <- 13
cs  <- 0.94
ss  <- 1
cb  <- 1

disco   <- "C"
marcas  <- "444"
dmina   <- "600000"
dmaxa   <- "1000000"
anguloa <- "_8_12_15_70_"
casoa   <- paste(sep="",disco,":/Gly/G_1_10_",marcas,anguloa,dmina,"_",dmaxa)
dminb   <- "600000"
dmaxb   <- "1000000"
angulob <- "_12_16_15_70_"
casob   <- paste(sep="",disco,":/Gly/G_1_10_",marcas,angulob,dminb,"_",dmaxb)
porcentaje_entrenar <- ""

datosentrenar <- read.table(paste(sep="",casoa,"_muestra_",porcentaje_entrenar,"tot_meA.csv"), header=TRUE, sep=";", na.strings="NA", dec=",", strip.white=TRUE)
ponderae      <- ifelse(datosentrenar$n == 0, 1, datosentrenar$n)
sumapon       <- sum(ponderae)
entrena       <- datosentrenar[,-c(5)]
Eentrena      <- as.numeric(entrena$E)
setDT(entrena)
coordenadasEntrena  <- model.matrix(~.+0, data = entrena[,-c("E"), with=F]) 
dentrena <- xgb.DMatrix(data = coordenadasEntrena, label = Eentrena) 

iter_ar <- abs(round((hasta_ar - desde_ar)/inc_ar, 0)) + 1
iter_de <- abs(round((hasta_de - desde_de)/inc_de, 0)) + 1
iter_r  <- abs(round((hasta_r - desde_r)/inc_r, 0)) + 1
iter_u  <- abs(round((u_max - u_min)/u_inc, 0)) + 1
pruebas <- iter_ar*iter_de*iter_r*iter_u

aciertos_mejorRU <- array(0,dim=c(iter_de))
mejoresRU       <- array(0,dim=c(iter_de))
fiacto_mejorRU  <- paste(sep="",casoa,"_muestra_",porcentaje_entrenar,"_cambiar_Ecuatoriales_resultados_Ar",desamb,".csv")
fimejorRU       <- paste(sep="",casoa,"_muestra_",porcentaje_entrenar,"_cambiar_Ecuatoriales_resultados_Mr",desamb,".csv")
fi              <- paste(sep="",casoa,"_muestra_",porcentaje_entrenar,"_cambiar_Ecuatoriales_log",desamb,".txt")

cat(sprintf("\nEntrenar: %s\nProbar  : %s\n",paste(sep="",casoa,"_muestra_",porcentaje_entrenar,"tot_meA.csv"),paste(sep="",casob,"_muestra_",porcentaje_entrenar,"tot_mpA.csv")), file=fi, append=TRUE, sep = "\n")
cat(sprintf("Ponderar: %8s",ponderar),    file=fi, append=TRUE, sep="\n")
cat(sprintf("nr      : %8d",nr),          file=fi, append=TRUE, sep="\n")
cat(sprintf("eta     : %8.3f",eta),       file=fi, append=TRUE, sep="\n")
cat(sprintf("md      : %8d",md),          file=fi, append=TRUE, sep="\n")
cat(sprintf("cs      : %8.3f",cs),        file=fi, append=TRUE, sep="\n")
cat(sprintf("ss      : %8d",ss),          file=fi, append=TRUE, sep="\n")
cat(sprintf("cb      : %8d",cb),          file=fi, append=TRUE, sep="\n")
cat("\nAscension recta en horas",         file=fi, append=TRUE, sep="\n")
cat(sprintf("  desde : %8.4f", desde_ar), file=fi, append=TRUE, sep="\n")
cat(sprintf("  hasta : %8.4f", hasta_ar), file=fi, append=TRUE, sep="\n")
cat(sprintf("  incr. : %8.5f", inc_ar),   file=fi, append=TRUE, sep="\n")
cat(sprintf("  iterac: %8d", iter_ar),    file=fi, append=TRUE, sep="\n")
cat("\nDeclinacion en grados",            file=fi, append=TRUE, sep="\n")
cat(sprintf("  desde : %8.4f", desde_de), file=fi, append=TRUE, sep="\n")
cat(sprintf("  hasta : %8.4f", hasta_de), file=fi, append=TRUE, sep="\n")
cat(sprintf("  incr. : %8.5f", inc_de),   file=fi, append=TRUE, sep="\n")
cat(sprintf("  iterac: %8d", iter_de),    file=fi, append=TRUE, sep="\n")
cat("\nDistancia en numero de celdas",    file=fi, append=TRUE, sep="\n")
cat(sprintf("  desde : %8.2f", desde_r),  file=fi, append=TRUE, sep="\n")
cat(sprintf("  hasta : %8.2f", hasta_r),  file=fi, append=TRUE, sep="\n")
cat(sprintf("  incr. : %8.5f", inc_r),    file=fi, append=TRUE, sep="\n")
cat(sprintf("  iterac: %8d", iter_r),     file=fi, append=TRUE, sep="\n")
cat("\nEscala",                           file=fi, append=TRUE, sep="\n")
cat(sprintf("  U min : %8.2f", u_min),    file=fi, append=TRUE, sep="\n")
cat(sprintf("  U max : %8.2f", u_max),    file=fi, append=TRUE, sep="\n")
cat(sprintf("  U inc : %8.2f", u_inc),    file=fi, append=TRUE, sep="\n")
cat(sprintf("  iterac: %8d\n", iter_u),   file=fi, append=TRUE, sep="\n")
cat(sprintf("Pruebas :%9.0f\n",pruebas),  file=fi, append=TRUE, sep="\n")

desde_ar <- desde_ar*pi/12
hasta_ar <- hasta_ar*pi/12
inc_ar   <- inc_ar  *pi/12 
desde_de <- desde_de*pi/180
hasta_de <- hasta_de*pi/180
inc_de   <- inc_de  *pi/180 

op <- options(digits.secs = 6)
momento <- Sys.time()
print(momento)
cat(sprintf("%s\n", momento), file=fi, append=TRUE, sep = "\n")
options(op)
ptm <- proc.time()

if (ponderar){
  logregobj <- function(preds, dtrain) {
    labels <- getinfo(dtrain, "label")
    preds  <- 1/(1 + exp(-preds))
    grad   <- (preds - labels) * ponderae
    hess   <- preds * (1 - preds) * ponderae
    return(list(grad = grad, hess = hess))
  }
  evalerror <- function(preds, dtrain) {
    labels <- getinfo(dtrain, "label")
    ppred  <- ifelse (preds > 0.5, 1, 0)
    dif    <- as.numeric(labels) - as.numeric(ppred)
    err    <- sum(dif * dif * ponderae)/sumapon
    return(list(metric = "error", value = err))
  }
  param <- list(
    objective    = logregobj
    ,eval_metric = evalerror
  )
} else {
  param <- list(
    objective    = "binary:logistic"
    ,eval_metric = "error"
  )
}
set.seed(1235) 
bst <- xgboost(subsample=ss, colsample_bylevel=cs, colsample_bytree=cb, eta=eta, max_depth=md, data=dentrena, params=param, nrounds=nr, print_every_n=1000)

datosprobar   <- read.table(paste(sep="",casob,"_muestra_",porcentaje_entrenar,"tot_mpA.csv"), header=TRUE, sep=";", na.strings="NA", dec=",", strip.white=TRUE)
prueba        <- datosprobar[,c(1,2,3)]
Eprueba       <- as.numeric(datosprobar$E)
cartesianas   <- data.matrix(prueba)
ecuatorialini <- CaE(cartesianas)
ecuatorial    <- ecuatorialini
setDT(prueba)
coordenadasPrueba <- model.matrix(~.+0, data = prueba[,, with=F])
dprueba <- xgb.DMatrix(data = coordenadasPrueba, label = Eprueba)

pred <- predict (bst, dprueba)
pred <- ifelse (pred > 0.5, 1, 0)
acto_ini <- 1 - as.numeric(sum(pred != Eprueba))/length(Eprueba)
cat("       i   Inc-AR   Inc-DE   Inc-r  unidad  Acierto", file=fi, append=TRUE, sep = "\n")
cat("  ------ -------- -------- ------- ------- --------", file=fi, append=TRUE, sep = "\n")
linea <- sprintf("  %6d %8.4f %8.4f %7.2f %7.3f %7.3f%s", 0, 0, 0, 0, 1.0, acto_ini*100, "%")
cat(sprintf("%s", linea), file=fi, append=TRUE, sep = "\n")

hasta_ar   <- desde_ar + inc_ar * iter_ar
hasta_de   <- desde_de + inc_de * iter_de
AR <- seq(desde_ar,hasta_ar,inc_ar)
DE <- seq(desde_de,hasta_de,inc_de)
aciertosXY_mejorRU <- expand.grid(x=AR,y=DE)
aciertosXY_mejorRU$z <- 0.5

da         <- 0
dd         <- 0
dr         <- 0
du         <- 0
acto_mejor <- 0
num_iter   <- 0
pb     <- winProgressBar(title = casoa, min = 1, max = iter_ar, width = 400)
ar_act <- desde_ar - inc_ar
for (i in 1:iter_ar){
  #flush.console()
  ar_act <- ar_act + inc_ar
  
  setWinProgressBar(pb, (i-1), title = paste(casoa," ", round((i-1)/iter_ar*100, 0), "% hecho"))
  
  ecuatorial[,1] <- ecuatorialini[,1] + ar_act
  de_act         <- desde_de - inc_de
  for (j in 1:iter_de){
    de_act <- de_act + inc_de
    
    ecuatorial[,2] <- ecuatorialini[,2] + de_act
    cartesianas    <- EaC(ecuatorial)
    
    # totalizar en r y U

    ru_mejoracierto <- 0
    r_act           <- desde_r - inc_r
    for (k in 1:iter_r){
      r_act  <- r_act + inc_r
      for (u in seq(u_min,u_max,u_inc)){
        num_iter <- num_iter + 1
        
        prueba$x <- (cartesianas[,1] + r_act) / u
        prueba$y <- (cartesianas[,2] + r_act) / u
        prueba$z <- (cartesianas[,3] + r_act) / u

        #prueba$x <- round((cartesianas[,1] + r_act) / u, 0)
        #prueba$y <- round((cartesianas[,2] + r_act) / u, 0)
        #prueba$z <- round((cartesianas[,3] + r_act) / u, 0)
        
        coordenadasPrueba <- model.matrix(~.+0, data = prueba[,, with=F])
        dprueba <- xgb.DMatrix(data = coordenadasPrueba, label = Eprueba)
        
        pred <- predict (bst, dprueba)
        pred <- ifelse (pred > 0.5, 1, 0)
        acto <- 1 - as.numeric(sum(pred != Eprueba))/length(Eprueba)
        
        if (acto > ru_mejoracierto){
          ru_mejoracierto     <- acto
          mejoresRU[j]        <- sprintf("%7.2f %7.3f", r_act, u)
          aciertos_mejorRU[j] <- ru_mejoracierto
          if (acto > acto_mejor){
            da         <- ar_act
            dd         <- de_act
            dr         <- r_act
            du         <- u
            acto_mejor <- acto
            cat(sprintf("  %6d %8.4f %8.4f %7.2f %7.3f %7.3f%s", num_iter, ar_act*12/pi, de_act*180/pi, r_act, u, acto_mejor*100, "%"), file=fi, append=TRUE, sep = "\n")
          }
        }
      }
    }
    nn <- (j - 1) * iter_ar + i
    aciertosXY_mejorRU[nn,"z"] <- ru_mejoracierto
  }
  write.table(aciertos_mejorRU, file=fiacto_mejorRU, append=TRUE, row.names=FALSE, col.names=FALSE, sep=";", dec=",", eol=";")
  cat("", file=fiacto_mejorRU, append=TRUE, sep="\n")
  write.table(mejoresRU, file=fimejorRU, append=TRUE, row.names=FALSE, col.names=FALSE, sep=";", dec=",", eol=";")
  cat("", file=fimejorRU, append=TRUE, sep="\n")
}
cat("  ------ -------- -------- ------- ------- --------", file=fi, append=TRUE, sep = "\n")
cat(sprintf("\n  Mejora Total : %7.3f%s", (acto_mejor - acto_ini)*100, "%"), file=fi, append=TRUE, sep = "\n")
close(pb)
write.table(aciertosXY_mejorRU, file=paste(sep="",casob,porcentaje_entrenar,"_cambiar_Ecuatoriales_grafico3d",desamb,".csv"), row.names=FALSE, col.names=TRUE, sep=";", dec=",")
if (require(lattice)){
  wireframe(
    z~x*y, data=aciertosXY_mejorRU,
    xlab="x", ylab="y",
    main="% Acierto",
    drape=TRUE,
    colorkey=TRUE,
    screen=list(z=-60, x=-60)
  )
}
if (require(rgl)){
  open3d()
  bg3d("white")
  material3d(col = "black")
  pastel <- 0.6
  alto   <- (aciertosXY_mejorRU$z - range(aciertosXY_mejorRU$z)[1]) / diff(range(aciertosXY_mejorRU$z))
  rojo   <- pastel + alto
  verde  <- pastel + 0.5
  azul   <- pastel + 1 - alto
  color  <- rgb(rojo, verde, azul, maxColorValue=1 + pastel)
  persp3d(AR,DE,aciertosXY_mejorRU$z,
    col=color,
    xlab="x", ylab="y", zlab="z", main="% Acierto",axes=TRUE,
    specular="black",
    scale=FALSE,
    ticktype = "detailed"
  )
  par3d(windowRect=c(50,50,800,800))
  view3d(theta=-20,phi=-30)
  snapshot3d(paste(sep="",casob,"_muestra_",porcentaje_entrenar,"cambiar_Ecuatoriales_grafico3d_a",desamb,".png"), top=TRUE)
}

ecuatorial[,1] <- ecuatorialini[,1] + da
ecuatorial[,2] <- ecuatorialini[,2] + dd
cartesianas    <- EaC(ecuatorial)
datosprobar$x  <- (cartesianas[,1] + dr) / du
datosprobar$y  <- (cartesianas[,2] + dr) / du
datosprobar$z  <- (cartesianas[,3] + dr) / du

#datosprobar$x  <- round((cartesianas[,1] + dr) / du, 0)
#datosprobar$y  <- round((cartesianas[,2] + dr) / du, 0)
#datosprobar$z  <- round((cartesianas[,3] + dr) / du, 0)

write.table(datosprobar, file = paste(sep="",casob,"_muestra_",porcentaje_entrenar,"cambiar_Ecuatoriales_muestraP",desamb,".csv"), row.names=FALSE, col.names=TRUE, sep=";", dec=",")

tpo <- proc.time() - ptm
op  <- options(digits.secs = 6)
momento <- Sys.time()
print(momento)
cat(sprintf("\n%s\n", momento), file=fi, append=TRUE, sep = "\n")
options(op)
print(tpo)
cat(sprintf("Tiempo. user: %7.1f  system: %7.1f  elapsed: %7.1f\n",as.numeric(tpo[1]),as.numeric(tpo[2]),as.numeric(tpo[3])), file=fi, append=TRUE, sep = "\n")
print("FIN")
