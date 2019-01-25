rm(list=ls())
library(xgboost)
library(data.table)

# ENCUADRE Y ENFOQUE

# Cambia las coordenadas cartesianasde de las observaciones de prueba, cambiando el origen y la escala.
# Se ENCUADRA cambiando el origen de coordenadas mediante una traslación.
# Se ENFOCA cambiando la escala mediante la aplicación de un factor '1/u_act'
# La unidad de rangos y variaciones es una celda.

desamb <- "_Lucas01"

x_min  <- +5
x_max  <- +7
x_inc  <- 0.05

y_min  <- +19
y_max  <- +22
y_inc  <- 0.05

z_min  <- -19
z_max  <- -16
z_inc  <- 0.05

u_min  <- 1.35
u_max  <- 1.50
u_inc  <- 0.005

tot_resultados <- TRUE     # sepadados por ';'
tot_salidas    <- FALSE    # igual que tot_resultados pero formateados y sin ';'
ponderar       <- TRUE

nr  <- 620
eta <- 0.084
md  <- 13
cs  <- 0.94
ss  <- 1
cb  <- 1

disco   <- "F"
marcas  <- "444"
dmina   <- "600000"
dmaxa   <- "1000000"
anguloa <- "_8_12_15_70_"
dminb   <- "600000"
dmaxb   <- "1000000"
angulob <- "_12_16_15_70_"

casoa   <- paste(sep="",disco,":/Gly/G_1_10_",marcas,anguloa,dmina,"_",dmaxa)
casob   <- paste(sep="",disco,":/Gly/G_1_10_",marcas,angulob,dminb,"_",dmaxb)
porcentaje_entrenar <- ""

datosentrenar <- read.table(paste(sep="",casoa,"_muestra_",porcentaje_entrenar,"tot_meA.csv"), header=TRUE, sep=";", na.strings="NA", dec=",", strip.white=TRUE)
ponderae      <- ifelse(datosentrenar$n == 0, 1, datosentrenar$n)
sumapon       <- sum(ponderae)
entrena       <- datosentrenar[,-c(5)]
Eentrena      <- as.numeric(entrena$E)
setDT(entrena)
coordenadasEntrena  <- model.matrix(~.+0, data=entrena [,-c("E"), with=F]) 
dentrena  <- xgb.DMatrix(data=coordenadasEntrena, label=Eentrena) 

iter_x  <- (x_max - x_min)/x_inc + 1
iter_y  <- (y_max - y_min)/y_inc + 1
iter_z  <- (z_max - z_min)/z_inc + 1
iter_u  <- abs(round((u_max - u_min)/u_inc, 0)) + 1
pruebas <- iter_x*iter_y*iter_z*iter_u

if (tot_resultados){
  aciertos_mejorZU <- array(0,dim=c(iter_y))
  mejoresZU        <- array(0,dim=c(iter_y))
  fiacto_mejorZU   <- paste(sep="",casoa,"_muestra_",porcentaje_entrenar,"_Encuadrar_Enfocar_resultados_Azu",desamb,".csv")
  fimejorZU        <- paste(sep="",casoa,"_muestra_",porcentaje_entrenar,"_Encuadrar_Enfocar_resultados_Mzu",desamb,".csv")
}
fi <- paste(sep="",casoa,"_muestra_",porcentaje_entrenar,"_Encuadrar_Enfocar_log",desamb,".txt")
cat(sprintf("\nEntrenar: %s\nProbar  : %s\n",paste(sep="",casoa,"_muestra_",porcentaje_entrenar,"tot_meA.csv"),paste(sep="",casob,"_muestra_",porcentaje_entrenar,"tot_mpA.csv")), file=fi, append=TRUE, sep="\n")
cat(sprintf("Ponderar   : %9s",ponderar),   file=fi, append=TRUE, sep="\n")
cat(sprintf("nr         : %9d",nr),         file=fi, append=TRUE, sep="\n")
cat(sprintf("eta        : %9.3f",eta),      file=fi, append=TRUE, sep="\n")
cat(sprintf("md         : %9d",md),         file=fi, append=TRUE, sep="\n")
cat(sprintf("cs         : %9.3f",cs),       file=fi, append=TRUE, sep="\n")
cat(sprintf("ss         : %9d",ss),         file=fi, append=TRUE, sep="\n")
cat(sprintf("cb         : %9d\n",cb),       file=fi, append=TRUE, sep="\n")
cat(sprintf("x min      : %9.2f", x_min),   file=fi, append=TRUE, sep="\n")
cat(sprintf("x max      : %9.2f", x_max),   file=fi, append=TRUE, sep="\n")
cat(sprintf("x inc      : %9.2f", x_inc),   file=fi, append=TRUE, sep="\n")
cat(sprintf("iteraciones: %9d\n", iter_x),  file=fi, append=TRUE, sep="\n")
cat(sprintf("y min      : %9.2f", y_min),   file=fi, append=TRUE, sep="\n")
cat(sprintf("y max      : %9.2f", y_max),   file=fi, append=TRUE, sep="\n")
cat(sprintf("y inc      : %9.2f", y_inc),   file=fi, append=TRUE, sep="\n")
cat(sprintf("iteraciones: %9d\n", iter_y),  file=fi, append=TRUE, sep="\n")
cat(sprintf("z min      : %9.2f", z_min),   file=fi, append=TRUE, sep="\n")
cat(sprintf("z max      : %9.2f", z_max),   file=fi, append=TRUE, sep="\n")
cat(sprintf("z inc      : %9.2f", z_inc),   file=fi, append=TRUE, sep="\n")
cat(sprintf("iteraciones: %9d\n", iter_z),  file=fi, append=TRUE, sep="\n")
cat(sprintf("U min      : %9.2f", u_min),   file=fi, append=TRUE, sep="\n")
cat(sprintf("U max      : %9.2f", u_max),   file=fi, append=TRUE, sep="\n")
cat(sprintf("U inc      : %9.2f", u_inc),   file=fi, append=TRUE, sep="\n")
cat(sprintf("iteraciones: %9d\n", iter_u),  file=fi, append=TRUE, sep="\n")
cat(sprintf("Pruebas    : %9.0f\n",pruebas),file=fi, append=TRUE, sep="\n")

if (tot_salidas){
  fitot <- paste(sep="",casoa,"_muestra_",porcentaje_entrenar,"_Encuadrar_Enfocar_resultados_TT",desamb,".csv")
  cat(sprintf("\nEntrenar: %s\nProbar  : %s\n",paste(sep="",casoa,"_muestra_",porcentaje_entrenar,"tot_meA.csv"),paste(sep="",casob,"_muestra_",porcentaje_entrenar,"tot_mpA.csv")), file=fitot, append=TRUE, sep="\n")
  cat(sprintf("Ponderar   : %9s",ponderar),   file=fitot, append=TRUE, sep="\n")
  cat(sprintf("nr         : %9d",nr),         file=fitot, append=TRUE, sep="\n")
  cat(sprintf("eta        : %9.3f",eta),      file=fitot, append=TRUE, sep="\n")
  cat(sprintf("md         : %9d",md),         file=fitot, append=TRUE, sep="\n")
  cat(sprintf("cs         : %9.3f",cs),       file=fitot, append=TRUE, sep="\n")
  cat(sprintf("ss         : %9d",ss),         file=fitot, append=TRUE, sep="\n")
  cat(sprintf("cb         : %9d\n",cb),       file=fitot, append=TRUE, sep="\n")
  cat(sprintf("x min      : %9.2f", x_min),   file=fitot, append=TRUE, sep="\n")
  cat(sprintf("x max      : %9.2f", x_max),   file=fitot, append=TRUE, sep="\n")
  cat(sprintf("x inc      : %9.2f", x_inc),   file=fitot, append=TRUE, sep="\n")
  cat(sprintf("iteraciones: %9d\n", iter_x),  file=fitot, append=TRUE, sep="\n")
  cat(sprintf("y min      : %9.2f", y_min),   file=fitot, append=TRUE, sep="\n")
  cat(sprintf("y max      : %9.2f", y_max),   file=fitot, append=TRUE, sep="\n")
  cat(sprintf("y inc      : %9.2f", y_inc),   file=fitot, append=TRUE, sep="\n")
  cat(sprintf("iteraciones: %9d\n", iter_y),  file=fitot, append=TRUE, sep="\n")
  cat(sprintf("z min      : %9.2f", z_min),   file=fitot, append=TRUE, sep="\n")
  cat(sprintf("z max      : %9.2f", z_max),   file=fitot, append=TRUE, sep="\n")
  cat(sprintf("z inc      : %9.2f", z_inc),   file=fitot, append=TRUE, sep="\n")
  cat(sprintf("iteraciones: %9d\n", iter_z),  file=fitot, append=TRUE, sep="\n")
  cat(sprintf("U min      : %9.2f", u_min),   file=fitot, append=TRUE, sep="\n")
  cat(sprintf("U max      : %9.2f", u_max),   file=fitot, append=TRUE, sep="\n")
  cat(sprintf("U inc      : %9.2f", u_inc),   file=fitot, append=TRUE, sep="\n")
  cat(sprintf("iteraciones: %9d\n", iter_u),  file=fitot, append=TRUE, sep="\n")
  cat(sprintf("Pruebas    : %9.0f\n",pruebas),file=fitot, append=TRUE, sep="\n")
}

op <- options(digits.secs=6)
momento <- Sys.time()
print(momento)
cat(sprintf("%s\n", momento), file=fi, append=TRUE, sep="\n")
options(op)
ptm <- proc.time()

if (ponderar){
  logregobj <- function(preds, dtrain) {
    labels <- getinfo(dtrain, "label")
    preds  <- 1/(1 + exp(-preds))
    grad   <- (preds - labels) * ponderae
    hess   <- preds * (1 - preds) * ponderae
    return(list(grad=grad, hess=hess))
  }
  evalerror <- function(preds, dtrain) {
    labels <- getinfo(dtrain, "label")
    ppred  <- ifelse (preds > 0.5, 1, 0)
    dif    <- as.numeric(labels) - as.numeric(ppred)
    err    <- sum(dif * dif * ponderae)/sumapon
    return(list(metric="error", value=err))
  }
  param <- list(
    objective   =logregobj
    ,eval_metric=evalerror
  )
} else {
  param <- list(
    objective   ="binary:logistic"
    ,eval_metric="error"
  )
}
set.seed(1235) 
bst <- xgboost(subsample=ss, colsample_bylevel=cs, colsample_bytree=cb, eta=eta, max_depth=md, data=dentrena, params=param, nrounds=nr, print_every_n=1000)

datosprobar <- read.table(paste(sep="",casob,"_muestra_",porcentaje_entrenar,"tot_mpA.csv"), header=TRUE, sep=";", na.strings="NA", dec=",", strip.white=TRUE)
prueba      <- datosprobar[,c(1,2,3)]
Eprueba     <- as.numeric(datosprobar$E)
np          <- length(Eprueba)
setDT(prueba)
pruebaini   <- prueba
coordenadasPrueba <- model.matrix(~.+0, data=prueba[,, with=F])
dprueba     <- xgb.DMatrix(data=coordenadasPrueba, label=Eprueba)
pred        <- predict (bst, dprueba)
pred        <- ifelse (pred > 0.5, 1, 0)
acto_ini    <- 1 - as.numeric(sum(pred != Eprueba))/length(Eprueba)
cat("despl-x  despl-y  despl-z   unidad    mejor", file=fi, append=TRUE, sep="\n")
cat("-------  -------  -------  -------  --------", file=fi, append=TRUE, sep="\n")
cat(sprintf("%7.2f  %7.2f  %7.2f  %7.3f  %7.3f%s", 0, 0, 0, 1.0, acto_ini*100, "%"), file=fi, append=TRUE, sep="\n")
if (tot_salidas){
  cat("despl-x  despl-y  despl-z   unidad  acierto", file=fitot, append=TRUE, sep="\n")
  cat("-------  -------  -------  -------  --------", file=fitot, append=TRUE, sep="\n")
  cat(sprintf("%7.2f  %7.2f  %7.2f  %7.3f  %7.5f", 0, 0, 0, 1.0, acto_ini), file=fitot, append=TRUE, sep="\n")
}

dx         <- 0
dy         <- 0
dz         <- 0
du         <- 0
acto_mejor <- 0

X <- seq(x_min,x_max,x_inc)
Y <- seq(y_min,y_max,y_inc)
aciertosXY_mejorZU   <- expand.grid(x=X,y=Y)
aciertosXY_mejorZU$z <- 0.5

xpy <- iter_x * iter_y
pb <- winProgressBar(title=casoa, min=0, max=xpy, width=400)

ix <- 0
for (x_act in X){
  ix <- ix + 1

  iy <- 0
  for (y_act in Y){
    iy <- iy + 1
    
    # totalizar en z y U
    
    zu_mejoracierto <- 0
    for (z_act in seq(z_min,z_max,z_inc)){
      for (u_act in seq(u_min,u_max,u_inc)){
        prueba$x <- (pruebaini$x + x_act) / u_act
        prueba$y <- (pruebaini$y + y_act) / u_act
        prueba$z <- (pruebaini$z + z_act) / u_act
        
        coordenadasPrueba <- model.matrix(~.+0, data=prueba[,, with=F])
        dprueba           <- xgb.DMatrix(data=coordenadasPrueba, label=Eprueba)
        pred              <- predict (bst, dprueba)
        pred              <- ifelse (pred > 0.5, 1, 0)
        acto              <- 1 - as.numeric(sum(pred != Eprueba))/np
        if (tot_salidas){
          cat(sprintf("%7.2f %7.2f %7.2f %7.3f %7.5f", x_act, y_act, z_act, u_act, acto), file=fitot, append=TRUE, sep="\n")
        }
        if (tot_resultados){
          if (acto > zu_mejoracierto){
            zu_mejoracierto      <- acto
            mejoresZU[iy]        <- sprintf("%7.2f %7.3f", z_act, u_act)
            aciertos_mejorZU[iy] <- zu_mejoracierto
          }
        }
        if (acto > acto_mejor){
          dx         <- x_act
          dy         <- y_act
          dz         <- z_act
          du         <- u_act
          acto_mejor <- acto
          if (x_act==x_min){
            lx <- "i"
          } else if (x_act==x_max){
            lx <- "s"
          } else{
            lx <- " "
          }
          if (y_act==y_min){
            ly <- "i"
          } else if (y_act==y_max){
            ly <- "s"
          } else{
            ly <- " "
          }
          if (z_act==z_min){
            lz <- "i"
          } else if (z_act==z_max){
            lz <- "s"
          } else{
            lz <- " "
          }
          if (u_act==u_min){
            lu <- "i"
          } else if (u_act==u_max){
            lu <- "s"
          } else{
            lu <- " "
          }
          cat(sprintf("%7.2f%s %7.2f%s %7.2f%s %7.3f%s %7.3f%s", x_act, lx, y_act, ly, z_act, lz, u_act, lu, acto_mejor*100, "%"), file=fi, append=TRUE, sep="\n")
        }
      }
    }
    if (tot_resultados){
      nn <- (iy - 1) * iter_x + ix
      aciertosXY_mejorZU[nn,"z"] <- zu_mejoracierto
    }
    nact <- (ix-1) * iter_y + iy
    setWinProgressBar(pb, nact, title=paste(casoa," ", round(100*nact/xpy, 0), "% hecho"))
  }
  if (tot_resultados){
    write.table(aciertos_mejorZU, file=fiacto_mejorZU, append=TRUE, row.names=FALSE, col.names=FALSE, sep=";", dec=",", eol=";")
    cat("", file=fiacto_mejorZU, append=TRUE, sep="\n")
    write.table(mejoresZU, file=fimejorZU, append=TRUE, row.names=FALSE, col.names=FALSE, sep=";", dec=",", eol=";")
    cat("", file=fimejorZU, append=TRUE, sep="\n")

  }
}
if (tot_salidas){
  cat("-------  -------  -------  -------  --------", file=fitot, append=TRUE, sep="\n")
}
cat("-------  -------  -------  -------  --------", file=fi, append=TRUE, sep="\n")
cat(sprintf("Mejora Total: %7.3f%s", (acto_mejor - acto_ini)*100, "%"), file=fi, append=TRUE, sep="\n")
close(pb)

if (tot_resultados){
  write.table(aciertosXY_mejorZU, file=paste(sep="",casob,"_muestra_",porcentaje_entrenar,"_Encuadrar_Enfocar_grafico3d",desamb,".csv"), row.names=FALSE, col.names=TRUE, sep=";", dec=",")
  if (require(lattice)){
    wireframe(
      z~x*y, data=aciertosXY_mejorZU,
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
    alto   <- (aciertosXY_mejorZU$z - range(aciertosXY_mejorZU$z)[1]) / diff(range(aciertosXY_mejorZU$z))
    rojo   <- pastel + alto
    verde  <- pastel + 0.5
    azul   <- pastel + 1 - alto
    color  <- rgb(rojo, verde, azul, maxColorValue=1 + pastel)
    persp3d(X,Y,aciertosXY_mejorZU$z,
      col=color,
      xlab="x", ylab="y", zlab="z", main="% Acierto",axes=TRUE,
      specular="black",
      scale=FALSE,
      ticktype = "detailed"
    )
    par3d(windowRect=c(50,50,800,800))
    view3d(theta=-20,phi=-30)
    snapshot3d(paste(sep="",casob,"_muestra_",porcentaje_entrenar,"_Encuadrar_Enfocar_grafico3d",desamb,".png"), top=TRUE)
  }
}

datosprobar$x <- (datosprobar$x + dx) / du
datosprobar$y <- (datosprobar$y + dy) / du
datosprobar$z <- (datosprobar$z + dz) / du
write.table(datosprobar, file=paste(sep="",casob,porcentaje_entrenar,"_Encuadrar_Enfocar_muestraP",desamb,".csv"), row.names=FALSE, col.names=TRUE, sep=";", dec=",")

tpo <- proc.time() - ptm
op  <- options(digits.secs=6)
momento <- Sys.time()
print(momento)
cat(sprintf("\n%s\n", momento), file=fi, append=TRUE, sep="\n")
options(op)
print(tpo)
cat(sprintf("Tiempo. user: %7.1f  system: %7.1f  elapsed: %7.1f\n",as.numeric(tpo[1]),as.numeric(tpo[2]),as.numeric(tpo[3])), file=fi, append=TRUE, sep="\n")
print("FIN")

