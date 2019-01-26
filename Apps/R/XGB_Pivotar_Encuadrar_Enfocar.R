rm(list=ls())
library(xgboost)
library(data.table)

# ENCUADRE Y ENFOQUE

# Cambia las coordenadas cartesianasde de las observaciones de prueba, cambiando la orientación, el origen y la escala.
# Se ENCUADRA:
#   pivotando (cambiando la AR y la declinación) sobre el centro geométrico de la muestra (x media, y media, z media). Esto es cuestionable.
#   cambiando el origen de coordenadas mediante una traslación.
# Se ENFOCA cambiando la escala mediante la aplicación de un factor '1/u'
# Una unidad de distancias y variaciones de las mismas es igual a una celda.

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

desamb <- "_Lucas01"

x_min  <- +5.90
x_max  <- +6.40
x_inc  <- 0.05

y_min  <- +20.30
y_max  <- +20.90
y_inc  <- 0.05

z_min  <- -18.00
z_max  <- -17.40
z_inc  <- 0.05

u_min  <- 1.40
u_max  <- 1.45
u_inc  <- 0.005

# horas

desde_ar <- -0.017777795
hasta_ar <- +0.007777781
inc_ar   <-  0.001111112

# grados

desde_de <- -0.176666669
hasta_de <- -0.043333333
inc_de   <-  0.016666667

# celdas

desde_r  <- +0
hasta_r  <- +0
inc_r    <- 1.0

tot_resultados <- TRUE     # sepadados por ';'
tot_salidas    <- FALSE    # igual que tot_resultados pero formateados y sin ';'
ponderar       <- TRUE

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

iter_x  <- abs(round((x_max - x_min)/x_inc,        0)) + 1
iter_y  <- abs(round((y_max - y_min)/y_inc,        0)) + 1
iter_z  <- abs(round((z_max - z_min)/z_inc,        0)) + 1
iter_u  <- abs(round((u_max - u_min)/u_inc,        0)) + 1
iter_ar <- abs(round((hasta_ar - desde_ar)/inc_ar, 0)) + 1
iter_de <- abs(round((hasta_de - desde_de)/inc_de, 0)) + 1
iter_r  <- abs(round((hasta_r - desde_r)/inc_r,    0)) + 1
pruebas <- iter_x*iter_y*iter_z*iter_u*iter_ar*iter_de*iter_r

if (tot_resultados){
  aciertos_mejor <- array(0,dim=c(iter_y))
  mejores        <- array(0,dim=c(iter_y))
  fiacto_mejor   <- paste(sep="",casoa,"_muestra_",porcentaje_entrenar,"Pivota_Encuadrar_Enfocar_resultados_Azu",desamb,".csv")
  fimejor        <- paste(sep="",casoa,"_muestra_",porcentaje_entrenar,"Pivota_Encuadrar_Enfocar_resultados_Mzu",desamb,".csv")
}
fi <- paste(sep="",casoa,"_muestra_",porcentaje_entrenar,"Pivota_Encuadrar_Enfocar_log",desamb,".txt")
cat(sprintf("\nEntrenar: %s\nProbar  : %s\n",paste(sep="",casoa,"_muestra_",porcentaje_entrenar,"tot_meA.csv"),paste(sep="",casob,"_muestra_",porcentaje_entrenar,"tot_mpA.csv")), file=fi, append=TRUE, sep="\n")
cat(sprintf("Ponderar   : %9s",ponderar),    file=fi, append=TRUE, sep="\n")
cat(sprintf("nr         : %9d",nr),          file=fi, append=TRUE, sep="\n")
cat(sprintf("eta        : %9.3f",eta),       file=fi, append=TRUE, sep="\n")
cat(sprintf("md         : %9d",md),          file=fi, append=TRUE, sep="\n")
cat(sprintf("cs         : %9.3f",cs),        file=fi, append=TRUE, sep="\n")
cat(sprintf("ss         : %9d",ss),          file=fi, append=TRUE, sep="\n")
cat(sprintf("cb         : %9d\n",cb),        file=fi, append=TRUE, sep="\n")
cat(sprintf("x min      : %9.2f", x_min),    file=fi, append=TRUE, sep="\n")
cat(sprintf("x max      : %9.2f", x_max),    file=fi, append=TRUE, sep="\n")
cat(sprintf("x inc      : %9.2f", x_inc),    file=fi, append=TRUE, sep="\n")
cat(sprintf("iteraciones: %9d\n", iter_x),   file=fi, append=TRUE, sep="\n")
cat(sprintf("y min      : %9.2f", y_min),    file=fi, append=TRUE, sep="\n")
cat(sprintf("y max      : %9.2f", y_max),    file=fi, append=TRUE, sep="\n")
cat(sprintf("y inc      : %9.2f", y_inc),    file=fi, append=TRUE, sep="\n")
cat(sprintf("iteraciones: %9d\n", iter_y),   file=fi, append=TRUE, sep="\n")
cat(sprintf("z min      : %9.2f", z_min),    file=fi, append=TRUE, sep="\n")
cat(sprintf("z max      : %9.2f", z_max),    file=fi, append=TRUE, sep="\n")
cat(sprintf("z inc      : %9.2f", z_inc),    file=fi, append=TRUE, sep="\n")
cat(sprintf("iteraciones: %9d\n", iter_z),   file=fi, append=TRUE, sep="\n")
cat(sprintf("U min      : %9.2f", u_min),    file=fi, append=TRUE, sep="\n")
cat(sprintf("U max      : %9.2f", u_max),    file=fi, append=TRUE, sep="\n")
cat(sprintf("U inc      : %9.2f", u_inc),    file=fi, append=TRUE, sep="\n")
cat(sprintf("iteraciones: %9d\n", iter_u),   file=fi, append=TRUE, sep="\n")
cat("Ascension recta en horas",              file=fi, append=TRUE, sep="\n")
cat(sprintf("desde      :%10.6f", desde_ar), file=fi, append=TRUE, sep="\n")
cat(sprintf("hasta      :%10.6f", hasta_ar), file=fi, append=TRUE, sep="\n")
cat(sprintf("incr.      :%10.6f", inc_ar),   file=fi, append=TRUE, sep="\n")
cat(sprintf("iteraciones: %9d\n", iter_ar),  file=fi, append=TRUE, sep="\n")
cat("Declinacion en grados",                 file=fi, append=TRUE, sep="\n")
cat(sprintf("desde      :%10.6f", desde_de), file=fi, append=TRUE, sep="\n")
cat(sprintf("hasta      :%10.6f", hasta_de), file=fi, append=TRUE, sep="\n")
cat(sprintf("incr.      :%10.6f", inc_de),   file=fi, append=TRUE, sep="\n")
cat(sprintf("iteraciones: %9d\n", iter_de),  file=fi, append=TRUE, sep="\n")
cat("Distancia en numero de celdas",         file=fi, append=TRUE, sep="\n")
cat(sprintf("desde      : %9.2f", desde_r),  file=fi, append=TRUE, sep="\n")
cat(sprintf("hasta      : %9.2f", hasta_r),  file=fi, append=TRUE, sep="\n")
cat(sprintf("incr.      : %9.5f", inc_r),    file=fi, append=TRUE, sep="\n")
cat(sprintf("iteraciones: %9d\n", iter_r),   file=fi, append=TRUE, sep="\n")
cat(sprintf("Pruebas    : %9.0f\n",pruebas), file=fi, append=TRUE, sep="\n")

if (tot_salidas){
  fitot <- paste(sep="",casoa,"_muestra_",porcentaje_entrenar,"Pivota_Encuadrar_Enfocar_resultados_TT",desamb,".csv")
  cat(sprintf("Ponderar   : %9s",ponderar),    file=fitot, append=TRUE, sep="\n")
  cat(sprintf("nr         : %9d",nr),          file=fitot, append=TRUE, sep="\n")
  cat(sprintf("eta        : %9.3f",eta),       file=fitot, append=TRUE, sep="\n")
  cat(sprintf("md         : %9d",md),          file=fitot, append=TRUE, sep="\n")
  cat(sprintf("cs         : %9.3f",cs),        file=fitot, append=TRUE, sep="\n")
  cat(sprintf("ss         : %9d",ss),          file=fitot, append=TRUE, sep="\n")
  cat(sprintf("cb         : %9d\n",cb),        file=fitot, append=TRUE, sep="\n")
  cat(sprintf("x min      : %9.2f", x_min),    file=fitot, append=TRUE, sep="\n")
  cat(sprintf("x max      : %9.2f", x_max),    file=fitot, append=TRUE, sep="\n")
  cat(sprintf("x inc      : %9.2f", x_inc),    file=fitot, append=TRUE, sep="\n")
  cat(sprintf("iteraciones: %9d\n", iter_x),   file=fitot, append=TRUE, sep="\n")
  cat(sprintf("y min      : %9.2f", y_min),    file=fitot, append=TRUE, sep="\n")
  cat(sprintf("y max      : %9.2f", y_max),    file=fitot, append=TRUE, sep="\n")
  cat(sprintf("y inc      : %9.2f", y_inc),    file=fitot, append=TRUE, sep="\n")
  cat(sprintf("iteraciones: %9d\n", iter_y),   file=fitot, append=TRUE, sep="\n")
  cat(sprintf("z min      : %9.2f", z_min),    file=fitot, append=TRUE, sep="\n")
  cat(sprintf("z max      : %9.2f", z_max),    file=fitot, append=TRUE, sep="\n")
  cat(sprintf("z inc      : %9.2f", z_inc),    file=fitot, append=TRUE, sep="\n")
  cat(sprintf("iteraciones: %9d\n", iter_z),   file=fitot, append=TRUE, sep="\n")
  cat(sprintf("U min      : %9.2f", u_min),    file=fitot, append=TRUE, sep="\n")
  cat(sprintf("U max      : %9.2f", u_max),    file=fitot, append=TRUE, sep="\n")
  cat(sprintf("U inc      : %9.2f", u_inc),    file=fitot, append=TRUE, sep="\n")
  cat(sprintf("iteraciones: %9d\n", iter_u),   file=fitot, append=TRUE, sep="\n")
  cat("Ascension recta en horas",              file=fitot, append=TRUE, sep="\n")
  cat(sprintf("desde      :%10.6f", desde_ar), file=fitot, append=TRUE, sep="\n")
  cat(sprintf("hasta      :%10.6f", hasta_ar), file=fitot, append=TRUE, sep="\n")
  cat(sprintf("incr.      :%10.6f", inc_ar),   file=fitot, append=TRUE, sep="\n")
  cat(sprintf("iteraciones: %9d\n", iter_ar),  file=fitot, append=TRUE, sep="\n")
  cat("Declinacion en grados",                 file=fitot, append=TRUE, sep="\n")
  cat(sprintf("desde      :%10.6f", desde_de), file=fitot, append=TRUE, sep="\n")
  cat(sprintf("hasta      :%10.6f", hasta_de), file=fitot, append=TRUE, sep="\n")
  cat(sprintf("incr.      :%10.6f", inc_de),   file=fitot, append=TRUE, sep="\n")
  cat(sprintf("iteraciones: %9d\n", iter_de),  file=fitot, append=TRUE, sep="\n")
  cat("Distancia en numero de celdas",         file=fitot, append=TRUE, sep="\n")
  cat(sprintf("desde      : %9.2f", desde_r),  file=fitot, append=TRUE, sep="\n")
  cat(sprintf("hasta      : %9.2f", hasta_r),  file=fitot, append=TRUE, sep="\n")
  cat(sprintf("incr.      : %9.5f", inc_r),    file=fitot, append=TRUE, sep="\n")
  cat(sprintf("iteraciones: %9d\n", iter_r),   file=fitot, append=TRUE, sep="\n")
  cat(sprintf("Pruebas    : %9.0f\n",pruebas), file=fitot, append=TRUE, sep="\n")
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

# Centrar

dxm = sum(datosprobar$x)/length(datosprobar$x)
dym = sum(datosprobar$y)/length(datosprobar$y)
dzm = sum(datosprobar$z)/length(datosprobar$z)
datosprobar$x <- datosprobar$x - dxm
datosprobar$y <- datosprobar$y - dym
datosprobar$z <- datosprobar$z - dzm

prueba        <- datosprobar[,c(1,2,3)]
Eprueba       <- as.numeric(datosprobar$E)
np            <- length(Eprueba)
cartesianas   <- data.matrix(prueba)
ecuatorialini <- CaE(cartesianas)
ecuatorial    <- ecuatorialini
setDT(prueba)

# Encuadre para probar el caso inicial

prueba$x <- prueba$x + dxm
prueba$y <- prueba$y + dym
prueba$z <- prueba$z + dzm

coordenadasPrueba <- model.matrix(~.+0, data=prueba[,, with=F])
dprueba     <- xgb.DMatrix(data=coordenadasPrueba, label=Eprueba)
pred        <- predict (bst, dprueba)
pred        <- ifelse (pred > 0.5, 1, 0)
acto_ini    <- 1 - as.numeric(sum(pred != Eprueba))/length(Eprueba)
cat("despl-x  despl-y  despl-z   unidad         Inc-AR      Inc-DE    Inc-r     mejor", file=fi, append=TRUE, sep="\n")
cat("-------  -------  -------  -------     ----------  ----------  -------  --------", file=fi, append=TRUE, sep="\n")
cat(sprintf("%7.2f  %7.2f  %7.2f  %7.3f     %10.6f  %10.6f  %7.2f  %7.3f%s", 0, 0, 0, 1.0, 0, 0, 0, acto_ini*100, "%"), file=fi, append=TRUE, sep="\n")
if (tot_salidas){
  cat("despl-x despl-y despl-z  unidad     Inc-AR     Inc-DE   Inc-r  acierto", file=fitot, append=TRUE, sep="\n")
  cat("------- ------- ------- ------- ---------- ---------- ------- --------", file=fitot, append=TRUE, sep="\n")
  cat(sprintf("%7.2f %7.2f %7.2f %7.3f %10.6f %10.6f %7.2f %7.5f", 0, 0, 0, 1.0, 0, 0, 0, acto_ini), file=fitot, append=TRUE, sep="\n")
}

dx         <- 0
dy         <- 0
dz         <- 0
du         <- 0
da         <- 0
dd         <- 0
dr         <- 0
acto_mejor <- 0

# Angulos en radianes

desde_ar <- desde_ar* pi/12
hasta_ar <- hasta_ar* pi/12
inc_ar   <- inc_ar  * pi/12 
desde_de <- desde_de* pi/180
hasta_de <- hasta_de* pi/180
inc_de   <- inc_de  * pi/180 

X <- seq(x_min,x_max,x_inc)
Y <- seq(y_min,y_max,y_inc)
AR <- seq(desde_ar,hasta_ar,inc_ar)
DE <- seq(desde_de,hasta_de,inc_de)
R  <- seq(desde_r,hasta_r,inc_r)

aciertosXY_mejor   <- expand.grid(x=X,y=Y)
aciertosXY_mejor$z <- 0.5

xpy <- iter_x * iter_y
pb <- winProgressBar(title=casoa, min=0, max=xpy, width=400)

ix <- 0
for (x_act in X){
  ix <- ix + 1
  iy <- 0
  for (y_act in Y){
    iy <- iy + 1
    
    # Totalizar en z, U, AR, 'de' y r
    
    mejoracierto <- 0
    for (z_act in seq(z_min,z_max,z_inc)){
      for (u_act in seq(u_min,u_max,u_inc)){
        for (ar_act in AR){
          
          # Pivotar AR
        
          ecuatorial[,1] <- ecuatorialini[,1] + ar_act
          for (de_act in DE){
          
            # Pivotar declinación
          
            ecuatorial[,2] <- ecuatorialini[,2] + de_act
            cartesianas    <- EaC(ecuatorial)
            for (r_act in R){
              
              # Desplazar r, encuadrar y enfocar
              
              vx <- r_act + dxm + x_act
              vy <- r_act + dym + y_act
              vz <- r_act + dzm + z_act
              prueba$x <- (cartesianas[,1] + vx) / u_act
              prueba$y <- (cartesianas[,2] + vy) / u_act
              prueba$z <- (cartesianas[,3] + vz) / u_act
              
              coordenadasPrueba <- model.matrix(~.+0, data = prueba[,, with=F])
              dprueba <- xgb.DMatrix(data = coordenadasPrueba, label = Eprueba)
        
              pred <- predict (bst, dprueba)
              pred <- ifelse (pred > 0.5, 1, 0)
              acto <- 1 - as.numeric(sum(pred != Eprueba))/np
              if (tot_salidas){
                cat(sprintf("%7.2f %7.2f %7.2f %7.3f %10.6f %10.6f %7.2f %7.5f", x_act, y_act, z_act, u_act, ar_act*12/pi, de_act*180/pi, r_act, acto), file=fitot, append=TRUE, sep="\n")
              }
              if (tot_resultados){
                if (acto > mejoracierto){
                  mejoracierto       <- acto
                  mejores[iy]        <- sprintf("%7.2f %7.3f %10.6f %10.6f %7.2f", z_act, u_act, ar_act*12/pi, de_act*180/pi, r_act)
                  aciertos_mejor[iy] <- mejoracierto
                }
              }
              if (acto > acto_mejor){
                dx         <- x_act
                dy         <- y_act
                dz         <- z_act
                du         <- u_act
                da         <- ar_act
                dd         <- de_act
                dr         <- r_act
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

                if (ar_act==desde_ar){
                  lar <- "i"
                } else if (ar_act==hasta_ar){
                  lar <- "s"
                } else{
                  lar <- " "
                }
                if (de_act==desde_de){
                  lde <- "i"
                } else if (de_act==hasta_de){
                  lde <- "s"
                } else{
                  lde <- " "
                }
                if (r_act==desde_r){
                  lr <- "i"
                } else if (r_act==hasta_r){
                  lr <- "s"
                } else{
                  lr <- " "
                }
                cat(sprintf("%7.2f%s %7.2f%s %7.2f%s %7.3f%s    %10.6f%s %10.6f%s %7.2f%s %7.3f%s", x_act, lx, y_act, ly, z_act, lz, u_act, lu, ar_act*12/pi, lar, de_act*180/pi, lde, r_act, lr, acto_mejor*100, "%"), file=fi, append=TRUE, sep="\n")
              }
            }
          }
        }
      }
    }
    if (tot_resultados){
      nn <- (iy - 1) * iter_x + ix
      aciertosXY_mejor[nn,"z"] <- mejoracierto
    }
    nact <- (ix-1) * iter_y + iy
    setWinProgressBar(pb, nact, title=paste(casoa," ", round(100*nact/xpy, 0), "% hecho"))
  }
  if (tot_resultados){
    write.table(aciertos_mejor, file=fiacto_mejor, append=TRUE, row.names=FALSE, col.names=FALSE, sep=";", dec=",", eol=";")
    cat("", file=fiacto_mejor, append=TRUE, sep="\n")
    write.table(mejores, file=fimejor, append=TRUE, row.names=FALSE, col.names=FALSE, sep=";", dec=",", eol=";")
    cat("", file=fimejor, append=TRUE, sep="\n")
  }
}
if (tot_salidas){
  cat("------- ------- ------- ------- ---------- ---------- ------- --------", file=fitot, append=TRUE, sep="\n")
}
cat("-------  -------  -------  -------     ----------  ----------  -------  --------", file=fi, append=TRUE, sep="\n")
cat(sprintf("Mejora Total: %7.3f%s", (acto_mejor - acto_ini)*100, "%"), file=fi, append=TRUE, sep="\n")
close(pb)

if (tot_resultados){
  write.table(aciertosXY_mejor, file=paste(sep="",casob,"_muestra_",porcentaje_entrenar,"Pivota_Encuadrar_Enfocar_grafico3d",desamb,".csv"), row.names=FALSE, col.names=TRUE, sep=";", dec=",")
  if (require(lattice)){
    wireframe(
      z~x*y, data=aciertosXY_mejor,
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
    alto   <- (aciertosXY_mejor$z - range(aciertosXY_mejor$z)[1]) / diff(range(aciertosXY_mejor$z))
    rojo   <- pastel + alto
    verde  <- pastel + 0.5
    azul   <- pastel + 1 - alto
    color  <- rgb(rojo, verde, azul, maxColorValue=1 + pastel)
    persp3d(X,Y,aciertosXY_mejor$z,
      col=color,
      xlab="x", ylab="y", zlab="z", main="% Acierto",axes=TRUE,
      specular="black",
      scale=FALSE,
      ticktype = "detailed"
    )
    par3d(windowRect=c(50,50,800,800))
    view3d(theta=-20,phi=-30)
    snapshot3d(paste(sep="",casob,"_muestra_",porcentaje_entrenar,"Pivota_Encuadrar_Enfocar_grafico3d_a",desamb,".png"), top=TRUE)
  }
}

ecuatorial[,1] <- ecuatorialini[,1] + da
ecuatorial[,2] <- ecuatorialini[,2] + dd
cartesianas    <- EaC(ecuatorial)
datosprobar$x  <- cartesianas[,1] + dr
datosprobar$y  <- cartesianas[,2] + dr
datosprobar$z  <- cartesianas[,3] + dr
datosprobar$x  <- (datosprobar$x + dxm + dx) / du
datosprobar$y  <- (datosprobar$y + dym + dy) / du
datosprobar$z  <- (datosprobar$z + dzm + dz) / du
write.table(datosprobar, file=paste(sep="",casob,porcentaje_entrenar,"Pivota_Encuadrar_Enfocar_muestraP",desamb,".csv"), row.names=FALSE, col.names=TRUE, sep=";", dec=",")

tpo <- proc.time() - ptm
op  <- options(digits.secs=6)
momento <- Sys.time()
print(momento)
cat(sprintf("\n%s\n", momento), file=fi, append=TRUE, sep="\n")
options(op)
print(tpo)
cat(sprintf("Tiempo. user: %7.1f  system: %7.1f  elapsed: %7.1f\n",as.numeric(tpo[1]),as.numeric(tpo[2]),as.numeric(tpo[3])), file=fi, append=TRUE, sep="\n")
print("FIN")

