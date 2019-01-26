rm(list = ls())
library(ggplot2)
library(xgboost)
library(caret)
library(data.table)

normaliza =  function(df, rangos){
  as.data.frame(
    Map(function(columna, rango){
      # En este caso, los atributos constantes los dejamos sin modificar
      if(is.numeric(columna) & rango[2]>rango[1]) 
        (columna-rango[1])/(rango[2]-rango[1])
      else columna},
      df,
      rangos)
  )
}

ponderar <- TRUE
nr  <- 620
eta <- 0.084
md  <- 13
cs  <- 0.94
ss  <- 1
cb  <- 1

disco   <- "C"
tipomta <- "tot"
marcas  <- "444"
dmina   <- "600000"
dmaxa   <- "1000000"
anguloa <- "_8_12_15_70_"
dminb   <- "600000"
dmaxb   <- "1000000"
angulob <- "_12_16_15_70_"

normalizar_estandar <- FALSE

normalizar_gly      <- FALSE
xnorme <-  0.0
ynorme <-  0.0
znorme <-  0.0
unorme <-  1.0
xnormp <-  0.0
ynormp <-  0.0
znormp <-  0.0
unormp <-  1.0

casoa   <- paste(sep="",disco,":/Gly/G_1_10_",marcas,anguloa,dmina,"_",dmaxa)
casob   <- paste(sep="",disco,":/Gly/G_1_10_",marcas,angulob,dminb,"_",dmaxb)
porcentaje_entrenar <- "60_"

datosentrenar <- read.table(paste(sep="",casoa,"_muestra_",porcentaje_entrenar,tipomta,"_meA.csv"), header=TRUE, sep=";", na.strings="NA", dec=",", strip.white=TRUE)
ponderae <- ifelse(datosentrenar$n == 0, 1, datosentrenar$n)
sumapon  <- sum(ponderae)
ponderae_n <- ponderae
sumapon_n <- sumapon
datosprobar   <- read.table(paste(sep="",casob,"_muestra_",porcentaje_entrenar,tipomta,"_mpA.csv"), header=TRUE, sep=";", na.strings="NA", dec=",", strip.white=TRUE)
if (normalizar_estandar){
  rangos = lapply(datosentrenar, function(columna) if(is.numeric(columna)) range(columna) else c(NA,NA))
  print(rangos)
  datosentrenar <- normaliza(datosentrenar, rangos)
  datosprobar   <- normaliza(datosprobar, rangos)
}
if (normalizar_gly){
  datosentrenar$x <- (datosentrenar$x - xnorme) / unorme
  datosentrenar$y <- (datosentrenar$y - ynorme) / unorme
  datosentrenar$z <- (datosentrenar$z - znorme) / unorme

  datosprobar$x   <- (datosprobar$x - xnormp) / unormp
  datosprobar$y   <- (datosprobar$y - ynormp) / unormp
  datosprobar$z   <- (datosprobar$z - znormp) / unormp
}
entrena  <- datosentrenar[,-c(5)]
Eentrena  <- as.numeric(entrena$E)
setDT(entrena)
coordenadasEntrena  <- model.matrix(~.+0, data = entrena [,-c("E"), with=F]) 
#dentrena  <- xgb.DMatrix(data = coordenadasEntrena, label = Eentrena, weight = ponderae) 
dentrena  <- xgb.DMatrix(data = coordenadasEntrena, label = Eentrena) 

prueba   <- datosprobar  [,-c(5)]
setDT(prueba)
Eprueba   <- as.numeric(prueba$E)
coordenadasPrueba   <- model.matrix(~.+0, data = prueba  [,-c("E"), with=F])
dprueba   <- xgb.DMatrix(data = coordenadasPrueba,  label = Eprueba)

if (ponderar){
  # La funcion objetivo debe retornar la primera y segunda derivada (gradiente y hessiano)
  # Esta es equivalente a 'binary:logistic' salvo la ponderación
  logregobj <- function(preds, dtrain) {
    labels <- getinfo(dtrain, "label")
    preds  <- 1/(1 + exp(-preds))
    grad   <- (preds - labels) * ponderae
    hess   <- preds * (1 - preds) * ponderae
    return(list(grad = grad, hess = hess))
  }
  #
  # 'eval_metric' no interviene en la optimización, sólo evalua a posteriori
  #
  # Como 'logregobj' devuelve las predicciones con la función logística ya aplicada hay
  # que reformular 'evalerror' porque la de fábrica la espera antes de la la función logística
  evalerror <- function(preds, dtrain) {
    labels <- getinfo(dtrain, "label")
    ppred  <- ifelse (preds > 0.5, 1, 0)
#    err <- as.numeric(sum(labels != (preds > 0)))/length(labels)
#    err <- as.numeric(sum(ppred != labels))/length(labels)
    dif <- as.numeric(labels) - as.numeric(ppred)
    err <- sum(dif * dif * ponderae)/sumapon
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

#watchlist <- list(train=dentrena, test=dprueba)
#bst <- xgb.train(subsample=ss, colsample_bylevel=cs, colsample_bytree=cb, eta=eta, max_depth=md, data=dentrena, params=param, nrounds=nr, print_every_n=100, watchlist=watchlist)

set.seed(1235) 
bst <- xgboost(subsample=ss, colsample_bylevel=cs, colsample_bytree=cb, eta=eta, max_depth=md, data=dentrena, params=param, nrounds=nr, print_every_n=1000)

importance_matrix <- xgb.importance(model = bst)
print(importance_matrix)
xgb.plot.importance(importance_matrix = importance_matrix)

#require(DiagrammeR)
#xgb.plot.tree(model = bst)

pred <- predict (bst, dentrena)
pred <- ifelse (pred > 0.5, 1, 0)
mc <- confusionMatrix(factor(pred,levels=0:1), factor(Eentrena,levels=0:1))
linea0 <- paste(sep="",
  "Hyperleda",
  sprintf(" nr:%6d",nr),
  sprintf(" eta:%6.3f",eta),
  sprintf(" cs:%5.2f",cs),
  sprintf(" md:%5.1f",md),
  sprintf(" Kappa:%6.3f",unname(mc$overall[2])),
  sprintf(" McnemarPValue:%6.3f",unname(mc$overall[7])),
  sprintf(" Acierto:%7.3f",unname(mc$overall[1])*100),"%",
  sprintf(" Entrena:%6d No:%6d Si:%6d ( %7.3f%s )",nrow(entrena),nrow(entrena)-sum(entrena$E),sum(entrena$E),sum(entrena$E)*100/nrow(entrena),"%"),
  sprintf(" Prueba: %6d No:%6d Si:%6d ( %7.3f%s )",nrow(prueba),nrow(prueba)-sum(prueba$E),sum(prueba$E),sum(prueba$E)*100/nrow(prueba),"%"),
  sprintf(" %6.0f",mc$table[1,1]),
  sprintf(" %6.0f",mc$table[2,2])
)
linea0c <- gsub('\\.',',',paste(sep="",
  "Hyperleda",
  sprintf(" %5.3f",eta),
  sprintf(" %5.2f",cs),
  sprintf(" %4.1f",md),
  sprintf(" %8.5f",unname(mc$overall[2])),
  sprintf(" %7.5f",unname(mc$overall[7])),
  sprintf(" %7.5f",unname(mc$overall[1])),
  sprintf(" %6d %6d %6d %7.5f",nrow(entrena),nrow(entrena)-sum(entrena$E),sum(entrena$E),sum(entrena$E)/nrow(entrena)),
  sprintf(" %6d %6d %6d %7.5f",nrow(prueba),nrow(prueba)-sum(prueba$E),sum(prueba$E),sum(prueba$E)/nrow(prueba)),
  sprintf(" %6d",mc$table[1,1]),
  sprintf(" %6d",mc$table[2,2])
))

pred <- predict (bst, dprueba)
pred <- ifelse (pred > 0.5, 1, 0)
mc <- confusionMatrix(factor(pred,levels=0:1), factor(Eprueba,levels=0:1))
linea1 <- paste(sep="",
  "Hyperleda",
  sprintf(" nr:%6d",nr),
  sprintf(" eta:%6.3f",eta),
  sprintf(" cs:%5.2f",cs),
  sprintf(" md:%5.1f",md),
  sprintf(" Kappa:%6.3f",unname(mc$overall[2])),
  sprintf(" McnemarPValue:%6.3f",unname(mc$overall[7])),
  sprintf(" Acierto:%7.3f",unname(mc$overall[1])*100),"%",
  sprintf(" Entrena:%6d No:%6d Si:%6d ( %7.3f%s )",nrow(entrena),nrow(entrena)-sum(entrena$E),sum(entrena$E),sum(entrena$E)*100/nrow(entrena),"%"),
  sprintf(" Prueba: %6d No:%6d Si:%6d ( %7.3f%s )",nrow(prueba),nrow(prueba)-sum(prueba$E),sum(prueba$E),sum(prueba$E)*100/nrow(prueba),"%"),
  sprintf(" %6.0f",mc$table[1,1]),
  sprintf(" %6.0f",mc$table[2,2])
)

linea1c <- gsub('\\.',',',paste(sep="",
  "Hyperleda",
  sprintf(" %5.3f",eta),
  sprintf(" %5.2f",cs),
  sprintf(" %4.1f",md),
  sprintf(" %8.5f",unname(mc$overall[2])),
  sprintf(" %7.5f",unname(mc$overall[7])),
  sprintf(" %7.5f",unname(mc$overall[1])),
  sprintf(" %6d %6d %6d %7.5f",nrow(entrena),nrow(entrena)-sum(entrena$E),sum(entrena$E),sum(entrena$E)/nrow(entrena)),
  sprintf(" %6d %6d %6d %7.5f",nrow(prueba),nrow(prueba)-sum(prueba$E),sum(prueba$E),sum(prueba$E)/nrow(prueba)),
  sprintf(" %6d",mc$table[1,1]),
  sprintf(" %6d",mc$table[2,2])
))
datosentrenar <- read.table(paste(sep="",casoa,"_muestra_",porcentaje_entrenar,"tot_az_meA.csv"), header=TRUE, sep=";", na.strings="NA", dec=",", strip.white=TRUE)

# Para muestras aleatorias es incoherente porque el número de galaxias 'porcentaje_entrenar' celda es siempre uno

ponderae <- ifelse(datosentrenar$n == 0, 1, datosentrenar$n)
sumapon  <- sum(ponderae)
ponderae_az <- ponderae
sumapon_az <- sumapon
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
    dif <- as.numeric(labels) - as.numeric(ppred)
    err <- sum(dif * dif * ponderae)/sumapon
    return(list(metric = "error", value = err))
  }
  param <- list(
    objective    = logregobj
    ,eval_metric = evalerror
  )
}
datosprobar   <- read.table(paste(sep="",casob,"_muestra_",porcentaje_entrenar,"tot_az_mpA.csv"), header=TRUE, sep=";", na.strings="NA", dec=",", strip.white=TRUE)
if (normalizar_estandar){
  rangosaz = lapply(datosentrenar, function(columna) if(is.numeric(columna)) range(columna) else c(NA,NA))
  print(rangosaz)
  datosentrenar = normaliza(datosentrenar, rangosaz)
  datosprobar = normaliza(datosprobar, rangosaz)
}
if (normalizar_gly){
  datosentrenar$x <- (datosentrenar$x - xnorme) / unorme
  datosentrenar$y <- (datosentrenar$y - ynorme) / unorme
  datosentrenar$z <- (datosentrenar$z - znorme) / unorme

  datosprobar$x   <- (datosprobar$x - xnormp) / unormp
  datosprobar$y   <- (datosprobar$y - ynormp) / unormp
  datosprobar$z   <- (datosprobar$z - znormp) / unormp
}
entrena <- datosentrenar[,-c(5)]
prueba  <- datosprobar  [,-c(5)]
setDT(entrena)
setDT(prueba)
Eentrena <- entrena$E 
Eprueba <- prueba$E
coordenadasEntrena <- model.matrix(~.+0, data = entrena[,-c("E"), with=F]) 
coordenadasPrueba  <- model.matrix(~.+0, data = prueba [,-c("E"), with=F])
Eentrena <- as.numeric(Eentrena)
Eprueba  <- as.numeric(Eprueba)
dentrena <- xgb.DMatrix(data = coordenadasEntrena, label = Eentrena) 
dprueba  <- xgb.DMatrix(data = coordenadasPrueba, label = Eprueba)
set.seed(1235) 
bstaz <- xgboost(subsample=ss, colsample_bylevel=cs, colsample_bytree=cb, eta=eta, max_depth=md, data=dentrena, params=param, nrounds=nr, print_every_n=1000)
pred <- predict (bstaz, dprueba)
pred <- ifelse (pred > 0.5, 1, 0)
mc <- confusionMatrix(factor(pred,levels=0:1), factor(Eprueba,levels=0:1))
linea2 <- paste(sep="",
  "Azar     ",
  sprintf(" nr:%6d",nr),
  sprintf(" eta:%6.3f",eta),
  sprintf(" cs:%5.2f",cs),
  sprintf(" md:%5.1f",md),
  sprintf(" Kappa:%6.3f",unname(mc$overall[2])),
  sprintf(" McnemarPValue:%6.3f",unname(mc$overall[7])),
  sprintf(" Acierto:%7.3f",unname(mc$overall[1])*100),"%",
  sprintf(" Entrena:%6d No:%6d Si:%6d ( %7.3f%s )",nrow(entrena),nrow(entrena)-sum(entrena$E),sum(entrena$E),sum(entrena$E)*100/nrow(entrena),"%"),
  sprintf(" Prueba: %6d No:%6d Si:%6d ( %7.3f%s )",nrow(prueba),nrow(prueba)-sum(prueba$E),sum(prueba$E),sum(prueba$E)*100/nrow(prueba),"%"),
  sprintf(" %6.0f",mc$table[1,1]),
  sprintf(" %6.0f",mc$table[2,2])
)
linea2c <- gsub('\\.',',',paste(sep="",
  "Azar     ",
  sprintf(" %5.3f",eta),
  sprintf(" %5.2f",cs),
  sprintf(" %4.1f",md),
  sprintf(" %8.5f",unname(mc$overall[2])),
  sprintf(" %7.5f",unname(mc$overall[7])),
  sprintf(" %7.5f",unname(mc$overall[1])),
  sprintf(" %6d %6d %6d %7.5f",nrow(entrena),nrow(entrena)-sum(entrena$E),sum(entrena$E),sum(entrena$E)/nrow(entrena)),
  sprintf(" %6d %6d %6d %7.5f",nrow(prueba),nrow(prueba)-sum(prueba$E),sum(prueba$E),sum(prueba$E)/nrow(prueba)),
  sprintf(" %6d",mc$table[1,1]),
  sprintf(" %6d",mc$table[2,2])
))
pred <- predict (bst, dprueba)
pred <- ifelse (pred > 0.5, 1, 0)
mc <- confusionMatrix(factor(pred,levels=0:1), factor(Eprueba,levels=0:1))
linea2bis <- paste(sep="",
  "H ->Azar ",
  sprintf(" nr:%6d",nr),
  sprintf(" eta:%6.3f",eta),
  sprintf(" cs:%5.2f",cs),
  sprintf(" md:%5.1f",md),
  sprintf(" Kappa:%6.3f",unname(mc$overall[2])),
  sprintf(" McnemarPValue:%6.3f",unname(mc$overall[7])),
  sprintf(" Acierto:%7.3f",unname(mc$overall[1])*100),"%",
  sprintf(" Entrena:%6d No:%6d Si:%6d ( %7.3f%s )",nrow(entrena),nrow(entrena)-sum(entrena$E),sum(entrena$E),sum(entrena$E)*100/nrow(entrena),"%"),
  sprintf(" Prueba: %6d No:%6d Si:%6d ( %7.3f%s )",nrow(prueba),nrow(prueba)-sum(prueba$E),sum(prueba$E),sum(prueba$E)*100/nrow(prueba),"%"),
  sprintf(" %6.0f",mc$table[1,1]),
  sprintf(" %6.0f",mc$table[2,2])
)
linea2cbis <- gsub('\\.',',',paste(sep="",
  "H ->Azar ",
  sprintf(" %5.3f",eta),
  sprintf(" %5.2f",cs),
  sprintf(" %4.1f",md),
  sprintf(" %8.5f",unname(mc$overall[2])),
  sprintf(" %7.5f",unname(mc$overall[7])),
  sprintf(" %7.5f",unname(mc$overall[1])),
  sprintf(" %6d %6d %6d %7.5f",nrow(entrena),nrow(entrena)-sum(entrena$E),sum(entrena$E),sum(entrena$E)/nrow(entrena)),
  sprintf(" %6d %6d %6d %7.5f",nrow(prueba),nrow(prueba)-sum(prueba$E),sum(prueba$E),sum(prueba$E)/nrow(prueba)),
  sprintf(" %6d",mc$table[1,1]),
  sprintf(" %6d",mc$table[2,2])
))


datosentrenar <- read.table(paste(sep="",casoa,"_muestra_",porcentaje_entrenar,"tot_azaj_meA.csv"), header=TRUE, sep=";", na.strings="NA", dec=",", strip.white=TRUE)

# Para muestras aleatorias es incoherente porque el número de galaxias porcentaje_entrenar celda es siempre uno

ponderae <- ifelse(datosentrenar$n == 0, 1, datosentrenar$n)
sumapon  <- sum(ponderae)
ponderae_azaj <- ponderae
sumapon_azaj <- sumapon
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
    dif <- as.numeric(labels) - as.numeric(ppred)
    err <- sum(dif * dif * ponderae)/sumapon
    return(list(metric = "error", value = err))
  }
  param <- list(
    objective    = logregobj
    ,eval_metric = evalerror
  )
}
datosprobar   <- read.table(paste(sep="",casob,"_muestra_",porcentaje_entrenar,"tot_azaj_mpA.csv"), header=TRUE, sep=";", na.strings="NA", dec=",", strip.white=TRUE)
if (normalizar_estandar){
  rangosazaj = lapply(datosentrenar, function(columna) if(is.numeric(columna)) range(columna) else c(NA,NA))
  print(rangosazaj)
  datosentrenar = normaliza(datosentrenar, rangosazaj)
  datosprobar = normaliza(datosprobar, rangosazaj)
}
if (normalizar_gly){
  datosentrenar$x <- (datosentrenar$x - xnorme) / unorme
  datosentrenar$y <- (datosentrenar$y - ynorme) / unorme
  datosentrenar$z <- (datosentrenar$z - znorme) / unorme

  datosprobar$x   <- (datosprobar$x - xnormp) / unormp
  datosprobar$y   <- (datosprobar$y - ynormp) / unormp
  datosprobar$z   <- (datosprobar$z - znormp) / unormp
}
entrena <- datosentrenar[,-c(5)]
prueba  <- datosprobar  [,-c(5)]
setDT(entrena)
setDT(prueba)
Eentrena <- entrena$E 
Eprueba <- prueba$E
coordenadasEntrena <- model.matrix(~.+0, data = entrena[,-c("E"), with=F]) 
coordenadasPrueba  <- model.matrix(~.+0, data = prueba [,-c("E"), with=F])
Eentrena <- as.numeric(Eentrena)
Eprueba  <- as.numeric(Eprueba)
dentrena <- xgb.DMatrix(data = coordenadasEntrena, label = Eentrena) 
dprueba  <- xgb.DMatrix(data = coordenadasPrueba, label = Eprueba)
set.seed(1235) 
bstazaj <- xgboost(subsample=ss, colsample_bylevel=cs, colsample_bytree=cb, eta=eta, max_depth=md, data=dentrena, params=param, nrounds=nr, print_every_n=1000)
pred <- predict (bstazaj, dprueba)
pred <- ifelse (pred > 0.5, 1, 0)
mc <- confusionMatrix(factor(pred,levels=0:1), factor(Eprueba,levels=0:1))
linea3 <- paste(sep="",
  "Azar ajus",
  sprintf(" nr:%6d",nr),
  sprintf(" eta:%6.3f",eta),
  sprintf(" cs:%5.2f",cs),
  sprintf(" md:%5.1f",md),
  sprintf(" Kappa:%6.3f",unname(mc$overall[2])),
  sprintf(" McnemarPValue:%6.3f",unname(mc$overall[7])),
  sprintf(" Acierto:%7.3f",unname(mc$overall[1])*100),"%",
  sprintf(" Entrena:%6d No:%6d Si:%6d ( %7.3f%s )",nrow(entrena),nrow(entrena)-sum(entrena$E),sum(entrena$E),sum(entrena$E)*100/nrow(entrena),"%"),
  sprintf(" Prueba: %6d No:%6d Si:%6d ( %7.3f%s )",nrow(prueba),nrow(prueba)-sum(prueba$E),sum(prueba$E),sum(prueba$E)*100/nrow(prueba),"%"),
  sprintf(" %6.0f",mc$table[1,1]),
  sprintf(" %6.0f",mc$table[2,2])
)
linea3c <- gsub('\\.',',',paste(sep="",
  "Azar ajus",
  sprintf(" %5.3f",eta),
  sprintf(" %5.2f",cs),
  sprintf(" %4.1f",md),
  sprintf(" %8.5f",unname(mc$overall[2])),
  sprintf(" %7.5f",unname(mc$overall[7])),
  sprintf(" %7.5f",unname(mc$overall[1])),
  sprintf(" %6d %6d %6d %7.5f",nrow(entrena),nrow(entrena)-sum(entrena$E),sum(entrena$E),sum(entrena$E)/nrow(entrena)),
  sprintf(" %6d %6d %6d %7.5f",nrow(prueba),nrow(prueba)-sum(prueba$E),sum(prueba$E),sum(prueba$E)/nrow(prueba)),
  sprintf(" %6d",mc$table[1,1]),
  sprintf(" %6d",mc$table[2,2])
))
pred <- predict (bst, dprueba)
pred <- ifelse (pred > 0.5, 1, 0)
mc <- confusionMatrix(factor(pred,levels=0:1), factor(Eprueba,levels=0:1))
linea3bis <- paste(sep="",
  "H ->AzarJ",
  sprintf(" nr:%6d",nr),
  sprintf(" eta:%6.3f",eta),
  sprintf(" cs:%5.2f",cs),
  sprintf(" md:%5.1f",md),
  sprintf(" Kappa:%6.3f",unname(mc$overall[2])),
  sprintf(" McnemarPValue:%6.3f",unname(mc$overall[7])),
  sprintf(" Acierto:%7.3f",unname(mc$overall[1])*100),"%",
  sprintf(" Entrena:%6d No:%6d Si:%6d ( %7.3f%s )",nrow(entrena),nrow(entrena)-sum(entrena$E),sum(entrena$E),sum(entrena$E)*100/nrow(entrena),"%"),
  sprintf(" Prueba: %6d No:%6d Si:%6d ( %7.3f%s )",nrow(prueba),nrow(prueba)-sum(prueba$E),sum(prueba$E),sum(prueba$E)*100/nrow(prueba),"%"),
  sprintf(" %6.0f",mc$table[1,1]),
  sprintf(" %6.0f",mc$table[2,2])
)
linea3cbis <- gsub('\\.',',',paste(sep="",
  "H ->AzarJ",
  sprintf(" %5.3f",eta),
  sprintf(" %5.2f",cs),
  sprintf(" %4.1f",md),
  sprintf(" %8.5f",unname(mc$overall[2])),
  sprintf(" %7.5f",unname(mc$overall[7])),
  sprintf(" %7.5f",unname(mc$overall[1])),
  sprintf(" %6d %6d %6d %7.5f",nrow(entrena),nrow(entrena)-sum(entrena$E),sum(entrena$E),sum(entrena$E)/nrow(entrena)),
  sprintf(" %6d %6d %6d %7.5f",nrow(prueba),nrow(prueba)-sum(prueba$E),sum(prueba$E),sum(prueba$E)/nrow(prueba)),
  sprintf(" %6d",mc$table[1,1]),
  sprintf(" %6d",mc$table[2,2])
))

# Total muestra sector b

ponderae    <- ponderae_n
sumapon     <- sumapon_n
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
    dif <- as.numeric(labels) - as.numeric(ppred)
    err <- sum(dif * dif * ponderae)/sumapon
    return(list(metric = "error", value = err))
  }
  param <- list(
    objective    = logregobj
    ,eval_metric = evalerror
  )
}
datosprobar <- read.table(paste(sep="",casob,"_muestra_tot_mzA.csv"), header=TRUE, sep=";", na.strings="NA", dec=",", strip.white=TRUE)
if (normalizar_estandar){
  datosprobar = normaliza(datosprobar, rangos)
}
if (normalizar_gly){
  datosprobar$x   <- (datosprobar$x - xnormp) / unormp
  datosprobar$y   <- (datosprobar$y - ynormp) / unormp
  datosprobar$z   <- (datosprobar$z - znormp) / unormp
}
prueba <- datosprobar[,-c(5)]
setDT(prueba)
Eprueba <- prueba$E
coordenadasprueba <- model.matrix(~.+0, data = prueba[,-c("E"), with=F])
Eprueba <- as.numeric(Eprueba)
dprueba <- xgb.DMatrix(data = coordenadasprueba, label = Eprueba)
pred <- predict (bst, dprueba)
pred <- ifelse (pred > 0.5, 1, 0)
mc <- confusionMatrix(factor(pred,levels=0:1), factor(Eprueba,levels=0:1))
lineat1 <- paste(sep="",
  "Hyperleda",
  sprintf(" nr:%6d",nr),
  sprintf(" eta:%6.3f",eta),
  sprintf(" cs:%5.2f",cs),
  sprintf(" md:%5.1f",md),
  sprintf(" Kappa:%6.3f",unname(mc$overall[2])),
  sprintf(" McnemarPValue:%6.3f",unname(mc$overall[7])),
  sprintf(" Acierto:%7.3f",unname(mc$overall[1])*100),"%",
  sprintf(" Entrena:%6d No:%6d Si:%6d ( %7.3f%s )",nrow(entrena),nrow(entrena)-sum(entrena$E),sum(entrena$E),sum(entrena$E)*100/nrow(entrena),"%"),
  sprintf(" Prueba: %6d No:%6d Si:%6d ( %7.3f%s )",nrow(prueba),nrow(prueba)-sum(prueba$E),sum(prueba$E),sum(prueba$E)*100/nrow(prueba),"%"),
  sprintf(" %6.0f",mc$table[1,1]),
  sprintf(" %6.0f",mc$table[2,2])
)
lineat1c <- gsub('\\.',',',paste(sep="",
  "Hyperleda",
  sprintf(" %5.3f",eta),
  sprintf(" %5.2f",cs),
  sprintf(" %4.1f",md),
  sprintf(" %8.5f",unname(mc$overall[2])),
  sprintf(" %7.5f",unname(mc$overall[7])),
  sprintf(" %7.5f",unname(mc$overall[1])),
  sprintf(" %6d %6d %6d %7.5f",nrow(entrena),nrow(entrena)-sum(entrena$E),sum(entrena$E),sum(entrena$E)/nrow(entrena)),
  sprintf(" %6d %6d %6d %7.5f",nrow(prueba),nrow(prueba)-sum(prueba$E),sum(prueba$E),sum(prueba$E)/nrow(prueba)),
  sprintf(" %6d",mc$table[1,1]),
  sprintf(" %6d",mc$table[2,2])
))

# Galaxias sector b

datosprobar <- read.table(paste(sep="",casob,"_galaxias_intA.csv"), header=TRUE, sep=";", na.strings="NA", dec=",", strip.white=TRUE)
if (normalizar_estandar){
  datosprobar = normaliza(datosprobar, rangos)
}
if (normalizar_gly){
  datosprobar$x   <- (datosprobar$x - xnormp) / unormp
  datosprobar$y   <- (datosprobar$y - ynormp) / unormp
  datosprobar$z   <- (datosprobar$z - znormp) / unormp
}
prueba <- datosprobar[,-c(5,6,7,8,9,10)]
setDT(prueba)
Eprueba <- prueba$E
coordenadasprueba <- model.matrix(~.+0, data = prueba[,-c("E"), with=F])
Eprueba <- as.numeric(Eprueba)
dprueba <- xgb.DMatrix(data = coordenadasprueba, label = Eprueba)
pred <- predict (bst, dprueba)

write.table(pred, file = paste(sep="",casob,"_PRED.csv"), row.names=FALSE, col.names=FALSE, sep=";", dec=",")
predex <- cbind(datosprobar[,c(6,7,8,9,10,1,2,3,5)], pred)
write.table(predex, file = paste(sep="",casob,"_galaxias_PRED_PROB.csv"), row.names=FALSE, col.names=FALSE, sep=";", dec=",")
predm <- ifelse (pred > 0.5, 2, 1)
predex <- cbind(datosprobar[,c(6,7,8,9,10,1,2,3,5)], predm)
write.table(predex, file = paste(sep="",casob,"_galaxias_PRED_CLASE.csv"), row.names=FALSE, col.names=FALSE, sep=";", dec=",")
predexf <- predex[predex$pred < 2,]
write.table(predexf, file = paste(sep="",casob,"_galaxias_PRED_FALLOS.csv"), row.names=FALSE, col.names=FALSE, sep=";", dec=",")

pred <- ifelse (pred > 0.5, 1, 0)
mc <- confusionMatrix(factor(pred,levels=0:1), factor(Eprueba,levels=0:1))
lineag1 <- paste(sep="",
  "Hyperleda",
  sprintf(" nr:%6d",nr),
  sprintf(" eta:%6.3f",eta),
  sprintf(" cs:%5.2f",cs),
  sprintf(" md:%5.1f",md),
  sprintf(" Kappa:%6.3f",unname(mc$overall[2])),
  sprintf(" McnemarPValue:%6.3f",unname(mc$overall[7])),
  sprintf(" Acierto:%7.3f",unname(mc$overall[1])*100),"%",
  sprintf(" Entrena:%6d No:%6d Si:%6d ( %7.3f%s )",nrow(entrena),nrow(entrena)-sum(entrena$E),sum(entrena$E),sum(entrena$E)*100/nrow(entrena),"%"),
  sprintf(" Prueba: %6d No:%6d Si:%6d ( %7.3f%s )",nrow(prueba),nrow(prueba)-sum(prueba$E),sum(prueba$E),sum(prueba$E)*100/nrow(prueba),"%"),
  sprintf(" %6.0f",mc$table[1,1]),
  sprintf(" %6.0f",mc$table[2,2])
)
lineag1c <- gsub('\\.',',',paste(sep="",
  "Hyperleda",
  sprintf(" %5.3f",eta),
  sprintf(" %5.2f",cs),
  sprintf(" %4.1f",md),
  sprintf(" %8.5f",unname(mc$overall[2])),
  sprintf(" %7.5f",unname(mc$overall[7])),
  sprintf(" %7.5f",unname(mc$overall[1])),
  sprintf(" %6d %6d %6d %7.5f",nrow(entrena),nrow(entrena)-sum(entrena$E),sum(entrena$E),sum(entrena$E)/nrow(entrena)),
  sprintf(" %6d %6d %6d %7.5f",nrow(prueba),nrow(prueba)-sum(prueba$E),sum(prueba$E),sum(prueba$E)/nrow(prueba)),
  sprintf(" %6d",mc$table[1,1]),
  sprintf(" %6d",mc$table[2,2])
))
ponderae <- ponderae_az
sumapon  <- sumapon_az
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
    dif <- as.numeric(labels) - as.numeric(ppred)
    err <- sum(dif * dif * ponderae)/sumapon
    return(list(metric = "error", value = err))
  }
  param <- list(
    objective    = logregobj
    ,eval_metric = evalerror
  )
}
pred <- predict (bstaz, dprueba)
pred <- ifelse (pred > 0.5, 1, 0)
mc <- confusionMatrix(factor(pred,levels=0:1), factor(Eprueba,levels=0:1))
lineag2 <- paste(sep="",
  "Azar     ",
  sprintf(" nr:%6d",nr),
  sprintf(" eta:%6.3f",eta),
  sprintf(" cs:%5.2f",cs),
  sprintf(" md:%5.1f",md),
  sprintf(" Kappa:%6.3f",unname(mc$overall[2])),
  sprintf(" McnemarPValue:%6.3f",unname(mc$overall[7])),
  sprintf(" Acierto:%7.3f",unname(mc$overall[1])*100),"%",
  sprintf(" Entrena:%6d No:%6d Si:%6d ( %7.3f%s )",nrow(entrena),nrow(entrena)-sum(entrena$E),sum(entrena$E),sum(entrena$E)*100/nrow(entrena),"%"),
  sprintf(" Prueba: %6d No:%6d Si:%6d ( %7.3f%s )",nrow(prueba),nrow(prueba)-sum(prueba$E),sum(prueba$E),sum(prueba$E)*100/nrow(prueba),"%"),
  sprintf(" %6.0f",mc$table[1,1]),
  sprintf(" %6.0f",mc$table[2,2])
)
lineag2c <- gsub('\\.',',',paste(sep="",
  "Azar     ",
  sprintf(" %5.3f",eta),
  sprintf(" %5.2f",cs),
  sprintf(" %4.1f",md),
  sprintf(" %8.5f",unname(mc$overall[2])),
  sprintf(" %7.5f",unname(mc$overall[7])),
  sprintf(" %7.5f",unname(mc$overall[1])),
  sprintf(" %6d %6d %6d %7.5f",nrow(entrena),nrow(entrena)-sum(entrena$E),sum(entrena$E),sum(entrena$E)/nrow(entrena)),
  sprintf(" %6d %6d %6d %7.5f",nrow(prueba),nrow(prueba)-sum(prueba$E),sum(prueba$E),sum(prueba$E)/nrow(prueba)),
  sprintf(" %6d",mc$table[1,1]),
  sprintf(" %6d",mc$table[2,2])
))
ponderae <- ponderae_azaj
sumapon  <- sumapon_azaj
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
    dif <- as.numeric(labels) - as.numeric(ppred)
    err <- sum(dif * dif * ponderae)/sumapon
    return(list(metric = "error", value = err))
  }
  param <- list(
    objective    = logregobj
    ,eval_metric = evalerror
  )
}
pred <- predict (bstazaj, dprueba)
pred <- ifelse (pred > 0.5, 1, 0)
mc <- confusionMatrix(factor(pred,levels=0:1), factor(Eprueba,levels=0:1))
lineag3 <- paste(sep="",
  "Azar ajus",
  sprintf(" nr:%6d",nr),
  sprintf(" eta:%6.3f",eta),
  sprintf(" cs:%5.2f",cs),
  sprintf(" md:%5.1f",md),
  sprintf(" Kappa:%6.3f",unname(mc$overall[2])),
  sprintf(" McnemarPValue:%6.3f",unname(mc$overall[7])),
  sprintf(" Acierto:%7.3f",unname(mc$overall[1])*100),"%",
  sprintf(" Entrena:%6d No:%6d Si:%6d ( %7.3f%s )",nrow(entrena),nrow(entrena)-sum(entrena$E),sum(entrena$E),sum(entrena$E)*100/nrow(entrena),"%"),
  sprintf(" Prueba: %6d No:%6d Si:%6d ( %7.3f%s )",nrow(prueba),nrow(prueba)-sum(prueba$E),sum(prueba$E),sum(prueba$E)*100/nrow(prueba),"%"),
  sprintf(" %6.0f",mc$table[1,1]),
  sprintf(" %6.0f",mc$table[2,2])
)
lineag3c <- gsub('\\.',',',paste(sep="",
  "Azar ajus",
  sprintf(" %5.3f",eta),
  sprintf(" %5.2f",cs),
  sprintf(" %4.1f",md),
  sprintf(" %8.5f",unname(mc$overall[2])),
  sprintf(" %7.5f",unname(mc$overall[7])),
  sprintf(" %7.5f",unname(mc$overall[1])),
  sprintf(" %6d %6d %6d %7.5f",nrow(entrena),nrow(entrena)-sum(entrena$E),sum(entrena$E),sum(entrena$E)/nrow(entrena)),
  sprintf(" %6d %6d %6d %7.5f",nrow(prueba),nrow(prueba)-sum(prueba$E),sum(prueba$E),sum(prueba$E)/nrow(prueba)),
  sprintf(" %6d",mc$table[1,1]),
  sprintf(" %6d",mc$table[2,2])
))

if (casoa != casob){
  
  # Galaxias sector a
  
  ponderae <- ponderae_n
  sumapon  <- sumapon_n
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
      dif <- as.numeric(labels) - as.numeric(ppred)
      err <- sum(dif * dif * ponderae)/sumapon
      return(list(metric = "error", value = err))
    }
    param <- list(
      objective    = logregobj
      ,eval_metric = evalerror
    )
  }
  datosprobar  <- read.table(paste(sep="",casoa,"_galaxias_intA.csv"), header=TRUE, sep=";", na.strings="NA", dec=",", strip.white=TRUE)
  if (normalizar_estandar){
    datosprobar = normaliza(datosprobar, rangos)
  }
  if (normalizar_gly){
    datosprobar$x   <- (datosprobar$x - xnormp) / unormp
    datosprobar$y   <- (datosprobar$y - ynormp) / unormp
    datosprobar$z   <- (datosprobar$z - znormp) / unormp
  }
  prueba <- datosprobar[,-c(5,6,7,8,9,10)]
  setDT(prueba)
  Eprueba <- prueba$E
  coordenadasprueba <- model.matrix(~.+0, data = prueba[,-c("E"), with=F])
  Eprueba <- as.numeric(Eprueba)
  dprueba <- xgb.DMatrix(data = coordenadasprueba, label = Eprueba)
  pred <- predict (bst, dprueba)
  
  write.table(pred, file = paste(sep="",casoa,"_PRED.csv"), row.names=FALSE, col.names=FALSE, sep=";", dec=",")
  predex <- cbind(datosprobar[,c(6,7,8,9,10,1,2,3,5)], pred)
  write.table(predex, file = paste(sep="",casoa,"_galaxias_PRED_PROB.csv"), row.names=FALSE, col.names=FALSE, sep=";", dec=",")
  predm <- ifelse (pred > 0.5, 2, 1)
  predex <- cbind(datosprobar[,c(6,7,8,9,10,1,2,3,5)], predm)
  write.table(predex, file = paste(sep="",casoa,"_galaxias_PRED_CLASE.csv"), row.names=FALSE, col.names=FALSE, sep=";", dec=",")
  predexf <- predex[predex$pred < 2,]
  write.table(predexf, file = paste(sep="",casoa,"_galaxias_PRED_FALLOS.csv"), row.names=FALSE, col.names=FALSE, sep=";", dec=",")

  pred <- ifelse (pred > 0.5, 1, 0)
  mc <- confusionMatrix(factor(pred,levels=0:1), factor(Eprueba,levels=0:1))
  lineagp1 <- paste(sep="",
    "Hyperleda",
    sprintf(" nr:%6d",nr),
    sprintf(" eta:%6.3f",eta),
    sprintf(" cs:%5.2f",cs),
    sprintf(" md:%5.1f",md),
    sprintf(" Kappa:%6.3f",unname(mc$overall[2])),
    sprintf(" McnemarPValue:%6.3f",unname(mc$overall[7])),
    sprintf(" Acierto:%7.3f",unname(mc$overall[1])*100),"%",
    sprintf(" Entrena:%6d No:%6d Si:%6d ( %7.3f%s )",nrow(entrena),nrow(entrena)-sum(entrena$E),sum(entrena$E),sum(entrena$E)*100/nrow(entrena),"%"),
    sprintf(" Prueba: %6d No:%6d Si:%6d ( %7.3f%s )",nrow(prueba),nrow(prueba)-sum(prueba$E),sum(prueba$E),sum(prueba$E)*100/nrow(prueba),"%"),
    sprintf(" %6.0f",mc$table[1,1]),
    sprintf(" %6.0f",mc$table[2,2])
  )
  lineagp1c <- gsub('\\.',',',paste(sep="",
    "Hyperleda",
    sprintf(" %5.3f",eta),
    sprintf(" %5.2f",cs),
    sprintf(" %4.1f",md),
    sprintf(" %8.5f",unname(mc$overall[2])),
    sprintf(" %7.5f",unname(mc$overall[7])),
    sprintf(" %7.5f",unname(mc$overall[1])),
    sprintf(" %6d %6d %6d %7.5f",nrow(entrena),nrow(entrena)-sum(entrena$E),sum(entrena$E),sum(entrena$E)/nrow(entrena)),
    sprintf(" %6d %6d %6d %7.5f",nrow(prueba),nrow(prueba)-sum(prueba$E),sum(prueba$E),sum(prueba$E)/nrow(prueba)),
    sprintf(" %6d",mc$table[1,1]),
    sprintf(" %6d",mc$table[2,2])
  ))
  ponderae <- ponderae_az
  sumapon  <- sumapon_az
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
      dif <- as.numeric(labels) - as.numeric(ppred)
      err <- sum(dif * dif * ponderae)/sumapon
      return(list(metric = "error", value = err))
    }
    param <- list(
      objective    = logregobj
      ,eval_metric = evalerror
    )
  }
  pred <- predict (bstaz, dprueba)
  pred <- ifelse (pred > 0.5, 1, 0)
  mc <- confusionMatrix(factor(pred,levels=0:1), factor(Eprueba,levels=0:1))
  lineagp2 <- paste(sep="",
    "Azar     ",
    sprintf(" nr:%6d",nr),
    sprintf(" eta:%6.3f",eta),
    sprintf(" cs:%5.2f",cs),
    sprintf(" md:%5.1f",md),
    sprintf(" Kappa:%6.3f",unname(mc$overall[2])),
    sprintf(" McnemarPValue:%6.3f",unname(mc$overall[7])),
    sprintf(" Acierto:%7.3f",unname(mc$overall[1])*100),"%",
    sprintf(" Entrena:%6d No:%6d Si:%6d ( %7.3f%s )",nrow(entrena),nrow(entrena)-sum(entrena$E),sum(entrena$E),sum(entrena$E)*100/nrow(entrena),"%"),
    sprintf(" Prueba: %6d No:%6d Si:%6d ( %7.3f%s )",nrow(prueba),nrow(prueba)-sum(prueba$E),sum(prueba$E),sum(prueba$E)*100/nrow(prueba),"%"),
    sprintf(" %6.0f",mc$table[1,1]),
    sprintf(" %6.0f",mc$table[2,2])
  )
  lineagp2c <- gsub('\\.',',',paste(sep="",
    "Azar     ",
    sprintf(" %5.3f",eta),
    sprintf(" %5.2f",cs),
    sprintf(" %4.1f",md),
    sprintf(" %8.5f",unname(mc$overall[2])),
    sprintf(" %7.5f",unname(mc$overall[7])),
    sprintf(" %7.5f",unname(mc$overall[1])),
    sprintf(" %6d %6d %6d %7.5f",nrow(entrena),nrow(entrena)-sum(entrena$E),sum(entrena$E),sum(entrena$E)/nrow(entrena)),
    sprintf(" %6d %6d %6d %7.5f",nrow(prueba),nrow(prueba)-sum(prueba$E),sum(prueba$E),sum(prueba$E)/nrow(prueba)),
    sprintf(" %6d",mc$table[1,1]),
    sprintf(" %6d",mc$table[2,2])
  ))
  ponderae <- ponderae_azaj
  sumapon  <- sumapon_azaj
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
      dif <- as.numeric(labels) - as.numeric(ppred)
      err <- sum(dif * dif * ponderae)/sumapon
      return(list(metric = "error", value = err))
    }
    param <- list(
      objective    = logregobj
      ,eval_metric = evalerror
    )
  }
  pred <- predict (bstazaj, dprueba)
  pred <- ifelse (pred > 0.5, 1, 0)
  mc <- confusionMatrix(factor(pred,levels=0:1), factor(Eprueba,levels=0:1))
  lineagp3 <- paste(sep="",
    "Azar ajus",
    sprintf(" nr:%6d",nr),
    sprintf(" eta:%6.3f",eta),
    sprintf(" cs:%5.2f",cs),
    sprintf(" md:%5.1f",md),
    sprintf(" Kappa:%6.3f",unname(mc$overall[2])),
    sprintf(" McnemarPValue:%6.3f",unname(mc$overall[7])),
    sprintf(" Acierto:%7.3f",unname(mc$overall[1])*100),"%",
    sprintf(" Entrena:%6d No:%6d Si:%6d ( %7.3f%s )",nrow(entrena),nrow(entrena)-sum(entrena$E),sum(entrena$E),sum(entrena$E)*100/nrow(entrena),"%"),
    sprintf(" Prueba: %6d No:%6d Si:%6d ( %7.3f%s )",nrow(prueba),nrow(prueba)-sum(prueba$E),sum(prueba$E),sum(prueba$E)*100/nrow(prueba),"%"),
    sprintf(" %6.0f",mc$table[1,1]),
    sprintf(" %6.0f",mc$table[2,2])
  )
  lineagp3c <- gsub('\\.',',',paste(sep="",
    "Azar ajus",
    sprintf(" %5.3f",eta),
    sprintf(" %5.2f",cs),
    sprintf(" %4.1f",md),
    sprintf(" %8.5f",unname(mc$overall[2])),
    sprintf(" %7.5f",unname(mc$overall[7])),
    sprintf(" %7.5f",unname(mc$overall[1])),
    sprintf(" %6d %6d %6d %7.5f",nrow(entrena),nrow(entrena)-sum(entrena$E),sum(entrena$E),sum(entrena$E)/nrow(entrena)),
    sprintf(" %6d %6d %6d %7.5f",nrow(prueba),nrow(prueba)-sum(prueba$E),sum(prueba$E),sum(prueba$E)/nrow(prueba)),
    sprintf(" %6d",mc$table[1,1]),
    sprintf(" %6d",mc$table[2,2])
  ))
}
{
  cat("Entrena:",casoa,"\n")
  cat("Prueba :",casob,"\n\n")
  cat(sprintf("Entrenar con a, probar con la muestra de entrenamiento\n%s\n%s\n\n",linea0,linea0c))
  cat(sprintf("Entrenar con a, probar con b\n%s\n%s\n%s\n\n",linea1,linea2,linea3))
  cat(sprintf("%s\n%s\n%s\n\n",linea1c,linea2c,linea3c))

  cat(sprintf("Entrenar con a HyperLeda, probar con b\n%s\n%s\n%s\n\n",linea1,linea2bis,linea3bis))
  cat(sprintf("%s\n%s\n%s\n\n",linea1c,linea2cbis,linea3cbis))
  
  cat(sprintf("Entrenar con a, probar galaxias b\n%s\n%s\n%s\n\n",lineag1,lineag2,lineag3))
  cat(sprintf("%s\n%s\n%s\n\n",lineag1c,lineag2c,lineag3c))
  cat(sprintf("Entrenar con a, probar muestra total b\n%s\n\n",lineat1))
  cat(sprintf("%s\n\n",lineat1c))
  if (casoa != casob){
    cat(sprintf("Galaxias sector a\n%s\n%s\n%s\n\n",lineagp1,lineagp2,lineagp3))
    cat(sprintf("%s\n%s\n%s\n\n",lineagp1c,lineagp2c,lineagp3c))
  }
}

