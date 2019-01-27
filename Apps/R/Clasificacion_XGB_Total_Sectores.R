rm(list = ls())
library(lattice)
library(ggplot2)
library(xgboost)
library(caret)
library(data.table)

desamb  <- "_Lucas"
casoini <-  1
casofin <- 40
disco   <- "C"

nr  <- 620
eta <- 0.084
md  <- 13
cs  <- 0.94

ponderar <- TRUE
salres   <- FALSE

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
    err <- as.numeric(sum(labels != (preds > 0)))/length(labels)
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

datos <- data.frame(
  celdas  =character(),
  galaxias=character(),
  stringsAsFactors=FALSE
)

datos[ 1,"celdas"] = paste(sep="",disco,":/Gly/G_1_10_448_8_16_15_70_550000_850000_muestra_tot_%s.csv")
datos[ 2,"celdas"] = paste(sep="",disco,":/Gly/G_1_10_456_8_16_15_70_500000_900000_muestra_tot_%s.csv")
datos[ 3,"celdas"] = paste(sep="",disco,":/Gly/G_1_10_458_8_16_15_70_450000_950000_muestra_tot_%s.csv")
datos[ 4,"celdas"] = paste(sep="",disco,":/Gly/G_1_10_458_8_16_15_70_400000_1000000_muestra_tot_%s.csv")
datos[ 5,"celdas"] = paste(sep="",disco,":/Gly/G_1_10_456_8_16_15_70_350000_1050000_muestra_tot_%s.csv")
datos[ 6,"celdas"] = paste(sep="",disco,":/Gly/G_1_10_450_8_16_15_70_300000_1100000_muestra_tot_%s.csv")
datos[ 7,"celdas"] = paste(sep="",disco,":/Gly/G_1_10_440_8_16_15_70_250000_1150000_muestra_tot_%s.csv")
datos[ 8,"celdas"] = paste(sep="",disco,":/Gly/G_1_10_432_8_16_15_70_200000_1200000_muestra_tot_%s.csv")

datos[ 9,"celdas"] = paste(sep="",disco,":/Gly/G_1_10_440_8_16_15_70_650000_950000_muestra_tot_%s.csv")
datos[10,"celdas"] = paste(sep="",disco,":/Gly/G_1_10_444_8_16_15_70_600000_1000000_muestra_tot_%s.csv")
datos[11,"celdas"] = paste(sep="",disco,":/Gly/G_1_10_444_8_16_15_70_550000_1050000_muestra_tot_%s.csv")
datos[12,"celdas"] = paste(sep="",disco,":/Gly/G_1_10_438_8_16_15_70_500000_1100000_muestra_tot_%s.csv")
datos[13,"celdas"] = paste(sep="",disco,":/Gly/G_1_10_434_8_16_15_70_450000_1150000_muestra_tot_%s.csv")
datos[14,"celdas"] = paste(sep="",disco,":/Gly/G_1_10_426_8_16_15_70_400000_1200000_muestra_tot_%s.csv")
datos[15,"celdas"] = paste(sep="",disco,":/Gly/G_1_10_418_8_16_15_70_350000_1250000_muestra_tot_%s.csv")
datos[16,"celdas"] = paste(sep="",disco,":/Gly/G_1_10_410_8_16_15_70_300000_1300000_muestra_tot_%s.csv")

datos[17,"celdas"] = paste(sep="",disco,":/Gly/G_1_10_434_8_16_15_70_750000_1050000_muestra_tot_%s.csv")
datos[18,"celdas"] = paste(sep="",disco,":/Gly/G_1_10_432_8_16_15_70_700000_1100000_muestra_tot_%s.csv")
datos[19,"celdas"] = paste(sep="",disco,":/Gly/G_1_10_424_8_16_15_70_650000_1150000_muestra_tot_%s.csv")
datos[20,"celdas"] = paste(sep="",disco,":/Gly/G_1_10_418_8_16_15_70_600000_1200000_muestra_tot_%s.csv")
datos[21,"celdas"] = paste(sep="",disco,":/Gly/G_1_10_410_8_16_15_70_550000_1250000_muestra_tot_%s.csv")
datos[22,"celdas"] = paste(sep="",disco,":/Gly/G_1_10_402_8_16_15_70_500000_1300000_muestra_tot_%s.csv")
datos[23,"celdas"] = paste(sep="",disco,":/Gly/G_1_10_394_8_16_15_70_450000_1350000_muestra_tot_%s.csv")
datos[24,"celdas"] = paste(sep="",disco,":/Gly/G_1_10_388_8_16_15_70_400000_1400000_muestra_tot_%s.csv")

datos[25,"celdas"] = paste(sep="",disco,":/Gly/G_1_10_420_8_16_15_70_850000_1150000_muestra_tot_%s.csv")
datos[26,"celdas"] = paste(sep="",disco,":/Gly/G_1_10_410_8_16_15_70_800000_1200000_muestra_tot_%s.csv")
datos[27,"celdas"] = paste(sep="",disco,":/Gly/G_1_10_400_8_16_15_70_750000_1250000_muestra_tot_%s.csv")
datos[28,"celdas"] = paste(sep="",disco,":/Gly/G_1_10_392_8_16_15_70_700000_1300000_muestra_tot_%s.csv")
datos[29,"celdas"] = paste(sep="",disco,":/Gly/G_1_10_386_8_16_15_70_650000_1350000_muestra_tot_%s.csv")
datos[30,"celdas"] = paste(sep="",disco,":/Gly/G_1_10_380_8_16_15_70_600000_1400000_muestra_tot_%s.csv")
datos[31,"celdas"] = paste(sep="",disco,":/Gly/G_1_10_376_8_16_15_70_550000_1450000_muestra_tot_%s.csv")
datos[32,"celdas"] = paste(sep="",disco,":/Gly/G_1_10_374_8_16_15_70_500000_1500000_muestra_tot_%s.csv")

datos[33,"celdas"] = paste(sep="",disco,":/Gly/G_1_10_386_8_16_15_70_950000_1250000_muestra_tot_%s.csv")
datos[34,"celdas"] = paste(sep="",disco,":/Gly/G_1_10_380_8_16_15_70_900000_1300000_muestra_tot_%s.csv")
datos[35,"celdas"] = paste(sep="",disco,":/Gly/G_1_10_376_8_16_15_70_850000_1350000_muestra_tot_%s.csv")
datos[36,"celdas"] = paste(sep="",disco,":/Gly/G_1_10_372_8_16_15_70_800000_1400000_muestra_tot_%s.csv")
datos[37,"celdas"] = paste(sep="",disco,":/Gly/G_1_10_368_8_16_15_70_750000_1450000_muestra_tot_%s.csv")
datos[38,"celdas"] = paste(sep="",disco,":/Gly/G_1_10_366_8_16_15_70_700000_1500000_muestra_tot_%s.csv")
datos[39,"celdas"] = paste(sep="",disco,":/Gly/G_1_10_362_8_16_15_70_650000_1550000_muestra_tot_%s.csv")
datos[40,"celdas"] = paste(sep="",disco,":/Gly/G_1_10_360_8_16_15_70_600000_1600000_muestra_tot_%s.csv")

datos[ 1,"galaxias"] = paste(sep="",disco,":/Gly/G_1_10_448_8_16_15_70_550000_850000_galaxias_intA.csv")
datos[ 2,"galaxias"] = paste(sep="",disco,":/Gly/G_1_10_456_8_16_15_70_500000_900000_galaxias_intA.csv")
datos[ 3,"galaxias"] = paste(sep="",disco,":/Gly/G_1_10_458_8_16_15_70_450000_950000_galaxias_intA.csv")
datos[ 4,"galaxias"] = paste(sep="",disco,":/Gly/G_1_10_458_8_16_15_70_400000_1000000_galaxias_intA.csv")
datos[ 5,"galaxias"] = paste(sep="",disco,":/Gly/G_1_10_456_8_16_15_70_350000_1050000_galaxias_intA.csv")
datos[ 6,"galaxias"] = paste(sep="",disco,":/Gly/G_1_10_450_8_16_15_70_300000_1100000_galaxias_intA.csv")
datos[ 7,"galaxias"] = paste(sep="",disco,":/Gly/G_1_10_440_8_16_15_70_250000_1150000_galaxias_intA.csv")
datos[ 8,"galaxias"] = paste(sep="",disco,":/Gly/G_1_10_432_8_16_15_70_200000_1200000_galaxias_intA.csv")

datos[ 9,"galaxias"] = paste(sep="",disco,":/Gly/G_1_10_440_8_16_15_70_650000_950000_galaxias_intA.csv")
datos[10,"galaxias"] = paste(sep="",disco,":/Gly/G_1_10_444_8_16_15_70_600000_1000000_galaxias_intA.csv")
datos[11,"galaxias"] = paste(sep="",disco,":/Gly/G_1_10_444_8_16_15_70_550000_1050000_galaxias_intA.csv")
datos[12,"galaxias"] = paste(sep="",disco,":/Gly/G_1_10_438_8_16_15_70_500000_1100000_galaxias_intA.csv")
datos[13,"galaxias"] = paste(sep="",disco,":/Gly/G_1_10_434_8_16_15_70_450000_1150000_galaxias_intA.csv")
datos[14,"galaxias"] = paste(sep="",disco,":/Gly/G_1_10_426_8_16_15_70_400000_1200000_galaxias_intA.csv")
datos[15,"galaxias"] = paste(sep="",disco,":/Gly/G_1_10_418_8_16_15_70_350000_1250000_galaxias_intA.csv")
datos[16,"galaxias"] = paste(sep="",disco,":/Gly/G_1_10_410_8_16_15_70_300000_1300000_galaxias_intA.csv")

datos[17,"galaxias"] = paste(sep="",disco,":/Gly/G_1_10_434_8_16_15_70_750000_1050000_galaxias_intA.csv")
datos[18,"galaxias"] = paste(sep="",disco,":/Gly/G_1_10_432_8_16_15_70_700000_1100000_galaxias_intA.csv")
datos[19,"galaxias"] = paste(sep="",disco,":/Gly/G_1_10_424_8_16_15_70_650000_1150000_galaxias_intA.csv")
datos[20,"galaxias"] = paste(sep="",disco,":/Gly/G_1_10_418_8_16_15_70_600000_1200000_galaxias_intA.csv")
datos[21,"galaxias"] = paste(sep="",disco,":/Gly/G_1_10_410_8_16_15_70_550000_1250000_galaxias_intA.csv")
datos[22,"galaxias"] = paste(sep="",disco,":/Gly/G_1_10_402_8_16_15_70_500000_1300000_galaxias_intA.csv")
datos[23,"galaxias"] = paste(sep="",disco,":/Gly/G_1_10_394_8_16_15_70_450000_1350000_galaxias_intA.csv")
datos[24,"galaxias"] = paste(sep="",disco,":/Gly/G_1_10_388_8_16_15_70_400000_1400000_galaxias_intA.csv")

datos[25,"galaxias"] = paste(sep="",disco,":/Gly/G_1_10_420_8_16_15_70_850000_1150000_galaxias_intA.csv")
datos[26,"galaxias"] = paste(sep="",disco,":/Gly/G_1_10_410_8_16_15_70_800000_1200000_galaxias_intA.csv")
datos[27,"galaxias"] = paste(sep="",disco,":/Gly/G_1_10_400_8_16_15_70_750000_1250000_galaxias_intA.csv")
datos[28,"galaxias"] = paste(sep="",disco,":/Gly/G_1_10_392_8_16_15_70_700000_1300000_galaxias_intA.csv")
datos[29,"galaxias"] = paste(sep="",disco,":/Gly/G_1_10_386_8_16_15_70_650000_1350000_galaxias_intA.csv")
datos[30,"galaxias"] = paste(sep="",disco,":/Gly/G_1_10_380_8_16_15_70_600000_1400000_galaxias_intA.csv")
datos[31,"galaxias"] = paste(sep="",disco,":/Gly/G_1_10_376_8_16_15_70_550000_1450000_galaxias_intA.csv")
datos[32,"galaxias"] = paste(sep="",disco,":/Gly/G_1_10_374_8_16_15_70_500000_1500000_galaxias_intA.csv")

datos[33,"galaxias"] = paste(sep="",disco,":/Gly/G_1_10_386_8_16_15_70_950000_1250000_galaxias_intA.csv")
datos[34,"galaxias"] = paste(sep="",disco,":/Gly/G_1_10_380_8_16_15_70_900000_1300000_galaxias_intA.csv")
datos[35,"galaxias"] = paste(sep="",disco,":/Gly/G_1_10_376_8_16_15_70_850000_1350000_galaxias_intA.csv")
datos[36,"galaxias"] = paste(sep="",disco,":/Gly/G_1_10_372_8_16_15_70_800000_1400000_galaxias_intA.csv")
datos[37,"galaxias"] = paste(sep="",disco,":/Gly/G_1_10_368_8_16_15_70_750000_1450000_galaxias_intA.csv")
datos[38,"galaxias"] = paste(sep="",disco,":/Gly/G_1_10_366_8_16_15_70_700000_1500000_galaxias_intA.csv")
datos[39,"galaxias"] = paste(sep="",disco,":/Gly/G_1_10_362_8_16_15_70_650000_1550000_galaxias_intA.csv")
datos[40,"galaxias"] = paste(sep="",disco,":/Gly/G_1_10_360_8_16_15_70_600000_1600000_galaxias_intA.csv")

if (salres){
  resultados <- data.frame(
    linea_cat   = character(),
    linea_az    = character(),
    linea_azaj  = character(),
    lineag_cat  = character(),
    lineag_az   = character(),
    lineag_azaj = character(),
    stringsAsFactors = FALSE
  )
}
fi  <- paste(sep="",disco,":/Gly/Total_Sectores_Celdas",desamb,".csv")
cat(sprintf("Ponderar   : %9s",ponderar),    file=fi, append=TRUE, sep="\n")
cat(sprintf("nr         : %9d",nr),          file=fi, append=TRUE, sep="\n")
cat(sprintf("eta        : %9.3f",eta),       file=fi, append=TRUE, sep="\n")
cat(sprintf("md         : %9d",md),          file=fi, append=TRUE, sep="\n")
cat(sprintf("cs         : %9.3f\n",cs),      file=fi, append=TRUE, sep="\n")
cat("---------------------------------- -------- ENTRENAR ---------- ----------- PROBAR -----------------------", file=fi, append=TRUE, sep="\n")
cat("Muestra    Kappa    PValue Acierto Obser. Vacias Ocupad %Ocupad Obser. Vacias Ocupad %Ocupad Vacias Ocupad", file=fi, append=TRUE, sep="\n")
cat("--------- -------- ------- ------- ------ ------ ------ ------- ------ ------ ------ ------- ------ ------", file=fi, append=TRUE, sep="\n")
fig <- paste(sep="",disco,":/Gly/Total_Sectores_Galaxias",desamb,".csv")
cat(sprintf("Ponderar   : %9s",ponderar),    file=fig, append=TRUE, sep="\n")
cat(sprintf("nr         : %9d",nr),          file=fig, append=TRUE, sep="\n")
cat(sprintf("eta        : %9.3f",eta),       file=fig, append=TRUE, sep="\n")
cat(sprintf("md         : %9d",md),          file=fig, append=TRUE, sep="\n")
cat(sprintf("cs         : %9.3f\n",cs),      file=fig, append=TRUE, sep="\n")
cat("---------------------------------- -------- ENTRENAR ---------- ----------- PROBAR -----------------------", file=fig, append=TRUE, sep="\n")
cat("Muestra    Kappa    PValue Acierto Obser. Vacias Ocupad %Ocupad Obser. Vacias Ocupad %Ocupad Vacias Ocupad", file=fig, append=TRUE, sep="\n")
cat("--------- -------- ------- ------- ------ ------ ------ ------- ------ ------ ------ ------- ------ ------", file=fig, append=TRUE, sep="\n")

op <- options(digits.secs=6)
momento <- Sys.time()
print(momento)
cat(sprintf("%s\n", momento), file=fi, append=TRUE, sep="\n")
options(op)
ptm <- proc.time()

casos   <- casofin - casoini + 1
pb      <- winProgressBar(title="Total sectores", min=0, max=casos, width=400)
casoact <- 0

for(i in casoini:casofin){
  casoact <- casoact + 1  

  # Celdas
  cat(sprintf(datos[i,"celdas"],""), file=fi, append=TRUE, sep="\n")
  datosentrenar <- read.table(sprintf(datos[i,"celdas"],"meA"), header=TRUE, sep=";", na.strings="NA", dec=",", strip.white=TRUE)
  ponderae      <- ifelse(datosentrenar$n == 0, 1, datosentrenar$n)
  entrena       <- datosentrenar[,c(1,2,3,4)]
  setDT(entrena)
  coordenadasEntrena  <- model.matrix(~.+0, data = entrena [,-c("E"), with=F]) 
  Eentrena      <- as.numeric(entrena$E)
  dentrena      <- xgb.DMatrix(data = coordenadasEntrena,  label = Eentrena) 
  rm(datosentrenar)
  rm(coordenadasEntrena)
  gc()
  datosprobar   <- read.table(sprintf(datos[i,"celdas"],"mpA"), header=TRUE, sep=";", na.strings="NA", dec=",", strip.white=TRUE)
  prueba        <- datosprobar  [,c(1,2,3,4)]
  setDT(prueba)
  coordenadasPrueba <- model.matrix(~.+0, data = prueba  [,-c("E"), with=F])
  Eprueba       <- as.numeric(prueba$E)
  dprueba       <- xgb.DMatrix(data = coordenadasPrueba,   label = Eprueba)
  rm(datosprobar)
  rm(coordenadasPrueba)
  gc()
  set.seed(1235) 
  bst           <- xgb.train(eta=eta, colsample_bylevel=cs, max_depth=md, data=dentrena, params=param, nrounds=nr, print_every_n=1000)
  pred          <- predict(bst, dprueba)
  pred          <- ifelse (pred > 0.5, 1, 0)
  mc            <- confusionMatrix(factor(pred,levels=0:1), factor(Eprueba,levels=0:1))
  if (salres){
    linea <- gsub('\\.',',',paste(sep="",
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
    resultados[i,"linea_cat"] = linea
  }
  cat(sprintf("Hyperleda %8.5f %7.5f %7.5f %6d %6d %6d %7.5f %6d %6d %6d %7.5f %6d %6d",unname(mc$overall[2]),unname(mc$overall[7]),unname(mc$overall[1]),nrow(entrena),nrow(entrena)-sum(entrena$E),sum(entrena$E),sum(entrena$E)/nrow(entrena),nrow(prueba),nrow(prueba)-sum(prueba$E),sum(prueba$E),sum(prueba$E)/nrow(prueba),mc$table[1,1],mc$table[2,2]), file=fi, append=TRUE, sep="\n")

  rm(pred)
  gc()
  datosentrenar <- read.table(sprintf(datos[i,"celdas"],"az_meA"), header=TRUE, sep=";", na.strings="NA", dec=",", strip.white=TRUE)
  ponderae      <- ifelse(datosentrenar$n == 0, 1, datosentrenar$n)
  entrena       <- datosentrenar[,c(1,2,3,4)]
  setDT(entrena)
  coordenadasEntrena  <- model.matrix(~.+0, data = entrena [,-c("E"), with=F]) 
  Eentrena      <- as.numeric(entrena$E)
  dentrena      <- xgb.DMatrix(data = coordenadasEntrena,  label = Eentrena) 
  rm(datosentrenar)
  rm(coordenadasEntrena)
  gc()
  datosprobar   <- read.table(sprintf(datos[i,"celdas"],"az_mpA"), header=TRUE, sep=";", na.strings="NA", dec=",", strip.white=TRUE)
  prueba        <- datosprobar  [,c(1,2,3,4)]
  setDT(prueba)
  coordenadasPrueba   <- model.matrix(~.+0, data = prueba  [,-c("E"), with=F])
  Eprueba       <- as.numeric(prueba$E)
  dprueba       <- xgb.DMatrix(data = coordenadasPrueba,   label = Eprueba)
  rm(datosprobar)
  rm(coordenadasPrueba)
  gc()
  set.seed(1235) 
  bstaz         <- xgb.train(eta = eta, colsample_bylevel = cs, max_depth = md, data = dentrena, params = param, nrounds = nr, print_every_n = 1000)
  pred          <- predict (bstaz, dprueba)
  pred          <- ifelse (pred > 0.5, 1, 0)
  mc            <- confusionMatrix(factor(pred,levels=0:1), factor(Eprueba,levels=0:1))
  if (salres){
    linea <- gsub('\\.',',',paste(sep="",
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
    resultados[i,"linea_az"] = linea
  }
  cat(sprintf("Azar      %8.5f %7.5f %7.5f %6d %6d %6d %7.5f %6d %6d %6d %7.5f %6d %6d",unname(mc$overall[2]),unname(mc$overall[7]),unname(mc$overall[1]),nrow(entrena),nrow(entrena)-sum(entrena$E),sum(entrena$E),sum(entrena$E)/nrow(entrena),nrow(prueba),nrow(prueba)-sum(prueba$E),sum(prueba$E),sum(prueba$E)/nrow(prueba),mc$table[1,1],mc$table[2,2]), file=fi, append=TRUE, sep="\n")
  rm(pred)
  gc()
  datosentrenar <- read.table(sprintf(datos[i,"celdas"],"azaj_meA"), header=TRUE, sep=";", na.strings="NA", dec=",", strip.white=TRUE)
  ponderae      <- ifelse(datosentrenar$n == 0, 1, datosentrenar$n)
  entrena       <- datosentrenar[,c(1,2,3,4)]
  setDT(entrena)
  coordenadasEntrena  <- model.matrix(~.+0, data = entrena [,-c("E"), with=F]) 
  Eentrena      <- as.numeric(entrena$E)
  dentrena      <- xgb.DMatrix(data = coordenadasEntrena,  label = Eentrena) 
  rm(datosentrenar)
  rm(coordenadasEntrena)
  gc()
  datosprobar   <- read.table(sprintf(datos[i,"celdas"],"azaj_mpA"), header=TRUE, sep=";", na.strings="NA", dec=",", strip.white=TRUE)
  prueba        <- datosprobar  [,c(1,2,3,4)]
  setDT(prueba)
  coordenadasPrueba   <- model.matrix(~.+0, data = prueba  [,-c("E"), with=F])
  Eprueba       <- as.numeric(prueba$E)
  dprueba       <- xgb.DMatrix(data = coordenadasPrueba,   label = Eprueba)
  rm(datosprobar)
  rm(coordenadasPrueba)
  gc()
  set.seed(1235) 
  bstazaj       <- xgb.train(eta = eta, colsample_bylevel = cs, max_depth = md, data = dentrena, params = param, nrounds = nr, print_every_n = 1000)
  pred          <- predict (bstazaj, dprueba)
  pred          <- ifelse (pred > 0.5, 1, 0)
  mc            <- confusionMatrix(factor(pred,levels=0:1), factor(Eprueba,levels=0:1))
  if (salres){
    linea <- gsub('\\.',',',paste(sep="",
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
    resultados[i,"linea_azaj"] = linea
  }
  cat(sprintf("Azar ajus %8.5f %7.5f %7.5f %6d %6d %6d %7.5f %6d %6d %6d %7.5f %6d %6d",unname(mc$overall[2]),unname(mc$overall[7]),unname(mc$overall[1]),nrow(entrena),nrow(entrena)-sum(entrena$E),sum(entrena$E),sum(entrena$E)/nrow(entrena),nrow(prueba),nrow(prueba)-sum(prueba$E),sum(prueba$E),sum(prueba$E)/nrow(prueba),mc$table[1,1],mc$table[2,2]), file=fi, append=TRUE, sep="\n")
  rm(pred)
  gc()

  # Galaxias
  
  cat(datos[i,"galaxias"], file=fig, append=TRUE, sep="\n")
  datosprobar <- read.table(datos[i,"galaxias"], header=TRUE, sep=";", na.strings="NA", dec=",", strip.white=TRUE)
  prueba      <- datosprobar  [,c(1,2,3,4)]
  setDT(prueba)
  coordenadasPrueba   <- model.matrix(~.+0, data = prueba  [,-c("E"), with=F])
  Eprueba     <- as.numeric(prueba$E)
  dprueba     <- xgb.DMatrix(data = coordenadasPrueba,   label = Eprueba)
  rm(datosprobar)
  rm(coordenadasPrueba)
  gc()
  pred        <- predict(bst, dprueba)
  pred        <- ifelse (pred > 0.5, 1, 0)
  mc          <- confusionMatrix(factor(pred,levels=0:1), factor(Eprueba,levels=0:1))
  if (salres){
    linea <- gsub('\\.',',',paste(sep="",
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
    resultados[i,"lineag_cat"] = linea
  }
  cat(sprintf("Hyperleda %8.5f %7.5f %7.5f %6d %6d %6d %7.5f %6d %6d %6d %7.5f %6d %6d",unname(mc$overall[2]),unname(mc$overall[7]),unname(mc$overall[1]),nrow(entrena),nrow(entrena)-sum(entrena$E),sum(entrena$E),sum(entrena$E)/nrow(entrena),nrow(prueba),nrow(prueba)-sum(prueba$E),sum(prueba$E),sum(prueba$E)/nrow(prueba),mc$table[1,1],mc$table[2,2]), file=fig, append=TRUE, sep="\n")
  rm(pred)
  pred        <- predict(bstaz, dprueba)
  pred        <- ifelse (pred > 0.5, 1, 0)
  mc          <- confusionMatrix(factor(pred,levels=0:1), factor(Eprueba,levels=0:1))
  if (salres){
    linea <- gsub('\\.',',',paste(sep="",
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
    resultados[i,"lineag_az"] = linea
  }
  cat(sprintf("Azar      %8.5f %7.5f %7.5f %6d %6d %6d %7.5f %6d %6d %6d %7.5f %6d %6d",unname(mc$overall[2]),unname(mc$overall[7]),unname(mc$overall[1]),nrow(entrena),nrow(entrena)-sum(entrena$E),sum(entrena$E),sum(entrena$E)/nrow(entrena),nrow(prueba),nrow(prueba)-sum(prueba$E),sum(prueba$E),sum(prueba$E)/nrow(prueba),mc$table[1,1],mc$table[2,2]), file=fig, append=TRUE, sep="\n")
  
  rm(pred)
  gc()
  pred    <- predict(bstazaj, dprueba)
  pred    <- ifelse (pred > 0.5, 1, 0)
  mc      <- confusionMatrix(factor(pred,levels=0:1), factor(Eprueba,levels=0:1))
  if (salres){
    linea <- gsub('\\.',',',paste(sep="",
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
    resultados[i,"lineag_azaj"] = linea
  }
  cat(sprintf("Azar ajus %8.5f %7.5f %7.5f %6d %6d %6d %7.5f %6d %6d %6d %7.5f %6d %6d",unname(mc$overall[2]),unname(mc$overall[7]),unname(mc$overall[1]),nrow(entrena),nrow(entrena)-sum(entrena$E),sum(entrena$E),sum(entrena$E)/nrow(entrena),nrow(prueba),nrow(prueba)-sum(prueba$E),sum(prueba$E),sum(prueba$E)/nrow(prueba),mc$table[1,1],mc$table[2,2]), file=fig, append=TRUE, sep="\n")

  rm(pred)
  gc()
  setWinProgressBar(pb, casoact, title=paste("Total sectores ", round(100*casoact/casos, 0), "% hecho"))
}
close(pb)

if (salres){
  for(i in casoini:casofin){
    l <- unlist(strsplit(datos[i,"celdas"], "[_]"))
    cat(sprintf("%s %s %s %s\n",l[4],l[9],l[10],resultados[i,"linea_cat"]))
    cat(sprintf("%s %s %s %s\n",l[4],l[9],l[10],resultados[i,"linea_az"]))
    cat(sprintf("%s %s %s %s\n\n",l[4],l[9],l[10],resultados[i,"linea_azaj"]))
  }
  for(i in casoini:casofin){
    l <- unlist(strsplit(datos[i,"celdas"], "[_]"))
    cat(sprintf("%s %s %s %s\n",l[4],l[9],l[10],resultados[i,"lineag_cat"]))
    cat(sprintf("%s %s %s %s\n",l[4],l[9],l[10],resultados[i,"lineag_az"]))
    cat(sprintf("%s %s %s %s\n\n",l[4],l[9],l[10],resultados[i,"lineag_azaj"]))
  }
}
tpo <- proc.time() - ptm
op  <- options(digits.secs=6)
momento <- Sys.time()
print(momento)
cat(sprintf("\n%s\n", momento), file=fi, append=TRUE, sep="\n")
options(op)
print(tpo)
cat(sprintf("Tiempo. user: %7.1f  system: %7.1f  elapsed: %7.1f\n",as.numeric(tpo[1]),as.numeric(tpo[2]),as.numeric(tpo[3])), file=fi, append=TRUE, sep="\n")
print("FIN")
