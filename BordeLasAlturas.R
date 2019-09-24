#Borde Las Alturas

setwd("C:/RAnalysis/Unmarked/distsamp/borde")

#Figura1. Densidad de C. crassidigitus en borde versus interior
library(ggplot2)

cant <- read.csv("cantidad.csv", sep = ",", header = T)
View(cant)

densi <- read.csv("Densidad_borde_interior.csv", sep = ",", header = T)
View(densi)
ggplot(densi, aes(x= B.i)) + geom_histogram()


bd <- read.csv("LRC_Ecotono1.csv", sep = "", header = T)
bd
head(bd)
class(bd)
summary(bd)
mean(bd$temp_transecto, na.rm = T)
#na.rm es para que no tenga en cuanta los NA cuando no hay datos, en el promedio pej
#la desv estandar indica que tanto los datos se alejan del promedio
summary(bd$temp_transecto)
boxplot(LRC ~ Ecotono, data = bd, xlab="Ecotone", ylab="LRC (mm)", notch=F)
boxplot(temp_transecto ~ transecto, data = bd, xlab="Transecto", ylab="Temperatura ambiente (??C)", 
        notch=F, outline=F)

ggplot(bd, aes(x=Ecotono, y=LRC, fill=Ecotono)) +
  geom_boxplot(alpha=0.4) +
  stat_summary(fun.y=mean, geom="point", shape=20, size=2, color="red", fill="red") +
  theme(legend.position="none") +
  scale_fill_brewer(palette="Set3")


#Distance sampling para C. crassidigitus en los transectos de borde de Las ALturas, Coton

library(unmarked)

library(AICcmodavg)

library(ggplot2)

dists <- read.csv("distance_crassidigitus.csv", sep=";", header=TRUE)
View(dists)
class(dists)
summary(dists$cod_transecto)
# Darle formato de factor a la columna que contiene las etiquetas de nuestras unidades de muestreo 
dists$cod_transecto <- as.factor(dists$cod_transecto)

levels(dists$cod_transecto) <- c(levels(dists$cod_transecto))

# se convierten los datos de frecuencias de distancias a rangos de distancias, en este caso cada 0.5 metros

cp = c(0, 1, 2, 3, 4, 5)

# organizamos el formato de nuestros datos con la funcion formatDistData

cdata <- formatDistData(dists, "dist_m", "cod_transecto", cp)
head(cdata)
class(cdata)
hist(cdata, xlab="distance (m)", main="", cex.lab=0.8, cex.axis=0.8)

# importamos los datos de las covariables desde el archivo "Atelopus_cov.csv", de manera que coincidan con la organizaci?n de los datos de los conteos en cada unidad de 
# muestreo por cada rango de distancia 

covs <- read.csv("covariables.csv", sep = "", header=TRUE)
head(covs)
class(covs)
scale_covs=scale(covs[ , 4:6])
head(scale_covs)
class(scale_covs)
scale_covs=as.data.frame(scale_covs)

library(fuzzySim) 

multicol(covs[ , 4:6])

###################Resultados###############
                Rsquared Tolerance      VIF
dist_de_borde 0.13781434 0.8621857 1.159843
temp          0.12658372 0.8734163 1.144929
Precipit      0.09973976 0.9002602 1.110790
############################################

# con la funcion unmarkedFrameDS organizamos nuestros datos para correrlos con la funci?n distamp

umf <- unmarkedFrameDS(y=as.matrix(cdata), siteCovs=scale_covs, survey="line", 
                       dist.breaks=c(0, 1, 2, 3, 4, 5), 
                       tlength = covs$longitud, unitsIn="m")
umf

head(umf)
# procedemos a ajustar nuestros datos a los modelos, iniciando con un modelo nulo
# con las diferentes funciones de distribuci?n halfnormal, hazard, exp y uniforme
#En distance el modelo nulo dice que la deteccion y la abundancia no varia, no depende de ninguna variable

hn_Null <- distsamp(~1~1, umf, keyfun="halfnorm", output="density", unitsOut="ha")


hz_Null <- distsamp(~1~1, umf, keyfun="hazard", output="density", unitsOut="ha")
##NOTA: [1] "Hessian is singular. Try using fewer covariates or providing starting values."

exp_Null <- distsamp(~1~1, umf, keyfun="exp", output="density", unitsOut="ha")

unf_Null <- distsamp(~1~1, umf, keyfun="uniform", output="density", unitsOut="ha")


# a continuaci?n probamos los modelos en los que la probabilidad de detecci?n sin covariables y 
#la densidad es explicada por si es borde o es interior
# con cada funci?n de distribuci?n
#El primero es la deteccion y el segundo es la abundancia

hn_Nulldb <- distsamp(~1~dist_de_borde, umf, keyfun="halfnorm", output="density", unitsOut="ha")

hz_Nulldb <- distsamp(~1~dist_de_borde, umf, keyfun="hazard", output="density", unitsOut="ha")

exp_Nulldb <- distsamp(~1~dist_de_borde, umf, keyfun="exp", output="density", unitsOut="ha")

unf_Nulldb <- distsamp(~1~dist_de_borde, umf, keyfun="uniform", output="density", unitsOut="ha")

# ajustamos los modelos en los que la probabilidad de detecci?n este explicada por borde o interior 
# y la densidad sin covariables con cada funci?n de distribuci?n

hn_dbNull <- distsamp(~dist_de_borde~1, umf, keyfun="halfnorm", output="density", unitsOut="ha")

hz_dbNull <- distsamp(~dist_de_borde~1, umf, keyfun="hazard", output="density", unitsOut="ha")

exp_dbNull <- distsamp(~dist_de_borde~1, umf, keyfun="exp", output="density", unitsOut="ha")

unf_dbNull <- distsamp(~dist_de_borde~1, umf, keyfun="uniform", output="density", unitsOut="ha")

# ajustamos los modelos en los que la probabilidad de detecci?n y 
# la densidad esten explicadas por borde o interior con cada funci?n de distribuci?n

hn_dbdb <- distsamp(~dist_de_borde~dist_de_borde, umf, keyfun="halfnorm", output="density", unitsOut="ha")

hz_dbdb <- distsamp(~dist_de_borde~dist_de_borde, umf, keyfun="hazard", output="density", unitsOut="ha")

exp_dbdb <- distsamp(~dist_de_borde~dist_de_borde, umf, keyfun="exp", output="density", unitsOut="ha")

unf_dbdb <- distsamp(~dist_de_borde~dist_de_borde, umf, keyfun="uniform", output="density", unitsOut="ha")


# ajustamos los modelos en los que la probabilidad de detecci?n este explicada por precipitacion 
# y la densidad sin covariables con cada funci?n de distribuci?n

hn_preNull <- distsamp(~Precipit~1, umf, keyfun="halfnorm", output="density", unitsOut="ha")

hz_preNull <- distsamp(~Precipit~1, umf, keyfun="hazard", output="density", unitsOut="ha")

exp_preNull <- distsamp(~Precipit~1, umf, keyfun="exp", output="density", unitsOut="ha")

unf_preNull <- distsamp(~Precipit~1, umf, keyfun="uniform", output="density", unitsOut="ha")

# ajustamos los modelos en los que la probabilidad de deteccion esta explicada por precipitacion y 
# la densidad esten explicadas por borde o interior con cada funci?n de distribuci?n

hn_predb <- distsamp(~Precipit~dist_de_borde, umf, keyfun="halfnorm", output="density", unitsOut="ha")

hz_predb <- distsamp(~Precipit~dist_de_borde, umf, keyfun="hazard", output="density", unitsOut="ha")

exp_predb <- distsamp(~Precipit~dist_de_borde, umf, keyfun="exp", output="density", unitsOut="ha")

unf_predb <- distsamp(~Precipit~dist_de_borde, umf, keyfun="uniform", output="density", unitsOut="ha")

# ajustamos los modelos en los que la probabilidad de detecci?n este explicada por temperatura 
# y la densidad sin covariables con cada funci?n de distribuci?n

hn_tempNull <- distsamp(~temp~1, umf, keyfun="halfnorm", output="density", unitsOut="ha")

hz_tempNull <- distsamp(~temp~1, umf, keyfun="hazard", output="density", unitsOut="ha")

exp_tempNull <- distsamp(~temp~1, umf, keyfun="exp", output="density", unitsOut="ha")

unf_tempNull <- distsamp(~temp~1, umf, keyfun="uniform", output="density", unitsOut="ha")

# ajustamos los modelos en los que la probabilidad de deteccion esta explicada por precipitacion y 
# la densidad esten explicadas por borde o interior con cada funci?n de distribuci?n

hn_tempdb <- distsamp(~temp~dist_de_borde, umf, keyfun="halfnorm", output="density", unitsOut="ha")

hz_tempdb <- distsamp(~temp~dist_de_borde, umf, keyfun="hazard", output="density", unitsOut="ha")

exp_tempdb <- distsamp(~temp~dist_de_borde, umf, keyfun="exp", output="density", unitsOut="ha")

unf_tempdb <- distsamp(~temp~dist_de_borde, umf, keyfun="uniform", output="density", unitsOut="ha")

#ahora ajustamos los modelos en los que la probabilidad de deteccion este explicada por
#PRECIPITACION + TEMP

hn_pretempNull <- distsamp(~Precipit+temp~1, umf, keyfun="halfnorm", output="density", unitsOut="ha")

hz_pretempNull <- distsamp(~Precipit+temp~1, umf, keyfun="hazard", output="density", unitsOut="ha")

exp_pretempNull <- distsamp(~Precipit+temp~1, umf, keyfun="exp", output="density", unitsOut="ha")

unf_pretempNull <- distsamp(~Precipit+temp~1, umf, keyfun="uniform", output="density", unitsOut="ha")


#ahora ajustamos los modelos en los que la probabilidad de deteccion este explicada por
#PRECIPITACION * TEMP

hn_preinttempNull <- distsamp(~Precipit*temp~1, umf, keyfun="halfnorm", output="density", unitsOut="ha")

hz_preinttempNull <- distsamp(~Precipit*temp~1, umf, keyfun="hazard", output="density", unitsOut="ha")

exp_preinttempNull <- distsamp(~Precipit*temp~1, umf, keyfun="exp", output="density", unitsOut="ha")

unf_preinttempNull <- distsamp(~Precipit*temp~1, umf, keyfun="uniform", output="density", unitsOut="ha")

#ahora ajustamos los modelos en los que la probabilidad de deteccion este explicada por
#PRECIPITACION + TEMP y la abundancia este explicada por distancia de borde

hn_pretempdb <- distsamp(~Precipit+temp~dist_de_borde, umf, keyfun="halfnorm", output="density", unitsOut="ha")

hz_pretempdb <- distsamp(~Precipit+temp~dist_de_borde, umf, keyfun="hazard", output="density", unitsOut="ha")

exp_pretempdb <- distsamp(~Precipit+temp~dist_de_borde, umf, keyfun="exp", output="density", unitsOut="ha")

unf_pretempdb <- distsamp(~Precipit+temp~dist_de_borde, umf, keyfun="uniform", output="density", unitsOut="ha")


#ahora ajustamos los modelos en los que la probabilidad de deteccion este explicada por
#PRECIPITACION * TEMP y la abundancia este explicada por distancia de borde

hn_preinttempdb <- distsamp(~Precipit*temp~dist_de_borde, umf, keyfun="halfnorm", output="density", unitsOut="ha")

hz_preinttempdb <- distsamp(~Precipit*temp~dist_de_borde, umf, keyfun="hazard", output="density", unitsOut="ha")

exp_preinttempdb <- distsamp(~Precipit*temp~dist_de_borde, umf, keyfun="exp", output="density", unitsOut="ha")

unf_preinttempdb <- distsamp(~Precipit*temp~dist_de_borde, umf, keyfun="uniform", output="density", unitsOut="ha")

#ahora ajustamos los modelos en los que la probabilidad de deteccion este explicada por
#TEMP + distancia al borde y la abundancia nulo

hn_tempdbNull <- distsamp(~temp+dist_de_borde~1, umf, keyfun="halfnorm", output="density", unitsOut="ha")

hz_tempdbNull <- distsamp(~temp+dist_de_borde~1, umf, keyfun="hazard", output="density", unitsOut="ha")

exp_tempdbNull <- distsamp(~temp+dist_de_borde~1, umf, keyfun="exp", output="density", unitsOut="ha")

unf_tempdbNull <- distsamp(~temp+dist_de_borde~1, umf, keyfun="uniform", output="density", unitsOut="ha")


#ahora ajustamos los modelos en los que la probabilidad de deteccion este explicada por
#TEMP * distancia al borde y la abundancia nulo

hn_tempintdbNull <- distsamp(~temp*dist_de_borde~1, umf, keyfun="halfnorm", output="density", unitsOut="ha")

hz_tempintdbNull <- distsamp(~temp*dist_de_borde~1, umf, keyfun="hazard", output="density", unitsOut="ha")

exp_tempintdbNull <- distsamp(~temp*dist_de_borde~1, umf, keyfun="exp", output="density", unitsOut="ha")

unf_tempintdbNull <- distsamp(~temp*dist_de_borde~1, umf, keyfun="uniform", output="density", unitsOut="ha")


#ahora ajustamos los modelos en los que la probabilidad de deteccion este explicada por
#TEMP + distancia al borde y la abundancia explicada por distancia al borde

hn_tempdbdb <- distsamp(~temp+dist_de_borde~dist_de_borde, umf, keyfun="halfnorm", output="density", unitsOut="ha")

hz_tempdbdb <- distsamp(~temp+dist_de_borde~dist_de_borde, umf, keyfun="hazard", output="density", unitsOut="ha")

exp_tempdbdb <- distsamp(~temp+dist_de_borde~dist_de_borde, umf, keyfun="exp", output="density", unitsOut="ha")

unf_tempdbdb <- distsamp(~temp+dist_de_borde~dist_de_borde, umf, keyfun="uniform", output="density", unitsOut="ha")

#ahora ajustamos los modelos en los que la probabilidad de deteccion este explicada por
#TEMP * distancia al borde y la abundancia explicada por distancia al borde

hn_tempintdbdb <- distsamp(~temp*dist_de_borde~dist_de_borde, umf, keyfun="halfnorm", output="density", unitsOut="ha")

hz_tempintdbdb <- distsamp(~temp*dist_de_borde~dist_de_borde, umf, keyfun="hazard", output="density", unitsOut="ha")

exp_tempintdbdb <- distsamp(~temp*dist_de_borde~dist_de_borde, umf, keyfun="exp", output="density", unitsOut="ha")

unf_tempintdbdb <- distsamp(~temp*dist_de_borde~dist_de_borde, umf, keyfun="uniform", output="density", unitsOut="ha")

# ajuste y selecci?n del modelo

#listar los modelos
Cand.models <- list(hn_Null, hz_Null, exp_Null, unf_Null, 
                    hn_Nulldb, hz_Nulldb, exp_Nulldb, unf_Nulldb, 
                    hn_dbNull, hz_dbNull, exp_dbNull, unf_dbNull, 
                    hn_dbdb, hz_dbdb, exp_dbdb, unf_dbdb, 
                    hn_preNull, hz_preNull, exp_preNull, unf_preNull, 
                    hn_predb, hz_predb, exp_predb, unf_predb, 
                    hn_tempNull, hz_tempNull, exp_tempNull, unf_tempNull,
                    hn_tempdb, hz_tempdb, exp_tempdb, unf_tempdb,
                    hn_pretempNull, hz_pretempNull, exp_pretempNull, unf_pretempNull,
                    hn_preinttempNull, hz_preinttempNull, exp_preinttempNull, unf_preinttempNull,
                    hn_pretempdb, hz_pretempdb, exp_pretempdb, unf_pretempdb,
                    hn_preinttempdb, hz_preinttempdb, exp_preinttempdb, unf_preinttempdb,
                    hn_tempdbNull, hz_tempdbNull, exp_tempdbNull, unf_tempdbNull,
                    hn_tempintdbNull, hz_tempintdbNull, exp_tempintdbNull, unf_tempintdbNull,
                    hn_tempdbdb, hz_tempdbdb, exp_tempdbdb, unf_tempdbdb,
                    hn_tempintdbdb, hz_tempintdbdb, exp_tempintdbdb, unf_tempintdbdb)

#nombrar los modelos
Modnames <- c("hn_Null", "hz_Null", "exp_Null", "unf_Null", 
              "hn_Nulldb", "hz_Nulldb", "exp_Nulldb", "unf_Nulldb", 
              "hn_dbNull", "hz_dbNull", "exp_dbNull", "unf_dbNull", 
              "hn_dbdb", "hz_dbdb", "exp_dbdb", "unf_dbdb", 
              "hn_preNull", "hz_preNull", "exp_preNull", "unf_preNull", 
              "hn_predb", "hz_predb", "exp_predb", "unf_predb", 
              "hn_tempNull", "hz_tempNull", "exp_tempNull", "unf_tempNull",
              "hn_tempdb", "hz_tempdb", "exp_tempdb", "unf_tempdb",
              "hn_pretempNull", "hz_pretempNull", "exp_pretempNull", "unf_pretempNull",
              "hn_preinttempNull", "hz_preinttempNull", "exp_preinttempNull", "unf_preinttempNull",
              "hn_pretempdb", "hz_pretempdb", "exp_pretempdb", "unf_pretempdb",
              "hn_preinttempdb", "hz_preinttempdb", "exp_preinttempdb", "unf_preinttempdb",
              "hn_tempdbNull", "hz_tempdbNull", "exp_tempdbNull", "unf_tempdbNull",
              "hn_tempintdbNull", "hz_tempintdbNull", "exp_tempintdbNull", "unf_tempintdbNull",
              "hn_tempdbdb", "hz_tempdbdb", "exp_tempdbdb", "unf_tempdbdb",
              "hn_tempintdbdb", "hz_tempintdbdb", "exp_tempintdbdb", "unf_tempintdbdb")

aictab(cand.set = Cand.models, modnames = Modnames,
       second.ord = TRUE, nobs = NULL, sort = TRUE)

#############################################################
#######################RESULTADOS############################
                   K   AICc Delta_AICc AICcWt Cum.Wt      LL
exp_tempNull       3 347.77       0.00   0.14   0.14 -170.43
hn_tempNull        3 348.56       0.78   0.10   0.24 -170.82
exp_pretempNull    4 349.05       1.27   0.08   0.32 -169.72
exp_tempdb         4 349.48       1.70   0.06   0.38 -169.94
hn_pretempNull     4 350.06       2.29   0.05   0.43 -170.23
hn_tempdb          4 350.09       2.31   0.05   0.47 -170.24
exp_tempdbNull     4 350.45       2.68   0.04   0.51 -170.43
exp_preNull        3 350.53       2.76   0.04   0.55 -171.80
hn_preNull         3 350.98       3.21   0.03   0.58 -172.03
exp_tempdbdb       5 351.00       3.22   0.03   0.60 -169.25
hn_tempdbNull      4 351.05       3.28   0.03   0.63 -170.73
exp_Nulldb         3 351.22       3.45   0.03   0.66 -172.15
exp_pretempdb      5 351.22       3.45   0.03   0.68 -169.36
hn_tempdbdb        5 351.34       3.57   0.02   0.71 -169.42
hz_tempNull        4 351.40       3.62   0.02   0.73 -170.90
exp_Null           2 351.48       3.71   0.02   0.75 -173.52
exp_predb          4 351.51       3.74   0.02   0.78 -170.96
exp_preinttempNull 5 351.53       3.76   0.02   0.80 -169.52
hn_Nulldb          3 351.80       4.03   0.02   0.82 -172.44
hn_predb           4 351.86       4.08   0.02   0.84 -171.13
hn_preinttempNull  5 351.96       4.18   0.02   0.86 -169.73
hn_pretempdb       5 351.99       4.22   0.02   0.87 -169.74
hn_Null            2 352.06       4.29   0.02   0.89 -173.81
hz_pretempNull     5 352.56       4.79   0.01   0.90 -170.03
exp_tempintdbNull  5 352.65       4.88   0.01   0.92 -170.08
exp_dbNull         3 352.68       4.91   0.01   0.93 -172.88
hz_tempdb          5 353.37       5.60   0.01   0.94 -170.44
exp_dbdb           4 353.79       6.01   0.01   0.94 -172.09
hn_tempintdbNull   5 353.86       6.09   0.01   0.95 -170.68
exp_tempintdbdb    6 353.94       6.17   0.01   0.96 -169.15
hn_dbNull          3 354.03       6.26   0.01   0.96 -173.56
exp_preinttempdb   6 354.05       6.27   0.01   0.97 -169.20
hn_dbdb            4 354.26       6.48   0.01   0.98 -172.33
hz_tempdbNull      5 354.29       6.52   0.01   0.98 -170.90
hn_preinttempdb    6 354.30       6.53   0.01   0.99 -169.32
hn_tempintdbdb     6 354.49       6.72   0.01   0.99 -169.42
hz_tempdbdb        6 354.68       6.91   0.00   1.00 -169.51
hz_pretempdb       6 355.18       7.40   0.00   1.00 -169.76
unf_tempintdbdb    2 373.11      25.33   0.00   1.00 -184.33
unf_tempdbdb       2 373.11      25.33   0.00   1.00 -184.33
unf_preinttempdb   2 373.11      25.33   0.00   1.00 -184.33
unf_pretempdb      2 373.11      25.33   0.00   1.00 -184.33
unf_tempdb         2 373.11      25.33   0.00   1.00 -184.33
unf_predb          2 373.11      25.33   0.00   1.00 -184.33
unf_dbdb           2 373.11      25.33   0.00   1.00 -184.33
unf_Nulldb         2 373.11      25.33   0.00   1.00 -184.33
unf_tempintdbNull  1 373.55      25.77   0.00   1.00 -185.70
unf_tempdbNull     1 373.55      25.77   0.00   1.00 -185.70
unf_preinttempNull 1 373.55      25.77   0.00   1.00 -185.70
unf_pretempNull    1 373.55      25.77   0.00   1.00 -185.70
unf_tempNull       1 373.55      25.77   0.00   1.00 -185.70
unf_preNull        1 373.55      25.77   0.00   1.00 -185.70
unf_dbNull         1 373.55      25.77   0.00   1.00 -185.70
unf_Null           1 373.55      25.77   0.00   1.00 -185.70
hz_preNull         4 374.71      26.94   0.00   1.00 -182.56
hz_predb           5 375.66      27.88   0.00   1.00 -181.58
hz_Nulldb          4 378.26      30.49   0.00   1.00 -184.33
hz_Null            3 378.33      30.55   0.00   1.00 -185.70
hz_preinttempNull  6 380.76      32.99   0.00   1.00 -182.56
hz_dbNull          4 381.00      33.23   0.00   1.00 -185.70
hz_dbdb            5 381.16      33.39   0.00   1.00 -184.33
hz_preinttempdb    7 382.59      34.81   0.00   1.00 -181.75
hz_tempintdbNull   6 387.05      39.28   0.00   1.00 -185.70
hz_tempintdbdb     7 387.76      39.98   0.00   1.00 -184.33
##############################################################

#exploramos los modelos con mejor ajuste, auqellos con deltaAICc menor a 2#

#asignamos nuevas covariables para la predicción utilizando las mismas covariables la 
#estructura modelo de distance nombrado "umf" utilizando siteCovs

nuevo<-siteCovs(umf)

class(nuevo)

#### Modelo exp_tempNull ####

densexp_tempNull= predict(exp_tempNull, type = "state", newdata = nuevo, 
                          appendData = T, level = 0.90)

densexp_tempNull

detexp_tempNull= predict(exp_tempNull, type = "det", newdata = nuevo, appendData = T, level = 0.90)

detexp_tempNull

qplot(temp, Predicted, data = detexp_tempNull, geom = "line",
      xlab = "Distancia borde", ylab = "deteccion") +
  geom_ribbon(aes(x = temp, ymin = lower, ymax = upper), alpha = 0.1) +
  theme_bw()

#### Modelo hn_tempNull ####

denshn_tempNull= predict(hn_tempNull, type = "state", newdata = nuevo, 
                          appendData = T, level = 0.90)

denshn_tempNull

dethn_tempNull= predict(hn_tempNull, type = "det", newdata = nuevo, 
                        appendData = T, level = 0.90)

dethn_tempNull

#### Modelo exp_pretempNull ####

densexp_pretempNull= predict(exp_pretempNull, type = "state", newdata = nuevo, 
                         appendData = T, level = 0.90)

densexp_pretempNull

detexp_pretempNull= predict(exp_pretempNull, type = "det", newdata = nuevo, 
                        appendData = T, level = 0.90)

detexp_pretempNull

#### Modelo exp_tempdb ####

densexp_tempdb= predict(exp_tempdb, type = "state", newdata = nuevo, 
                             appendData = T, level = 0.90)

densexp_tempdb

plot(densexp_temdb)

qplot(dist_de_borde, Predicted, data = densexp_tempdb, geom = "line",
      xlab = "Distancia borde", ylab = "Densidad") +
      geom_ribbon(aes(x = dist_de_borde, ymin = lower, ymax = upper), alpha = 0.1) +
      theme_bw()


detexp_tempdb= predict(exp_tempdb, type = "det", newdata = nuevo, 
                            appendData = T, level = 0.90)

detexp_tempdb

qplot(temp, Predicted, data = detexp_tempdb, geom = "line",
      xlab = "Distancia borde", ylab = "deteccion") +
  geom_ribbon(aes(x = temp, ymin = lower, ymax = upper), alpha = 0.1) +
  theme_bw()


modellistavg=fitList(exp_tempNull=exp_tempNull, hn_tempNull=hn_tempNull, 
                  exp_pretempNull=exp_pretempNull, exp_tempdb=exp_tempdb)


densidadavg <- predict(modellistavg, type = "state", newdata = nuevo, 
                       appendData = T, level = 0.90)
densidadavg

qplot(dist_de_borde, Predicted, data = densidadavg, geom = "line",
      xlab = "Distancia borde", ylab = "Densidad") +
  geom_ribbon(aes(x = dist_de_borde, ymin = lower, ymax = upper), alpha = 0.1) +
  theme_bw()

deteccionavg <- predict(modellistavg, type = "det",new=nuevo,appendData=T)
deteccionavg

qplot(temp, Predicted, data = deteccionavg, geom = "line",
      xlab = "Temperatura °C", ylab = "Deteccion") +
  geom_ribbon(aes(x = temp, ymin = lower, ymax = upper), alpha = 0.1) +
  theme_bw()

qplot(Precipit, Predicted, data = deteccionavg, geom = "line",
      xlab = "Temperatura °C", ylab = "Deteccion") +
  geom_ribbon(aes(x = Precipit, ymin = lower, ymax = upper), alpha = 0.1) +
  theme_bw()
