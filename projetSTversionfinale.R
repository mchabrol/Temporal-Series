rm(list = ls())



#install.packages("zoo")

#install.packages("tseries")

#install.packages("forecast")

#install.packages("ellipse")

#install.packages("ellipsis")

#install.packages("car")



library(zoo)

library(tseries)

library(fUnitRoots)



## importation des données


path <- "/Users/marionchabrol/Library/Mobile Documents/com~apple~CloudDocs/Documents/ENSAE/serie temp/"

setwd(path) #definit l'espace de travail (working directory ou "wd")
getwd() #affiche le wd
list.files() #liste les elements du wd

datafile <- "projet_charpente.csv" 
data <- read.csv(datafile,sep=";")


#mise en forme des dates

dates_char <- as.character(data$date)

dates_char[1] #date début

tail(dates_char,1) #date fin

dates <- as.yearmon(seq(from=1990+0/12, to=2023+1/12, by=1/12)) #définition des dates

valeur <- zoo(data$valeur, order.by=dates)


#on retire les 4 dernières valeurs (cela formera notre échantillon test)

T <- length(valeur)

test4 <- valeur[(T-3):T] #échantillon test

valeur <- valeur[1:(T-4)] #supprime les 4 dernieres valeurs

dates <- dates[1:(T-4)]


#visualisation de la série
monthplot(valeur)

plot(valeur)



#2

##différenciation

dvaleur <- diff(valeur,1)


#graphe des séries
par(mfrow=c(1,2))

plot(valeur)

plot(dvaleur)


#acf des séries
acf(ts(valeur))

acf(ts(dvaleur))


#pacf des séries
pacf(ts(valeur))

pacf(ts(dvaleur))





#Test de stationnarité ADF

Qtests <- function(series, k, fitdf=0) {
  
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    
    return(c("lag"=l,"pval"=pval))
    
  })
  
  return(t(pvals))
  
}


#tests ADF jusqu'à des residus non autocorrelés

adfTest_valid <- function(series, kmax, adftype){
  
  k <- 0
  
  noautocorr <- 0
  
  while (noautocorr==0){
    
    cat(paste0("ADF with ",k," lags: residuals OK? "))
    
    adf <- adfTest(series, lags=k, type=adftype)
    
    pvals <- Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))[,2]
    
    if (sum(pvals<0.05,na.rm=T)==0) { #il faut avoir des pvals > 0.05
      
      noautocorr <- 1; cat("OK \n")
      
    } else cat("nope \n")
    
    k <- k+1
    
  }
  
  return(adf)
  
}


#regression valeur sur les dates pour connaitre le type de test adf à faire

summary(lm(valeur ~ dates))

adf <- adfTest_valid(valeur,24,adftype="ct")

adf

#p_valeur > 0.05, on ne rejette pas l'hypothèse nulle de non stationnarité pour adf

kpss.test (valeur , null ="Trend")

#p_valeur < 0.05, on rejette l'hypothèse nulle de stationnarité pour kpss


#regression dvaleur (série différenciée) sur les dates pour connaitre le type de test adf à faire

summary(lm(dvaleur ~ dates[-1]))

#

adf <- adfTest_valid(dvaleur,24,"nc") #ok pour la série différenciée

adf

#p_valeur < 0.05, on rejette l'hypothèse nulle de non stationnarité

kpss.test(dvaleur , null ="Level")

#p_valeur > 0.05, on ne rejette pas l'hypothèse nulle de stationnarité



# Test de Breusch-Pagan (homoscedasticite)

lmtest :: bptest(lm(dvaleur ~ seq(1,length(dvaleur))))

# On ne rejette pas l ' homoscedasticite de X au niveau 5%



#3

#visualisation des données

plot(cbind(valeur,dvaleur))



###Partie 2

##4

#on a retenu dvaleur

par(mfrow=c(1,2))

#acf et pacf 
acf(ts(dvaleur), 24);pacf(ts(dvaleur), 24) #on regarde jusqu'à deux ans de retard


# estimation des ordres maximum pmax et qmax

pmax=3;qmax=2



signif <- function(estim){ #fonction de test des significations individuelles des coefficients
  
  coef <- estim$coef #recup coeff
  
  se <- sqrt(diag(estim$var.coef)) #calcul standard errors
  
  t <- coef/se #test de student
  
  pval <- (1-pnorm(abs(t)))*2
  
  return(rbind(coef,se,pval))
  
}



Qtests <- function(series, k, fitdf=0) { #series = ARMA, k=horizon temporel avec 24 (2 ans, données mensuelles)
  
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    
    return(c("lag"=l,"pval"=pval))
    
  })
  
  return(t(pvals))
  
}


#test de la fonction Qtest sur un premier modèle arma(3,2)
arma32 <- arima(dvaleur,c(3,0,2), include.mean=F)

Qtests(arma32$residuals, 24, 5) #5 = p+q = 3+2



arimafit <- function(estim){ #estim = modèle en tant que paramètres
  
  adjust <- round(signif(estim),3) #enlève les p_valeurs avec mille chiffres
  
  pvals <- Qtests(estim$residuals,24,length(estim$coef)-1)
  
  pvals <- matrix(apply(matrix(1:24,nrow=6),2,function(c) round(pvals[c,],3)),nrow=6)
  
  colnames(pvals) <- rep(c("lag", "pval"),4)
  
  cat("tests de nullite des coefficients :\n")
  
  print(adjust)
  
  cat("\n tests d'absence d'autocorrelation des residus : \n")
  
  print(pvals)
  
}


#on tests tous les modèles ARMA(p,q) avec p ≤ p_max et q ≤ q_max
estim <- arima(dvaleur,c(3,0,2), include.mean=F); arimafit(estim)

#pas bien ajusté, valide



estim <- arima(dvaleur,c(3,0,1), include.mean=F); arimafit(estim)

#pas bien ajusté, valide



estim <- arima(dvaleur,c(3,0,0), include.mean=F); arimafit(estim)

#bien ajusté, valide

ar3 <- arima(dvaleur,c(3,0,0), include.mean=F)



estim <- arima(dvaleur,c(2,0,1), include.mean=F); arimafit(estim)

#bien ajusté, valide

arma21 <- arima(dvaleur,c(2,0,1), include.mean=F)



estim <- arima(dvaleur,c(2,0,0), include.mean=F); arimafit(estim)

#bien ajusté, pas valide



estim <- arima(dvaleur,c(1,0,2), include.mean=F); arimafit(estim)

#pas bien ajusté, valide



estim <- arima(dvaleur,c(1,0,1), include.mean=F); arimafit(estim)

#bien ajusté, valide

arma11 <- arima(dvaleur,c(1,0,1), include.mean=F)



estim <- arima(dvaleur,c(1,0,0), include.mean=F); arimafit(estim)

#bien ajusté, pas valide



estim <- arima(dvaleur,c(0,0,2), include.mean=F); arimafit(estim)

#bien ajusté, valide

ma2 <- arima(dvaleur,c(0,0,2), include.mean=F)



estim <- arima(dvaleur,c(0,0,1), include.mean=F); arimafit(estim)

#bien ajusté, pas valide


#modèles sélectionnés

models <- c("ar3", "arma11", "arma21", "ma2"); names(models) <- models

apply(as.matrix(models),1, function(m) c("AIC"=AIC(get(m)), "BIC"=BIC(get(m))))



#meilleur modèle MA(2)



#5



#il faut regarder si les racines des modules du polynome sont supérieures à 1



ma2$coef

racine_ma2 <- polyroot(c(1, ma2$coef[1], ma2$coef[2]))



Mod(racine_ma2[1])

Mod(racine_ma2[2])



#Le modèle est bien causal (c'est un MA), on peut donc passer à l'ARIMA.

#Comme on a différencié une fois, d = 1. Le modèle est inversible.

#On obtient donc un ARIMA(0,1,2).



arima012 <- arima (valeur, c(0,1,2), include.mean=F)





# Partie 3 : Previsions

#6

#on vérifie que les résidus sont des bruits blancs gaussiens

tsdiag(arima012) 

qqnorm(arima012$residuals )

hist(arima012$residuals,40, freq = F)

mu<-mean(arima012$residuals)

sigma<-sd(arima012$residuals)

x<-seq(-40,40)

y<-dnorm(x,mu,sigma)

lines(x,y,lwd=0.5,col="blue")

residus <- resid(arima012)

#test de Ljung Box
Qtests(residus, 24, fitdf=2)

#Extraction des coefs du modele et de la variance des residus

arima012$coef

phi_1 <- as.numeric(arima012$coef[1])

phi_2 <- as.numeric(arima012$coef[2])

sigma2 <- as.numeric(arima012$sigma2)

phi_1

phi_2

sigma2



# Question 8

XT1 = predict(arima012, n.ahead=2)$pred[1]

XT2 = predict (arima012, n.ahead=2)$pred[2]

XT1

XT2



# On cherche d'abord a tracer le region de confiance univariee pour la serie originale a 95%.

install.packages("forecast")

library(forecast)

fore = forecast(arima012, h=4, level=95)

par(mfrow=c(1 ,1))

plot(fore,col=1,fcol=2,xlim =c(2020, 2025),ylim =c(25,200), shaded=TRUE,xlab="Temps",ylab="Valeur",
     
     main="Prevision pour la serie originale")

par(new=T)

plot(test4,col='#00FF00',fcol=2,xlim =c(2020, 2025),ylim =c(25,200), shaded=TRUE,xlab="Temps",ylab="Valeur",
     
     main="Prevision pour la serie originale")

#Ensuite , on represente la region de confiance bivariee a 95%.

library(ellipse)

arma = arima0(x, order=c(0,1,2))

Sigma <- matrix (c(sigma2,phi_1*sigma2,phi_1*sigma2,(phi_1)^2*sigma2 + sigma2),ncol =2)

inv_Sigma <- solve(Sigma)

plot ( XT1 , XT2 , xlim =c(0, 200) , ylim =c(50,200) , xlab =" Prevision de X(T+1) ",
       
       ylab ="Prevision de X(T+2)" , main =" Region de confiance bivariee 95%")

lines(ellipse(Sigma, centre=c(XT1,XT2)), type="l", col="red", xlab="Xt+1",ylab="Xt+2", main="Ellipse de confiance pour (Xt+1,Xt+2)")

abline(h=XT1, v=XT2)