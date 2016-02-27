#TP de DYnamique des Population
#Christelle LOPES
#A rendre le 29/02/16 Ã  12h Ã  christelle.lopes@univ-lyon1.fr
#Overleaf :
# https://www.overleaf.com/4385545vxmdws#/13073637/

################### Nettoyage et chargement des packages ############################

rm(list=ls())
library(deSolve)
library(phaseR)
library(rootSolve)

################### Equations du modele ############################

#Format N°1 : rootSolve.
EquationsF1 = function(t,y,parms){
  dNdP = list( c( y[1]*(2 - y[1]) - y[1]*y[2]/(1 + y[1]), y[1]*y[2]/(1 + y[1]) - parms[1]*y[2]))
}
#Format N°2 : phaseR.
EquationsF2 <- function(t, y, parameters){
  dy <- numeric(2)
  dy[1] <- y[1]*(2 - y[1]) - y[1]*y[2]/(1 + y[1])
  dy[2] <- y[1]*y[2]/(1 + y[1]) - parameters[1]*y[2]
  list(dy)
}

################### Chroniques ############################

#Mise en évidence de la bifurcation : on constate la disparition progressive d'un cycle limite
#Attractant autour de l'équilibre non trivia, qui lui perd son instabilité avec 
#l'augmentation de mu.

phases <- function(){
  other = TRUE
  par(mfrow = c(2,2))
  for (mu in c(0.3,0.3,0.4,0.5)){
    title = paste("mu = ", mu)
    
    x.lim = c(0,2.5)
    y.lim = c(0,4)
  
    flow <- flowField(EquationsF2,x.lim,y.lim,parameters=mu,points=15,add=FALSE, 
                      xlab = "N", ylab = "P", main = title)
    null <- nullclines(EquationsF2,x.lim,y.lim,parameters=mu,add=TRUE)
    if (!other){
    init = c(1,4)
      traj <- trajectory(EquationsF2,y0 = init,parameters=mu,add=TRUE,t.end = 600)
      init = c(2.5,0.5)
      traj <- trajectory(EquationsF2,y0 = init,parameters=mu,add=TRUE,t.end = 600)
    }
    if (other){
      init = c(0.7,2)
      traj <- trajectory(EquationsF2,y0 = init,parameters=mu,add=TRUE,t.end = 600)
      other = FALSE
    }
  }
}

phases()

################### Equilibres et stabilite ############################

#Renvoit les trois équilibre du système ainsi que leur stabilité : 
#Stabilite : -1, Instabilite = 1.
#Structure du data.frame renvoyé : x y stabilité.
#Accès à la nature du point d'équilibre non trivial par resultat[,3][3].
stability = function(mu){
  d <- data.frame("x" = c(0,2,mu/(1-mu)), "y" = c(0,0,2 + (mu-2*mu*mu)/((1-mu)*(1-mu))), "Stability" = rep(1,3))
  for (i in 1:3){
    x = c(d[i,][1],d[i,][2])
    x = c(x[[1]],x[[2]])
    jacob = jacobian.full(x,EquationsF1,parms=mu)
    if ((det(jacob)>0)&& (jacob[1,1]+jacob[2,2])<0){
      d[,3][i] = -1
    }
  }
  return(d)
}

###################### Portrait de phases ########################
#On incrémente mu et on trace les équilibres et leur stabilité. 
#Traçons deux portraits de phase : celui des densités de population, puis celui plus concis 
# de la stabilité du point d'équilibre non trivial.
bifurcation <- function(){
  par(mfrow=c(1,1))
  
  # Premier portrait de phase
  museq <- seq(0,2/3,by=0.001)
  plot(0,xlim=range(museq),ylim=c(-0.1,3),type="n",xlab="mu",ylab="N*",main="Diagramme de bifurcation")
  for (mu in museq) {
    st <- stability(mu)
    points(rep(mu,nrow(st)),st[,2],pch=22,col=c("darkred","black","lightpink")[st[,3]+2],bg =c("darkred","black","lightpink")[st[,3]+2])
    points(rep(mu,nrow(st)),st[,1],pch=22,col=c("darkblue","black","lightblue")[st[,3]+2],bg =c("darkblue","black","lightblue")[st[,3]+2])  
  }
  legend("topleft",pch=22,pt.cex=2,c("N* stable","N* unstable","P* stable","P* unstable"),
         col=c("darkblue","lightblue","darkred","lightpink"),pt.bg=c("darkblue","lightblue","darkred","lightpink"))
  
  # Second portrait de phase - mettant en valeur l'apparition du cycle limite
  x=seq(0,2/3,by=0.001)
  y=rep(0,length(x))
  plot(x,y,type='l',lty=2,lwd=2,main = "Diagramme de bifurcation simplifié",xlab="mu",ylab=" ",ylim=c(-0.1,0.8))
  x2=seq(1/3,2/3,by=0.001)
  y2=rep(0,length(x2))
  lines(x2,y2,type='l',lwd=2)
  x3=seq(0,1/3,by=0.001)
  y3 = sqrt(-x3+1/3)
  lines(x3,y3,type='l',lty=1,lwd=2)
}

bifurcation()