rm(list=ls())
library(deSolve)
library(phaseR)
library(rootSolve)

#3
Equilibres = function(e){
  e1<-c(0,0)
  e2 = c(2,2-4*e)
  e3 = c(1/e,0)
  if (e == 0){ return(rbind(e1,e2))}
  if(e>=0.5){ return(rbind(e1,e3))}
  else {return(rbind(e1,e2,e3))}
}
#Modele de Bazikin
Bazikin = function(t,y,parms){
  dNdP = list(c(y[1]-y[1]*y[2]/(1+0.5*y[1])-parms[1]*y[1]*y[1], -y[2]*(1-y[1]/(1+0.5*y[1]))))
}

#Jacobienne pour epsilon = 0
print("Jacobienne pour epsilon = 0")
e = 1/4
parms = c(e)

#Renvoit en ligne pour chaque equilibre : x, y, stabilité
#Stabilité : -1, Instabilité = 1
Stability = function(e){
  Equi = Equilibres(e)
  n = length(Equi[,1])
  Stab = cbind(Equi,rep(1,n))
  for (i in 1:n){
    jacob = jacobian.full(Equi[i,],Bazikin,parms=e)
    if ((det(jacob)>0)&& (jacob[1,1]+jacob[2,2])<0){
      Stab[i,3]= -1
    }
  }
}

Stability(e)
####################"Chroniques##########################""
parms = c(0)
times=seq(0,200,0.25)
init = c(1,10)
sol1 = lsoda(init,times,Bazikin,parms=parms)
init = c(20,5)
sol2 = lsoda(init,times,Bazikin,parms=parms)

plot(sol1[,1],sol1[,2],type="l",col = 2,xlab="Temps",ylab="Densité de Population",main = "Chroniques - e = 0")
lines(sol1[,1],sol1[,3],type="l",col = 2,lty=2)
lines(sol2[,1],sol2[,2],type="l",col = 1)
lines(sol2[,1],sol2[,3],type="l",col = 1,lty=2)
legend("topright",legend = c("x : X0=(1,10)","y : X0=(1,10)","x : X0=(20,5)","y : X0=(20,5)"),col=c(2,2,1,1),lty=c(1,2,1,2))
#En l'absence d'une capacité limite 1/e, les proies et les prédateurs croissent exponentiellement.

####################"Portraits de Phase##########################""


BazikinPhase <- function(t, y, parameters){
  dy <- numeric(2)
  dy[1] <- y[1]-y[1]*y[2]/(1+0.5*y[1])-parameters[1]*y[1]*y[1]
  dy[2] <- -y[2]*(1-y[1]/(1+0.5*y[1]))
  list(dy)
}
x.lim = c(0,35)
y.lim = c(0,50)
flow <- flowField(BazikinPhase,x.lim,y.lim,parameters=parms,points=15,add=FALSE)
null <- nullclines(BazikinPhase,x.lim,y.lim,parameters=parms,add=TRUE)
init = c(1,10)
traj <- trajectory(BazikinPhase,y0 = init,parameters=parms,add=TRUE,t.end = 600)
init = c(20,5)
traj <- trajectory(BazikinPhase,y0 = init,parameters=parms,add=TRUE,t.end = 600)


#Cas de epsilon non nul
#Les conditions d'application du Th de PAH sont réunies, et
Stability(1/6)
#Instabilité du point d'equilibre non trivial a la bifurcation, d'ou une bifurcation
#de Hopf sous-critique.
x.lim = c(0,10)
y.lim = c(0,10)
flow <- flowField(BazikinPhase,x.lim,y.lim,parameters=parms,points=15,add=FALSE)
null <- nullclines(BazikinPhase,x.lim,y.lim,parameters=parms,add=TRUE)
init = c(2,1.4)

for (e in seq(0,1/6,0.01)){
  traj <- trajectory(BazikinPhase,y0 = init,parameters=e,add=TRUE,t.end = 200)
}

#Diagramme de Phase 
er = seq(0,1,0.01)
plot(0,xlim=range(er),ylim=c(0,10),type="n",main="Diagramme de bifurcation")
for (e in er ){
  s = Stability(e)
  n = nrow(s)
  for (i in 1:n){
    points(rep(e,nrow(s)),s[,1],col = c("darkblue","lightgrey","black")[s[,3]+2],pch=19)
  }
}

