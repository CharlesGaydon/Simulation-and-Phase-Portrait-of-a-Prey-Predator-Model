# Population Dynamic
# Author : Charles GAYDON
# Last edited : Second Semester 2016 (and 30/10/2017 for cleaning/translation).
# This work is a simulation of a Prey-Predator dynamic, based on a simple non geographical model,
# that took place in 2016.

################### Cleaning and loading ############################

rm(list=ls())
library(deSolve)
library(phaseR)
library(rootSolve)

################### Equations of the model ############################

#Format N1 : rootSolve.
EquationsF1 = function(t,y,parms){
  dNdP = list( c( y[1]*(2 - y[1]) - y[1]*y[2]/(1 + y[1]), y[1]*y[2]/(1 + y[1]) - parms[1]*y[2]))
}
#Format N2 : phaseR.
EquationsF2 <- function(t, y, parameters){
  dy <- numeric(2)
  dy[1] <- y[1]*(2 - y[1]) - y[1]*y[2]/(1 + y[1])
  dy[2] <- y[1]*y[2]/(1 + y[1]) - parameters[1]*y[2]
  list(dy)
}

################### Phase portrait ############################

phase <- function (mu, CI){
  title = paste("Phase portrait  - mu = ", mu," - CI= (", CI[1],",",CI[2],")")
    
  x.lim = c(0,2.5)
  y.lim = c(0,4)
  
  flow <- flowField(EquationsF2,x.lim,y.lim,parameters=mu,points=15,add=FALSE, 
                      xlab = "Preys", ylab = "Predator", main = title)
  null <- nullclines(EquationsF2,x.lim,y.lim,parameters=mu,add=TRUE)
  traj <- trajectory(EquationsF2,y0 = CI,parameters=mu,add=TRUE,t.end = 600)
}
################### Time-series ############################

chronique <- function (mu, init1, init2){
  temps=seq(0,200,by=0.1)
  s1=lsoda(y=init1,times=temps,func=EquationsF1,parms=mu)
  s2=lsoda(y=init2,times=temps,func=EquationsF1,parms=mu)

  plot(temps, s1[,2],col="dodgerblue",xlim=c(0,200), ylim=c(0,max(init1[2], init1[2])+3),type='l', ylab="Densit? de population", xlab='Temps', 
       main=paste("Chroniques - mu = ",round(mu,2)," - C1 = (", init1[1],",",init1[2],
                  ") et C2 = (", init2[1],",",init2[2],")"))
  lines(temps, s1[,3], col="indianred", type='l', lty = 1)
  lines(temps, s2[,2], col="dodgerblue", type='l', lty = 2)
  lines(temps, s2[,3], col="indianred", type='l', lty = 2)
  legend('topright',legend=c("Preys under C1",
                             "Predator under C1", 
                             "Preys under C2", 
                             "Predator under C2"),
                             lwd=1,col=c("indianred","dodgerblue","indianred",
                                         "dodgerblue"),lty=c(1,1,2,2), cex = 0.7 )
}
#chronique(0.3,c(0.5,2.2),c(4,5)) cool!!

################### Visualisation ############################

Graphiques <- function(){
  layout(matrix(c(1,2),nrow=1), widths=c(2,1))
  chronique(0.3,c(0.5,2.2),c(4,5))
  phase(0.3, c(4,5))
  chronique(0.4,c(0.5,2.2),c(4,5))
  phase(0.4, c(4,5))
  chronique(0.5,c(0.5,2.2),c(4,5))
  phase(0.5, c(4,5))
}

Graphiques()

# We highlight a bifurcation : a limit cycle slowly disappears. 
# There is an attracter around the non trivial equilibrium. The latter lost its instability when mu increases.

################### Equilibrium and stability ############################

#Returns three equilibrium and their stability :
# Stability == -1 ; Instability == 1.
#Structure of the dataframe returned is : x y stability
# To acces non trivial equilibrium : res[,3][,3]

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

###################### Phases portrait ########################

# We increment mu and draw equilibrium and their stability ;
# Two phase portrait : population density, and stability of non trivial equilibrium.

bifurcation <- function(){
  par(mfrow=c(1,1))
  
  # First phase portrait
  museq <- seq(0,2/3,by=0.001)
  plot(0,xlim=range(museq),ylim=c(-0.1,3),type="n",xlab="mu",ylab="N*",main="Bifurcation diagram")
  for (mu in museq) {
    st <- stability(mu)
    points(rep(mu,nrow(st)),st[,2],pch=22,col=c("darkred","black","lightpink")[st[,3]+2],bg =c("darkred","black","lightpink")[st[,3]+2])
    points(rep(mu,nrow(st)),st[,1],pch=22,col=c("darkblue","black","lightblue")[st[,3]+2],bg =c("darkblue","black","lightblue")[st[,3]+2])  
  }
  legend("topleft",pch=22,pt.cex=2,c("N* stable","N* unstable","P* stable","P* unstable"),
         col=c("darkblue","lightblue","darkred","lightpink"),pt.bg=c("darkblue","lightblue","darkred","lightpink"))
  
  # Second phase portrait - limit cycle is highlighted
  x=seq(0,2/3,by=0.001)
  y=rep(0,length(x))
  plot(x,y,type='l',lty=2,lwd=2,main = "Bifurcation diagram simplified",xlab="mu",ylab=" ",ylim=c(-0.1,0.8))
  x2=seq(1/3,2/3,by=0.001)
  y2=rep(0,length(x2))
  lines(x2,y2,type='l',lwd=2)
  x3=seq(0,1/3,by=0.001)
  y3 = sqrt(-x3+1/3)
  lines(x3,y3,type='l',lty=1,lwd=2)
}

bifurcation()