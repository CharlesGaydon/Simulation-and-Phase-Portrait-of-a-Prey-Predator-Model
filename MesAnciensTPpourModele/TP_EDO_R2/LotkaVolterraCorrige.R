rm(list=ls())
library(deSolve)
library(phaseR)

init = c(1,10)
times = seq(0,1000,0.1)

r=0.2
a=0.2
e=1
m=0.4
K=10
parms= c(r,a,e,m,K)
par(mfrow=c(1,2))

plot(1:10,1:10,xlim=c(0,20),ylim=c(0,20),col = 0)
derivativ <- function(t, y, parameters){
  dy <- numeric(2)
  dy[1] <- y[1]*parameters[1]*(1 - y[1]/parameters[5])-parameters[2]*y[1]*y[2]
  dy[2] <- parameters[3]*parameters[2]*y[1]*y[2] - parameters[4]*y[2]
  list(dy)
}
x.lim = c(0,4)
y.lim = c(0,4)
flow <- flowField(derivativ,x.lim,y.lim,parameters=parms,points=15,add=FALSE)
null <- nullclines(derivativ,x.lim,y.lim,parameters=parms,add=TRUE)
traj <- trajectory(derivativ,y0 = init,parameters=parms,add=TRUE,t.end = 600)


#Chroniques
Eq1 = function(t,y,parms){
  dNdP = list(c(y[1]*parms[1]*(1 - y[1]/parms[5])-parms[2]*y[1]*y[2], parms[3]*parms[2]*y[1]*y[2] - parms[4]*y[2]))
}
sol1 = lsoda(init,times,Eq1,parms=parms)
sol1

plot(sol1[,1],sol1[,2],main = "Chroniques Proie-Prédateur",type = "l",xlim=c(0,300),xlab="Temps",ylab="Densité de Population")
lines(sol1[,1],sol1[,3],col=2,lty=2)
legend("topright",legend=c("N","P"),col=c(1,2),lty=1)



#install.packages("rootSolve")
library(rootSolve)
###Equilibres
#Renvoit sur chaque ligne les equilibres et en 3e colonne 4 si stable, 1 sinon
equi <-function(parms){
  r=parms[1]
  a=parms[2]
  e=parms[3]
  m=parms[4]
  K=parms[5]
  
  if (K > m/(e*a)){
    
    y=list(c(0,0),c(K,0),c(m/(e*a),(r/a)*(1-m/(K*e*a))))    
    x=list(c(0,0,1),c(K,0,1),c(m/(e*a),(r/a)*(1-m/(K*e*a)),1))
    EQUI=matrix(ncol=3,nrow=3)
    for (i in 1:3){
      EQUI[i,]=x[[i]]
      jacob = jacobian.full(y[[i]],Eq1,parms=parms)
      if ((det(jacob)>0) && (jacob[1,1]+jacob[2,2])<0){
        EQUI[i,3]= 4
      }
    }
  } else {  
    y=list(c(0,0),c(K,0))   
    x=list(c(0,0,1),c(K,0,1))
    EQUI = matrix(nrow = 2,ncol = 3)
    
    for (i in 1:2){
      EQUI[i,]=x[[i]]
      jacob = jacobian.full(y[[i]],Eq1,parms=parms)
      if ((det(jacob)>0 )&& (jacob[1,1]+jacob[2,2])<0){
        EQUI[i,3]= 4
      }
    }
  }
  return(EQUI)
}

#Bifurcations
r=0.2
a=0.2
e=1
m=0.4

Krange = seq(0,10,by=0.25)
plot(0,xlim=range(Krange),ylim=range(Krange),type="n",main="Diagramme de bifurcation")
for (k in Krange){
    equilibres <- equi(c(r,a,e,m,k))
    points(rep(k,nrow(equilibres)),equilibres[,1],col = c("darkblue","lightgrey","black")[equilibres[,3]+2])
    points(rep(k,nrow(equilibres)),equilibres[,2],col = c("darkblue","lightgrey","black")[equilibres[,3]+2])
    #points(k,equilibres[i,2],col = c)

}

