rm(list=ls())
library(deSolve)
library(phaseR)
#Paramêtres
r=0.8
a=0.2
e=0.2
m=0.4
#Système Lotka Volterra- Reponse Holling I
Eq1 = function(t,x,parms){
  dNdP = list(c(parms[1]*x[1]-parms[2]*x[1]*x[2], parms[3]*parms[2]*x[1]*x[2] - parms[4]*x[2]))
}
parms=c(r,a,e,m)
#1
init= c(0.1,0.1)
times = seq(0,1000,0.1)
sol1 = lsoda(init,times,Eq1,parms=parms)
sol1
plot(sol1[,2],sol1[,3], type = "l",
     xlab = "Densité de population (N)",
     main = "Modèle de MLotka-Volterra",
     ylab = "Densité de population (P)",
     col = 1,
     lty = 1,
     ylim=c(0,45)
)

arrows(sol1[5,2],sol1[5,3],sol1[35,2],sol1[35,3])
#2
init= c(40,0.1)
sol1 = lsoda(init,times,Eq1,parms=c(r,a,e,m))
lines(sol1[,2],sol1[,3],col=2,lty =1,lwd=2)
#3
init= c(20,30)
sol1 = lsoda(init,times,Eq1,parms=c(r,a,e,m))
lines(sol1[,2],sol1[,3],col=3,lty =1,lwd=2)
arrows(sol1[5,2],sol1[5,3],sol1[25,2],sol1[25,3])

for (i in seq(0,120,20) ){
  for (j in seq(0,40,10)){
    dN=Eq1(c(1),c(i,j),c(r,a,e,m))[[1]][1]
    dP=Eq1(c(1),c(i,j),c(r,a,e,m))[[1]][2]
    print(dN)
    arrows(i,j,i+dN*0.1,j+dP*0.1)
  }
}

legend("bottom",legend = c("x0=(0.1,0.1)","x0=(40,0.1)","x0=(100,30)"),col=c(1,2,3))

par(mfrow=c(1,2))

derivative <- function(t, y, parameters){
  dy <- numeric(2)
  dy[1] <- parameters[1]*y[1]-parameters[2]*y[1]*y[2]
  dy[2] <- parameters[3]*parameters[2]*y[1]*y[2] - parameters[4]*y[2]
  list(dy)
}
x.lim = c(0,40)
y.lim = c(0,101)
flow <- flowField(derivative,x.lim,y.lim,parameters=parms,points=15,add=FALSE)
null <- nullclines(derivative,x.lim,y.lim,parameters=parms,add=TRUE)
traj <- trajectory(derivative,y0 = c(10,20),parameters=parms,add=TRUE,t.end = 50)
traj <- trajectory(derivative,y0 = c(20,10),parameters=parms,add=TRUE,t.end = 50)
#Chroniques
par(mfrow=c(1,1))
plot(sol1[,1],sol1[,2],xlim=c(0,75),col = 1,type = "l")
lines(sol1[,1],sol1[,3],xlim=c(0,75),col = 3,type = "l",lwd = 2)


