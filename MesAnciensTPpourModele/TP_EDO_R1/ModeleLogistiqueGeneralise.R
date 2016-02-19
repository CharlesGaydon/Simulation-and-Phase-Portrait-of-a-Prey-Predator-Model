rm(list = ls())

#Verhulst
#parms = c(r,K,theta)
#x = c(N)
Eq2 = function(t,x,parms){
  dx1 = parms[1]*x[1]*(1-(x[1]/parms[2])^parms[3])
  list(dx1)
}


#1
temps = seq(0,30,by=0.1)
initi=c(50)
r1 = 0.5
K1 = 300
theta=0.5
solution1 = lsoda(initi,temps,Eq2,parms=c(r1,K1,theta))
#2
initi=c(50)
r2 = 0.5
K2 = 300
theta=1
solution2 = lsoda(initi,temps,Eq2,parms=c(r2,K2,theta))
#2
initi=c(50)
r3 = 0.5
K3 = 300
theta=2
solution3 = lsoda(initi,temps,Eq2,parms=c(r3,K3,theta))
#4
initi=c(50)
r4 = 0.5
K4 = 300
theta=10
solution4 = lsoda(initi,temps,Eq2,parms=c(r4,K4,theta))
#5
initi=c(50)
r5 = -0.5
K5 = 300
theta=0.5
solution5 = lsoda(initi,temps,Eq2,parms=c(r5,K5,theta))


#Tracer les courbes

plot(solution1[,1],solution1[,2],ylim=c(0,350), type = "l",
     xlab = "Temps",
     main = "Modèle de Malthus pour différentes valeurs de r",
     ylab = "Densité de population (N)",
     col = 1,
     lty = 1
     )
lines(solution2[,1],solution2[,2],ylim=c(0,350),col=2,lty =1,lwd=2)
lines(solution3[,1],solution3[,2],ylim=c(0,350),col=3,lty = 3,lwd = 2)
lines(solution4[,1],solution4[,2],ylim=c(0,350),col=4,lty = 4,lwd = 2)
lines(solution5[,1],solution5[,2],ylim=c(0,350),col=5,lty = 5)


legend("bottomright",legend=c("N0=50,r=0.5,theta=0.5","N0=50,r=1,theta=1","N0=50,r=2,theta=2",
                           "N0=400,r=0.5,theta=10","N0=50,r=0.5,theta=0.5"),col=1:5,lty=c(1,1,3,4,5),lwd=c(1,2,2,2,1)) 

#lty donne le type
#lwd donne l'épaisseur


#Tracons dN/dt

par(mfrow=c(1,2))     #partition de la fenêtre graphique
#on peut aussi faire x11() qui ouvre une deuxieme fenetre "Device"
curve(0.5*x*(1-(x/300)^2),xlim=c(-100,400),ylim=c(-10,70),main = "r = 0.5",ylab="dN/dt",xlab ="N") 
#trace premier argument qui est une fonction prenant x en parametre
abline(h=0,lty=2)
curve(-0.5*x*(1-(x/300)^2),xlim=c(-100,400),ylim=c(-70,10),main = "r = -0.5",ylab="dN/dt",xlab ="N")
abline(h=0,lty=2)
