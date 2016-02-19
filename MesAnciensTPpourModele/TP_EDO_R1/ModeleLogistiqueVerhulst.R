rm(list = ls())

#Verhulst
#parms = c(r,K)
#x = c(N)
Eq2 = function(t,x,parms){
  dx1 = parms[1]*x[1]*(1-x[1]/parms[2])
  list(dx1)
}
#1
temps = seq(0,30,by=0.1)
initi=c(50)
r1 = 0.5
K1 = 300
solution1 = lsoda(initi,temps,Eq2,parms=c(r1,K1))
#2
initi=c(50)
r2 = 1
K2 = 300
solution2 = lsoda(initi,temps,Eq2,parms=c(r2,K2))
#2
initi=c(50)
r3 = 2
K3 = 300
solution3 = lsoda(initi,temps,Eq2,parms=c(r3,K3))
#4
initi=c(400)
r4 = 0.5
K4 = 300
solution4 = lsoda(initi,temps,Eq2,parms=c(r4,K4))
#5
initi=c(400)
r5 = 0.5
K5 = 400
solution5 = lsoda(initi,temps,Eq2,parms=c(r5,K5))
#6
initi=c(50)
r6 = 0.5
K6 = 400
solution6= lsoda(initi,temps,Eq2,parms=c(r6,K6))
#7
initi=c(50)
r7 = -0.5
K7= 300
solution7 = lsoda(initi,temps,Eq2,parms=c(r7,K7))

#Tracer les courbes

plot(solution1[,1],solution1[,2],ylim=c(0,500), type = "l",
     xlab = "Temps",
     main = "Modèle de Malthus pour différentes valeurs de r",
     ylab = "Densité de population (N)",
     col = 1,
     lty = 1
     )
lines(solution2[,1],solution2[,2],ylim=c(0,500),col=2,lty =1)
lines(solution3[,1],solution3[,2],ylim=c(0,500),col=3,lty = 1,lwd = 2)
lines(solution4[,1],solution4[,2],ylim=c(0,500),col=4,lty = 1,lwd = 2)
lines(solution5[,1],solution5[,2],ylim=c(0,500),col=5,lty = 5)
lines(solution6[,1],solution6[,2],ylim=c(0,500),col=6,lty = 6)
lines(solution7[,1],solution7[,2],ylim=c(0,500),col=7,lty = 7)

legend("bottomright",legend=c("N0=50,r=0.5,K=300","N0=50,r=1,K=300","N0=50,r=2,K=300",
                           "N0=400,r=0.5,K=300","N0=50,r=0.5,K=400","N0=50,r=0.5,K=400",
                           "N0=50,r=-0.5,K=300"),col=1:7,lty=c(1,1,1,1,5,6,7),lwd=c(1,1,2,2,1,1,1)) 

#lty donne le type
#lwd donne l'épaisseur
