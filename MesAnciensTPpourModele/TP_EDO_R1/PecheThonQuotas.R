m(list=ls())

#parms = c(r,K,Q)
#Valeurs positives par defaut : r = 2, K = 300
Eq3 = function(t,x,parms){
  dx1 = parms[1]*x[1]*(1-(x[1]/parms[2])) - parms[3]
  list(dx1)
}

temps=seq(0,2,by=0.1)
init= c(100)
r = 2
K = 300
Q=0.5

simulation1 = lsoda(init,temps,Eq3,parms=c(r,K,Q))

plot(temps,simulation1[,2],ylim=c(0,300),
     main = "PecheQuotas",
     ylab = "Densité de population (N)",
     lty=1,col=1,type= "l")
Q=10
simulation2 = lsoda(init,temps,Eq3,parms=c(r,K,Q))
lines(temps,simulation2[,2],lty = 2,col=2)
Q=100
simulation3 = lsoda(init,temps,Eq3,parms=c(r,K,Q))
lines(temps,simulation3[,2],lty = 2,col=3)
Q=150
simulation4 = lsoda(init,temps,Eq3,parms=c(r,K,Q))
lines(temps,simulation4[,2],lty = 2,col=4)

legend("bottomright",legend=c("Q=0.5","Q=10","Q=100","Q=150"),col=1:4,lty=1:4)


x11()
Q=10
curve(r*x*(1-x/K)-Q,xlim=c(-5,400),ylim=c(-100,200),main = "Q=10",ylab="dN/dt",xlab ="N",col = 1)
Q=100
curve(r*x*(1-x/K)-Q,xlim=c(-5,400),ylim=c(-100,200),main = "Q=100",ylab="dN/dt",xlab ="N",col = 2,add=TRUE)
Q=150
curve(r*x*(1-x/K)-Q,xlim=c(-5,400),ylim=c(-100,200),main = "Q=150",ylab="dN/dt",xlab ="N",col = 3,add=TRUE)
abline(h=0,lwd=3)
