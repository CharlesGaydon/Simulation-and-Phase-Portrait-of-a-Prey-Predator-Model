rm(list=ls())

#set working place : windows
setwd("H:/Travail Scolaire/INSA/3?me ann?e/EDO/Lopes")
#en salle de cours
setwd("/run/media/alhostis/DD Adele/Travail Scolaire/INSA/3ème année/EDO/Lopes")

#install the packages and input them
library(rootSolve)
library(deSolve)
library(phaseR)
library(scatterplot3d)


#Création fonctions
Holling1=function(t,x,parameters)
{
  dx1=x[1]*(2-x[1])-(x[1]*x[2])/(1+x[1])
  dx2=(x[1]*x[2])/(1+x[1])-parameters[1]*x[2]
  list(c(dx1,dx2))
}
#Construction de la même fonction pour le flowfield:
Holling2=function(t,y,parameters)
{
  dy1=y[1]*(2-y[1])-(y[1]*y[2])/(1+y[1]) 
  dy2=(y[1]*y[2])/(1+y[1])-parameters[1]*y[2]
  list(c(dy1,dy2))
}


#Ressortir les équilibres
Equilibrium<-function(mu){
  p1<-c(0,0)
  p2<-c(2,0) 
  p3<-c(mu/(1-mu),(2-3*mu)/(1-mu)^2)
  if(mu>0 && mu<(2/3)){eq<-rbind(p1,p2,p3)}
  return(eq)
}


#Ressortir la stabilité
Stability<-function(mu){
  Eig<-vector(length=nrow(Eq))
  Del<-vector(length=nrow(Eq))
  for (i in 1:nrow(Eq)){
    Jac<-jacobian.full(Eq[i,], func=Holling1, time = 0, parms = c(mu))
    Tr<-Jac[1,1]+Jac[2,2]
    Det<-Jac[1,1]*Jac[2,2]-Jac[1,2]*Jac[2,1]
    delta<-Tr^(2)-4*Det
    if (Det<0 && Tr!=0){stab<-"PS"}  #Point Selle
    if (Tr==0){stab<-"C"}   #Centres
    if(delta>0 && Det>0 && Tr<0){stab<- "NAS"}  #Noeud A. Stable
    if(delta<0 && Det>0 && Tr<0){stab<- "FAS"}  #Foyer A. Stable
    if (delta>0 && Det>0 && Tr>0){stab<- "NI"}  #Noeud Instable
    if (delta<0 && Det>0 && Tr>0){stab<- "FI"}  #Foyer Instable
    Eig[i]<-stab
    Del[i]<-delta
    print (paste("Jacobienne du point d'équilibre ",i))
    print(Jac)
  }
  print("Et ici leur stabilité:")
  return(cbind(Eq,Eig))
}

#Fonction graphique
Graphique<-function(mu, CI_1, CI_2){
  par(mfrow=c(1,3),cex.main=1.7)
  Eq<-Equilibrium(mu)
  temps=seq(0,1000,by=0.1)
  sln1=lsoda(y=CI_1,times=temps,func=Holling1,parms=c(mu))
  sln2=lsoda(y=CI_2,times=temps,func=Holling1,parms=c(mu))
  plot(sln1[,3]~sln1[,2],col="dodgerblue3",xlim=c(0,max(CI_1[1],CI_2[1])+1), ylim=c(0,max(CI_1[2],CI_2[2])+1),type='l', ylab="Prédateurs", xlab='Proies', 
       main=paste("Portrait de phase du modèle Holling pour mu = ",round(mu,4)))
  lines(sln2[,3]~sln2[,2],lty=2,col="indianred")
  legend('topright',legend=c("Simulation sous CI_1",
                             "Simulation sous CI_2",
                             "Origine simulation sous CI_1",
                             "Origine simulation sous CI_2",
                             "Points d'équilibres"),
         lwd=2,col=c("dodgerblue","indianred","dodgerblue","indianred","black")
         ,lty=c(1,2,NA,NA,NA),pch=c(NA,NA,4,4,4))
  flowField(deriv=Holling2, x.lim=c(0,25),y.lim=c(0,30),
            parameters=c(mu),points=25)
  
  trajectory(deriv=Holling1, y0 =c(30,15), n = NULL, t.start = 0,
             t.end=200, t.step = 0.01,parameters = c(mu), 
             system = "two.dim", colour = "dodgerblue",lwd=1,lty=2,
             add = TRUE)  
  #Ajout des pts d'eq et des pts de CI sur le schéma
  points(c(CI_1[1],CI_2[1]),c(CI_1[2],CI_2[2]),col=c("dodgerblue",
                                                     "indianred"),
         lwd=3,pch=4)
  for (i in 1:nrow(Eq)){
    points(Eq[i,1],Eq[i,2],col=1,lwd=3,pch=4)
  }
  #Et maintenant, les chroniques
  plot(temps,sln1[,2],type='l',col="indianred",ylab="Densité de pop", ylim=c(0,max(CI_1[1],CI_2[1])+1), xlim=c(0,100),
       xlab='Temps',main=paste("Chroniques du modèle Holling pour mu = ", 
                  round(mu,4)))
  lines(temps,sln1[,3],lty=1,col="dodgerblue")
  lines(temps,sln2[,2],lty=2,col="indianred")
  lines(temps,sln2[,3],lty=2,col="dodgerblue")
  legend('topright',legend=c("Proies sous CI_1",
                             "Prédateurs sous CI_1", 
                             "Proies sous CI_2", 
                             "Prédateurs sous CI_2"),
         lwd=1,col=c("indianred","dodgerblue","indianred",
                     "dodgerblue"),lty=c(1,1,2,2) )
  
  #Graphe en 3D
  x = c(sln1[,2],sln2[,2])
  y = c(sln1[,3],sln2[,3])
  t = c(sln1[,1],sln2[,1])
  color <-c(rep("indianred", length(sln1[,2])),rep("dodgerblue3", length(sln1[,2])))
  stp=scatterplot3d(t,x,y,xlab="Temps",xlim=c(0,100),ylab="Pred",zlab="Proies",lwd=2, color,main=paste("Evolution des populationsen dans le temps ; mu = ", round(mu,4)))
  legend("topright", legend = c("Simulation sous CI_1","Simulation sous CI_2"),col =  c("dodgerblue3","indianred"), lty=c(1,1),lwd=c(2,2))
  stp
  
  par(mfrow=c(1,1),cex.main=1)
}

#Fixons les CI pour nos simulations
CI_1 = c(4,1)
CI_2 = c(1,10)

#2e
mu=1/6
CI_1 = c(4,5)
CI_2 = c(1,2)
Graphique(mu,CI_1,CI_2)

#3e
mu = 1/3
Graphique(mu,CI_1,CI_2)

#4e
mu=2/3-0.2
Graphique(mu,CI_1,CI_2)

#\\\\\\\\\\\\\\\\\\\\\\\\\\_

#Diagrammes de Bifurcation

Equi <- function(mu) {
  p1<-c(0,0)
  p2<-c(2,0)
  p3<-c(mu/(1-mu),(2-3*mu)/(1-mu)^2)
  rbind(p1,p2,p3)
}
stability <- function(mu) {
  Eq <- Equi(mu)
  Eigen<-vector(length=nrow(Eq))
  for (i in 1:nrow(Eq)){
    Jac<-jacobian.full(Eq[i,], func=Holling1, time = 0, parms = c(mu))
    Tr <- Jac[1,1]+Jac[2,2]
    Det <- det(Jac)
    if (Det<0){stab <- 1}
    else{if (Tr<0){stab <- -1}
         else{stab <- 1}}
    Eigen[i]<-stab}
  return(cbind(Eq,Eigen))
}
museq <- seq(0,2/3,by=0.001)
plot(0,xlim=range(museq),ylim=c(-0.1,3),type="n",xlab="mu",ylab="N*",main="Diagramme de bifurcation")
for (mu in museq) {
  st <- stability(mu)
  #je me sers de l'indice de stabilité pour choisir la "bonne" couleur: si indice=-1 (stable) -1 + 2 = 1 et il prend la première couleur de mon vecteur col=foncé
  #si indice=1 (instable) 1+2=3 et il prend la troisième couleur de mon vecteur col=clair
  #je trace d'abord N, puis P
  points(rep(mu,nrow(st)),st[,2],pch=22,col=c("darkred","black","lightpink")[st[,3]+2],bg =c("darkred","black","lightpink")[st[,3]+2])
  points(rep(mu,nrow(st)),st[,1],pch=22,col=c("darkblue","black","lightblue")[st[,3]+2],bg =c("darkblue","black","lightblue")[st[,3]+2])  
}
legend("topleft",pch=22,pt.cex=2,c("N* stable","N* unstable","P* stable","P* unstable"),
       col=c("darkblue","lightblue","darkred","lightpink"),pt.bg=c("darkblue","lightblue","darkred","lightpink"))


#Construction notre diagramme de bifurcation.
x=seq(0,2/3,by=0.001)
y=rep(0,length(x))
x2=seq(1/3,2/3,by=0.001)
y2=rep(0,length(x2))
plot(x,y,type='l',lty=2,lwd=2,main = "diagramme de bifurcation",xlab="mu",ylab=" ",ylim=c(-0.1,0.8))
lines(x2,y2,type='l',lwd=2)
x3=seq(0,1/3,by=0.001)
y3 = sqrt(-x3+1/3)
lines(x3,y3,type='l',lty=1,lwd=2)

