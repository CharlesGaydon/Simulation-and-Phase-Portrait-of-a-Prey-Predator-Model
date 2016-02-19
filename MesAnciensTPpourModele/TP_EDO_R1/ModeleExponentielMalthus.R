rm(list=ls())    #nettoie l'espace de travail


#notre eq differentielle
#Les variables sont des vecteurs!
Eq1=function(t,x,parms){    
  dx1 = parms[1]*x[1]      #parms[1] = r ici 
  list(dx1)                #renvoit forcément une liste
}
#pas besoin de return car la fonction crée l'objet list et 
#lsoda attend cette liste


temps = seq(0,100,by=0.1)
init1=c(50)
r1=0.1
solution1=lsoda(y=init1,times=temps,func=Eq1,parms=c(r1))
r2=0
solution2=lsoda(y=init1,times=temps,func=Eq1,parms=c(r2))
r3=-0.1
solution3=lsoda(y=init1,times=temps,func=Eq1,parms=c(r3))

plot(temps,solution1[,2],ylim=c(0,120),type="l",
ylab = "Densité de population (N)",
main = "Modèle de Malthus pour différentes valeurs de r")
lines(temps,solution2[,2],ylim=c(0,120),col=2)      #col est la couleur
lines(temps,solution3[,2],ylim=c(0,120),col=3)

legend("topright",legend=c("r=0.1","r=0","r=0.1"),col=c(1,2,3),lty=1)  #topleft place la légende
#lty applique la couleur a la legende. Mais apparemment pour le type de ligne!