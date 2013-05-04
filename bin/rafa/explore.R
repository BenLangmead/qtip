load("for_rafa.Rda")

p=repCor.s.50/repTot.s.50
tt=repTot.s.50
x=as.numeric(rownames(p))
y=as.numeric(colnames(p))
y[y==10000]<-100
plot(x,rowMeans(p,na.rm=TRUE))
plot(y,colMeans(p,na.rm=TRUE))

##after data exploration noticed that
###if y>25 then it is pretty much 1, so let's get rid of them
K=24
repCor.s.50[,y==K]<-rowSums(repCor.s.50[,y>=K])
repCor.s.50<-repCor.s.50[,y<=K]
repTot.s.50[,y==K]<-rowSums(repTot.s.50[,y>=K])
repTot.s.50<-repTot.s.50[,y<=K]
p=repCor.s.50/repTot.s.50
tt=repTot.s.50
x=as.numeric(rownames(p))
y=as.numeric(colnames(p))



tmp=expand.grid(x,y)
X=tmp[,1]
Y=tmp[,2]
fit=loess(as.vector(p)~X*Y,span=1/2,degree=1)
z=matrix(predict(fit,newdata=data.frame(X=X,Y=Y)),nrow(p),ncol(p))
z[z>1]<-1;z[z<0]<-0

###z is the smoothed map. upi can use predict above to create map

image(-x,y,z,col=rev(brewer.pal(9,"Blues")))

###look at countour of function of x stratified by y
cols=rgb(seq(.9,.1,len=ncol(z)),0,seq(.1,.9,len=ncol(z)))
matplot(x,z,col=cols,lty=1,type="l")

###look at countour of function of y stratified by x
cols=rgb(seq(.9,.1,len=nrow(z)),0,seq(.1,.9,len=nrow(z)))
matplot(y,t(z),col=cols,lty=1,type="l")




####other stuff i looked at
if(FALSE){
  
mypar(5,5)
smooth<-matrix(NA,nrow(p),ncol(p))
for(i in 1:ncol(p)){
  plot(x,p[,i],ylim=c(0,1))
  fit=loess(p[,i]~x,span=3/4,weight=sqrt(tt[,i]),degree=1)
  lines(fit$x,fit$fitted,col="red")
  smooth[,i]<-predict(fit,newdata=data.frame(x=x))
}

mypar(6,5)
#smooth<-matrix(NA,nrow(p),ncol(p))
for(i in 1:nrow(p)){
  plot(y,p[i,],ylim=c(0,1))
  fit=loess(p[i,]~y,span=3/4,weight=sqrt(tt[i,]),degree=1)
  lines(fit$x,fit$fitted,col="red")
 # smooth[,i]<-predict(fit,newdata=data.frame(x=x))
}


