library(sp)
library(spdep)
library(raster)
library(geoR)
library(SpatialTools)
library(rgl)
library(rgdal)
library(gstat)
library(maptools)
library(MonteCarlo)
library(snowfall)
library(ggplot2)
library(rpart)
library(randomForest)

plume.creator=function (a){
library(gstat)
library(sp)
library(earth)
library(randomForest)
try({
method=sample(1:3, 1)
rel.diff=0.5
points=100
vol.plume=1000
#a=round(runif(1, min=2, max=20), 0)
#b=round(runif(1, min=2, max=20), 0)
rel.diff=0.5
points=100
vol.plume=1000
n=20
m=1
res=1
scale=10

val.x=runif(10, min=m, max=n)
val.y=runif(10, min=m, max=n)
val.z=runif(10, min=m, max=n)
val.c=runif(10, min=1300, max=10000)
dat=cbind(val.c, val.x, val.y, val.z)
colnames(dat)=c("conc", "x", "y", "z")
dat=data.frame(dat)

x.min=min(dat$x)
x.max=max(dat$y)
y.min=min(dat$y)
y.max=max(dat$y)
z.min=min(dat$z)
z.max=max(dat$z)


point1=c(5, 0, 0, z.min)
point2=c(5, 0, (y.max+y.max/2), z.min)
point3=c(5, (x.max+x.max/2), 0, z.min)
point4=c(5, (x.max+x.max/2), (y.max+y.max/2), z.min)
point5=c(5, 0, 0, (z.max+z.max/2))
point6=c(5, 0, (y.max+y.max/2), (z.max+z.max/2))
point7=c(5, (x.max+x.max/2), 0, (z.max+z.max/2))
point8=c(5, (x.max+x.max/2), (y.max+y.max/2), (z.max+z.max/2))
point9=c(5, ((x.max+x.min)/2), ((y.max+y.min)/2),((y.max+y.min)/2))

ref_points=rbind(point1, point2, point3, point4, point5, point6, point7, point8, point9)
colnames(ref_points)=c('conc', 'x', 'y', 'z')


x=seq(0, (x.max+y.max/2), by=res)
y=seq(0, (y.max+y.max/2), by=res)
z=seq(0, (z.max+z.max/2), by=res)

grid1=expand.grid(y, z)
colnames(grid1)=c("y","z")
grid1$x=0
grid1=grid1[,c(3,1,2)]

grid2=expand.grid(x, z)
colnames(grid2)=c("x","z")
grid2$y=0
grid2=grid2[,c(1,3,2)]


grid3=grid1
grid4=grid2

grid3$x=(x.max+x.max/2)
grid4$y=(y.max+y.max/2)

grid1=grid1[sample(nrow(grid1)),]
grid2=grid2[sample(nrow(grid2)),]
grid3=grid3[sample(nrow(grid3)),]
grid4=grid4[sample(nrow(grid4)),]

grid1=grid1[1:20,]
grid2=grid2[1:20,]
grid3=grid3[1:20,]
grid4=grid4[1:20,]

grid1=rbind(grid1, grid2, grid3, grid4)

temp=grid1
temp$conc=5

temp=temp[,c(4,1,2,3)]

colnames(temp)=colnames(dat)
dat=rbind(dat, temp)

x=seq(min(dat$x), max(dat$x), by=res)
y=seq(min(dat$y), max(dat$y), by=res)
z=seq(min(dat$z), max(dat$z), by=res)
prd.loc=expand.grid(x=x,y=y, z=z)
prd.loc=as.matrix(prd.loc)

if (method==1) {
prd.idw=prd.loc
prd.idw=data.frame(prd.idw)
coordinates(prd.idw)=~x+y+z
dat.idw=dat
coordinates(dat.idw)=~x+y+z
pred.idw=idw(conc~1, dat.idw, prd.idw, idp=sample(1:4, 1))
pred=pred.idw@data$var1.pred
} else if (method==2) {
#
#behrens et al method Spatial modelling with Euclidean distance fields and machine learning
distances=dist2(as.matrix(dat[,-1]), as.matrix(ref_points[,-1]))
conc=dat$conc
input_model=data.frame(cbind(conc, distances))
input_model=cbind(input_model, dat[,-1])
model.pred=earth(conc~., data=input_model)
distances2=dist2(as.matrix(ref_points[,-1]), as.matrix(prd.loc))
distances2=t(distances2)
distances2=cbind(distances2, prd.loc)
colnames(distances2)=colnames(input_model[,-1])
pred=predict(model.pred, newdata=distances2)
} else{
prd.idw=prd.loc
prd.idw=data.frame(prd.idw)
coordinates(prd.idw)=~x+y+z
dat.idw=dat
coordinates(dat.idw)=~x+y+z
range=((x.max+y.max))/sample(2:10, 1)
vario_model=sample(2:4, 1)
vario_model=vgm()[vario_model, 1]
m <- vgm(10, vario_model, range)
m=fit.variogram.reml(conc~1, dat.idw, model = m)
pred.idw=krige(conc~1, dat.idw, prd.idw, m)
pred=pred.idw@data$var1.pred
}

pred[pred<5]=5
k=cbind(pred, prd.loc)


dat2=cbind(pred, prd.loc)
dat2=data.frame(dat2)

dat2=subset(dat2, pred>1300)

dat2$z=round(dat2$z/scale,1)


#plot3d(dat2$x, dat2$y, dat2$z, aspect=1)
prd.loc.z=prd.loc[,3]/scale
vol.percent=nrow(dat2)/nrow(prd.loc)
vol.box=(max(x)-min(x))*(max(y)-min(y))*(max(prd.loc.z)-min(prd.loc.z))
vol.plume=vol.box*vol.percent

#kriging
z.val=round((z.max+z.max/2), 0)
x1=seq(0, (x.max+y.max/2), by=a)
y1=seq(0, (y.max+y.max/2), by=a)
z1=seq(0, z.val, by=a)
samples=expand.grid(x1, y1, z1)
colnames(samples)=c("x", "y", "z")
test=merge(samples, k, by=c("x", "y", "z"))
names=paste(test$x, test$y, test$z, sep="_")
test=aggregate(test[,4], list(test$x, test$y, test$z), mean)
colnames(test)=c("x", "y", "z", "pred")
input=test
prd.loc1=expand.grid(x, y)
colnames(prd.loc1)=c('x', "y")

distances=dist2(as.matrix(input[,-4]), as.matrix(ref_points[,-1]))
pred=test$pred
input_model=data.frame(cbind(pred, distances))
input_model=cbind(input_model, test[,-4])
model.pred=randomForest(pred~., data=input_model)
distances2=dist2(as.matrix(ref_points[,-1]), as.matrix(prd.loc))
distances2=t(distances2)
distances2=cbind(distances2, prd.loc)
colnames(distances2)=colnames(input_model[,-1])
pred=predict(model.pred, newdata=distances2)

maps=cbind(pred, prd.loc)
maps=data.frame(maps)
maps=subset(maps, pred>1300)
maps$z=maps$z/scale


points.total=nrow(prd.loc)
pred.vol.percent=nrow(maps)/points.total
pred.vol.plume=vol.box*pred.vol.percent
points=nrow(input)

diff=abs(vol.plume-pred.vol.plume)
rel.diff=abs(vol.plume-pred.vol.plume)/vol.plume
})
return(list("rel.diff"=rel.diff, 'vol.plume'=vol.plume, "points"=points))
}

a_grid=seq(2, 10, by=1)
param_list=list("a"=a_grid)

est=MonteCarlo(plume.creator, 10000, param_list, ncpus=1)

est$results

results=est$results$rel.diff
points.res=est$results$points
vol.plume=est$results$vol.plume




res.mean=apply(results, 1, function(x) mean(x))
points.mean=apply(points.res, 1, function(x) mean(x))
error.mean=qnorm(0.975)*res.mean/sqrt(ncol(results))

res.95=res.mean+error.mean
res.5=res.mean-error.mean

q95=apply(results, 1, function(x) quantile(x, 0.90))
q5=apply(results, 1, function(x) quantile(x, 0.10))


spacing=a_grid/20

fit=lm(log(q95)~log(points.mean))
fit.m=lm(log(res.mean)~log(points.mean))
fit.5=lm(log(q5)~log(points.mean))

points=10:900

q95.fit=exp(2.1214-0.5588*log(points))
mean.fit=exp(0.7055-0.4613*log(points))
q5.fit=exp(-0.8407-0.6204*log(points))


library(ggplot2)
dat=cbind(res.mean, spacing, points.mean, q95, q5)
dat=data.frame(dat)
dat.fit=cbind(q95.fit, mean.fit, q5.fit, points)
dat.fit=data.frame(dat.fit)

space.plot=ggplot(dat, aes(x=spacing, y=res.mean)) + 
  geom_smooth() + geom_errorbar(ymin=q5, ymax=q95) + ylim(0, 1.0) + 
  xlab("Borehole Spacing (%)") + ylab("Relative Error (%)")

points.plot=ggplot(dat, aes(x=points.mean, y=res.mean)) + 
  geom_point() + geom_errorbar(ymin=q5, ymax=q95) + ylim(0, 1.0) + xlab("Number of Data Points") + ylab("Relative Error (%)")


space.plot=ggplot(dat, aes(x=spacing, y=res.mean)) + 
  geom_smooth(se=FALSE) + geom_smooth(aes(x=spacing, y=q95, color="red"), se=FALSE) + 
  geom_smooth(aes(x=spacing, y=q5,color="red"), se=FALSE) +
  ylim(0, 1.0) + xlab("Borehole Spacing (%)") +    
  theme(legend.position="none")+
  ylab("Relative Error (%)") 

points.plot=ggplot(dat, aes(x=points.mean, y=res.mean)) + 
  geom_smooth(se=FALSE, span=0.3) + geom_smooth(aes(x=points.mean, y=q95, color="red"),span=0.3, se=FALSE) + 
  geom_smooth(aes(x=points.mean, y=q5,color="red"), se=FALSE) +
  ylim(0, 1.0) + xlab("Number of Samples") +    
  theme(legend.position="none")+
  ylab("Relative Error (%)") 

points.plot=ggplot(dat.fit, aes(x=points, y=mean.fit)) + 
  geom_line(col="blue") + geom_ribbon(aes(ymin=q5.fit, ymax=q95.fit), linetype=2, alpha=0.2) +
  ylim(0, 1.0) + xlab("Number of Samples") +
  xlim(20, 875) +
  geom_point(aes(x=points.mean, y=res.mean), data=dat) +
  geom_point(aes(x=points.mean, y=q95), shape=21, data=dat) +
  geom_point(aes(x=points.mean, y=q5), shape=21, data=dat) +
  theme(legend.position="none")+
  ylab("Relative Error (%)") 


#space.plot
points.plot

results.t=t(results)
colnames(results.t)=spacing

#optimization
q95
q5
res.mean

v=100


fit=lm(log(q95)~log(points.mean))
fit.m=lm(log(res.mean)~log(points.mean))
fit.5=lm(log(q5)~log(points.mean))

samples=1:3000
samples=log(samples)

pred=exp(2.8012-0.8937*samples)
pred=exp(pred)

plot(q95~points.mean)
lines(pred~exp(samples), col="red")
  
#s = samples
cost95=function(x) {a*x + b*v*(exp(2.1214-0.5588*log(x)))}
cost=function(x) {a*x + b*v*(exp(0.7055-0.4613*log(x)))}

#optim.100=optimize(cost, lower=1, upper=10000, maximum=FALSE)
a=10
b=120

results.optim.10=vector("list")
results.optim95.10=vector("list")
cost.optim.10=vector('list')
cost.optim95.10=vector('list')
for (i in 1:10000){
  v=i
  temp=optimize(cost, lower=1, upper=10000, maximum=FALSE)
  temp95=optimize(cost95, lower=1, upper=10000, maximum=FALSE)
  r1=temp$minimum
  r2=temp$objective
  results.optim.10[[i]]=r1
  cost.optim.10[[i]]=r2
  r1=temp95$minimum
  r2=temp95$objective
  results.optim95.10[[i]]=r1
  cost.optim95.10[[i]]=r2
  }

#$5 per sample
a=5
b=120

results.optim.5=vector("list")
results.optim95.5=vector("list")
cost.optim.5=vector('list')
cost.optim95.5=vector('list')
for (i in 1:10000){
  v=i
  temp=optimize(cost, lower=1, upper=10000, maximum=FALSE)
  temp95=optimize(cost95, lower=1, upper=10000, maximum=FALSE)
  r1=temp$minimum
  r2=temp$objective
  results.optim.5[[i]]=r1
  cost.optim.5[[i]]=r2
  r1=temp95$minimum
  r2=temp95$objective
  results.optim95.5[[i]]=r1
  cost.optim95.5[[i]]=r2
}

#$7.5 per sample
a=7.5
b=120

results.optim.7.5=vector("list")
results.optim95.7.5=vector("list")
cost.optim.7.5=vector('list')
cost.optim95.7.5=vector('list')
for (i in 1:10000){
  v=i
  temp=optimize(cost, lower=1, upper=10000, maximum=FALSE)
  temp95=optimize(cost95, lower=1, upper=10000, maximum=FALSE)
  r1=temp$minimum
  r2=temp$objective
  results.optim.7.5[[i]]=r1
  cost.optim.7.5[[i]]=r2
  r1=temp95$minimum
  r2=temp95$objective
  results.optim95.7.5[[i]]=r1
  cost.optim95.7.5[[i]]=r2
}

#100 per sample
a=100
b=120

results.optim.100=vector("list")
results.optim95.100=vector("list")
cost.optim.100=vector('list')
cost.optim95.100=vector('list')
for (i in 1:10000){
  v=i
  temp=optimize(cost, lower=1, upper=10000, maximum=FALSE)
  temp95=optimize(cost95, lower=1, upper=10000, maximum=FALSE)
  r1=temp$minimum
  r2=temp$objective
  results.optim.100[[i]]=r1
  cost.optim.100[[i]]=r2
  r1=temp95$minimum
  r2=temp95$objective
  results.optim95.100[[i]]=r1
  cost.optim95.100[[i]]=r2
}

#$5 per sample
a=50
b=120

results.optim.50=vector("list")
results.optim95.50=vector("list")
cost.optim.50=vector('list')
cost.optim95.50=vector('list')
for (i in 1:10000){
  v=i
  temp=optimize(cost, lower=1, upper=10000, maximum=FALSE)
  temp95=optimize(cost95, lower=1, upper=10000, maximum=FALSE)
  r1=temp$minimum
  r2=temp$objective
  results.optim.50[[i]]=r1
  cost.optim.50[[i]]=r2
  r1=temp95$minimum
  r2=temp95$objective
  results.optim95.50[[i]]=r1
  cost.optim95.50[[i]]=r2
}

#$75 per sample
a=75
b=120

results.optim.75=vector("list")
results.optim95.75=vector("list")
cost.optim.75=vector('list')
cost.optim95.75=vector('list')
for (i in 1:10000){
  v=i
  temp=optimize(cost, lower=1, upper=10000, maximum=FALSE)
  temp95=optimize(cost95, lower=1, upper=10000, maximum=FALSE)
  r1=temp$minimum
  r2=temp$objective
  results.optim.75[[i]]=r1
  cost.optim.75[[i]]=r2
  r1=temp95$minimum
  r2=temp95$objective
  results.optim95.75[[i]]=r1
  cost.optim95.75[[i]]=r2
}


#unlist and merge results

volume=1:10000
results.optim.10=unlist(results.optim.10)
cost.optim.10=unlist(cost.optim.10)
results.optim95.10=unlist(results.optim95.10)
cost.optim95.10=unlist(cost.optim95.10)

results.optim.5=unlist(results.optim.5)
cost.optim.5=unlist(cost.optim.5)
results.optim95.5=unlist(results.optim95.5)
cost.optim95.5=unlist(cost.optim95.5)

results.optim.7.5=unlist(results.optim.7.5)
cost.optim.7.5=unlist(cost.optim.7.5)
results.optim95.7.5=unlist(results.optim95.7.5)
cost.optim95.7.5=unlist(cost.optim95.7.5)

results.optim.100=unlist(results.optim.100)
cost.optim.100=unlist(cost.optim.100)
results.optim95.100=unlist(results.optim95.100)
cost.optim95.100=unlist(cost.optim95.100)

results.optim.50=unlist(results.optim.50)
cost.optim.50=unlist(cost.optim.50)
results.optim95.50=unlist(results.optim95.50)
cost.optim95.50=unlist(cost.optim95.50)

results.optim.75=unlist(results.optim.75)
cost.optim.75=unlist(cost.optim.75)
results.optim95.75=unlist(results.optim95.75)
cost.optim95.75=unlist(cost.optim95.75)


dat.res=cbind(results.optim.10, results.optim95.10, results.optim.5, results.optim95.5, results.optim.7.5, results.optim95.7.5,
          results.optim.100, results.optim95.100, results.optim.50, results.optim95.50, results.optim.75, results.optim95.75,
          volume)
dat.res=data.frame(dat.res)


optim.plot=ggplot(dat.res, aes(x=volume, y=results.optim.7.5)) + 
  geom_smooth(se=FALSE, span=0.3) + 
  geom_ribbon(aes(ymin=results.optim.10, ymax=results.optim.5), linetype=2, alpha=0.2) +
  geom_smooth(aes(x=volume, y=results.optim.75, col="red"), se=FALSE, span=0.3) + 
  geom_ribbon(aes(ymin=results.optim.100, ymax=results.optim.50), linetype=2, alpha=0.2) +
  xlab(expression("Plume Volume"~(m^3))) +    
  theme(legend.position="none")+
  ylab("Optimal Number of Samples") 


dat.cost=cbind(cost.optim.10, cost.optim95.10, cost.optim.5, cost.optim95.5, cost.optim.7.5, cost.optim95.7.5,
          cost.optim.100, cost.optim95.100, cost.optim.50, cost.optim95.50, cost.optim.75, cost.optim95.75,
          volume)

dat.cost=data.frame(dat.cost)

cost.plot=ggplot(dat.cost, aes(x=volume, y=cost.optim.7.5)) + 
  geom_smooth(se=FALSE, span=0.3) +
  geom_ribbon(aes(ymin=cost.optim.10, ymax=cost.optim.5), linetype=2, alpha=0.2) +
  geom_smooth(aes(x=volume, y=cost.optim.75, col="red"), se=FALSE, span=0.3) + 
  geom_ribbon(aes(ymin=cost.optim.100, ymax=cost.optim.50), linetype=2, alpha=0.2) +
  xlab(expression("Plume Volume"~(m^3))) +    
  theme(legend.position="none")+
  ylab("Total Data and Plume Error Cost ($)") 

#use geom_ribbon to add confidence envelopes around data

optim.plot
cost.plot

fit.cost=lm(log(results.optim)~log(volume))
plot(cost.optim~volume, type="l")

#optimal error plot
error.optim.7.5=exp(0.7055-0.4613*log(results.optim.7.5))
error.optim.5=exp(0.7055-0.4613*log(results.optim.5))
error.optim.10=exp(0.7055-0.4613*log(results.optim.10))

error.optim.75=exp(0.7055-0.4613*log(results.optim.75))
error.optim.50=exp(0.7055-0.4613*log(results.optim.50))
error.optim.100=exp(0.7055-0.4613*log(results.optim.100))



dat.error=cbind(error.optim.5, error.optim.7.5, error.optim.10, error.optim.50, error.optim.75, error.optim.100, volume)
dat.error=data.frame(dat.error)
dat.error=dat.error[10:100,]


error.plot=ggplot(dat.error, aes(x=volume, y=error.optim.7.5)) + 
  geom_smooth(se=FALSE, span=0.3) +
  geom_ribbon(aes(ymin=error.optim.10, ymax=error.optim.5), linetype=2, alpha=0.2) +
  geom_smooth(aes(x=volume, y=error.optim.75, col="red"), se=FALSE, span=0.3) + 
  geom_ribbon(aes(ymin=error.optim.100, ymax=error.optim.50), linetype=2, alpha=0.2) +
  xlab(expression("Plume Volume"~(m^3))) +    
  theme(legend.position="none")+
  ylab("Optimal Error (%)") 

error.plot
