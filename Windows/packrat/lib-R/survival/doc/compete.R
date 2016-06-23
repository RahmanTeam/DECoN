### R code from vignette source 'compete.Rnw'

###################################################
### code chunk number 1: compete.Rnw:20-25
###################################################
options(continue="  ", width=60)
options(SweaveHooks=list(fig=function() par(mar=c(4.1, 4.1, .3, 1.1))))
pdf.options(pointsize=10) #text in graph about the same as regular text
options(contrasts=c("contr.treatment", "contr.poly")) #ensure default
require("survival")


###################################################
### code chunk number 2: sfig1
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mar=c(.1, .1, .1, .1))
frame()
par(usr=c(0,100,0,100))
# first figure
xx <- c(0, 10, 10, 0)
yy <- c(0, 0, 10, 10)
polygon(xx +10, yy+70)
polygon(xx +30, yy+70)
arrows( 22, 75, 28, 75, length=.1)
text(c(15, 35), c(75,75), c("Alive", "Dead"))

# second figure
polygon(xx +60, yy+70)  
for (j in c(55, 70, 85)) {
    polygon(xx +80, yy+j)
    arrows(72, (5*75 +j+5)/6, 78, (100+j*5)/6, length=.1)
}
text(c(65, 85,85,85), c(70,55,70,85)+5, c("A", "D1", "D2", "D3")) 

# third figure
polygon(xx+20, yy+25)
for (j in c(15,35)) {
    polygon(xx +40, yy+j)
    arrows(32, (5*30 +j+4)/6, 38, (54+j*5)/6, length=.1)
}
arrows(38, 2+(55 + 35*5)/6, 32, 2+ (150 + 40)/6, length=.1)
arrows(45, 33, 45, 27, length=.1)
text(c(25, 45,45), c(30, 20, 40), c("Health", "Death", "Illness"))


###################################################
### code chunk number 3: mgus1
###################################################
getOption("SweaveHooks")[["fig"]]()
oldpar <- par(mfrow=c(1,2))
hist(mgus2$age, nclass=30, main='', xlab="Age")
with(mgus2, tapply(age, sex, mean))

mfit1 <- survfit(Surv(futime, death) ~ sex, data=mgus2)
mfit1
plot(mfit1, col=c(1,2), xscale=12, mark.time=FALSE, lwd=2,
     xlab="Years post diagnosis", ylab="Survival")
legend(6, .8, c("female", "male"), col=1:2, lwd=2, bty='n')
par(oldpar)


###################################################
### code chunk number 4: mgus2
###################################################
getOption("SweaveHooks")[["fig"]]()
etime <- with(mgus2, ifelse(pstat==0, futime, ptime))
event <- with(mgus2, ifelse(pstat==0, 2*death, 1))
event <- factor(event, 0:2, labels=c("censor", "pcm", "death"))
table(event)

mfit2 <- survfit(Surv(etime, event) ~ sex, data=mgus2)
mfit2
plot(mfit2, col=c(1,1,2,2), lty=c(2,1,2,1),
     xscale=12, mark.time=FALSE, lwd=2, 
     xlab="Years post diagnosis", ylab="Prevalence")
legend(20, .6, c("death:female", "death:male", "pcm:female", "pcm:male"), 
       col=c(1,2,1,2), lty=c(1,1,2,2), lwd=2, bty='n')


###################################################
### code chunk number 5: mgus3
###################################################
getOption("SweaveHooks")[["fig"]]()
pcmbad <- survfit(Surv(etime, pstat) ~ sex, data=mgus2)
plot(pcmbad[2], mark.time=FALSE, lwd=2, fun="event", conf=FALSE, xscale=12,
     xlab="Years post diagnosis", ylab="Fraction with PCM")
lines(mfit2[2,1], lty=2, lwd=2, mark.time=FALSE, conf=FALSE, xscale=12)
legend(0, .28, c("Males, PCM, incorrect curve", "Males, PCM, competing risk"),
       col=1, lwd=2, lty=c(1,2), bty='n')


###################################################
### code chunk number 6: mgus4
###################################################
ptemp <- with(mgus2, ifelse(ptime==futime & pstat==1, ptime-.1, ptime))
newdata <- tmerge(mgus2, mgus2,  id=id, death=event(futime, death))
newdata <- tmerge(newdata, mgus2, id, pcm = event(ptemp, pstat))
newdata <- tmerge(newdata, newdata, id, enum=cumtdc(tstart))
with(newdata, table(death, pcm))


###################################################
### code chunk number 7: mgus4g
###################################################
getOption("SweaveHooks")[["fig"]]()
temp <- with(newdata, ifelse(death==1, 2, pcm))
newdata$event <- factor(temp, 0:2, labels=c("censor", "pcm", "death"))  
mfit3 <- survfit(Surv(tstart, tstop, event) ~ sex, data=newdata, id=id)
plot(mfit3[,1], mark.time=FALSE, col=1:2, lty=1, lwd=2,
     xscale=12,
     xlab="Years post MGUS diagnosis", ylab="Prevalence of PCM")
legend(4, .04, c("female", "male"), lty=1, col=1:2, lwd=2, bty='n') 


###################################################
### code chunk number 8: mgus5
###################################################
getOption("SweaveHooks")[["fig"]]()
d2 <- with(newdata, ifelse(enum==2, 4, as.numeric(event)))
e2 <- factor(d2, labels=c("censor", "pcm", "death w/o pcm", 
                          "death after pcm"))
mfit4 <- survfit(Surv(tstart, tstop, e2) ~ sex, data=newdata, id=id)
plot(mfit2[2,], lty=c(2,1),
     xscale=12, mark.time=FALSE, lwd=2, 
     xlab="Years post diagnosis", ylab="Prevalence")
lines(mfit4[2,3], mark.time=FALSE, xscale=12, col=2, lty=2, lwd=2,
      conf=FALSE)

legend(15, .5, c("male:death w/o pcm", "male: ever pcm", 
                 "male: death after pcm"), col=c(1,1,2), lty=c(1,2,2), 
             lwd=2, bty='n')


###################################################
### code chunk number 9: cfit1
###################################################
mtemp <- mgus2
mtemp$age <- mtemp$age/10   #age in decades (easier coefficients)
mtemp$etime <- etime
mtemp$event <- event

options(show.signif.stars = FALSE)  # display intelligence
cfit2 <- coxph(Surv(futime, death) ~ age + sex + mspike, data=mtemp)
cfit2


###################################################
### code chunk number 10: cfit2
###################################################
cfit1 <- coxph(Surv(ptime, pstat) ~ age + sex + mspike, mtemp)
cfit1
quantile(mgus2$mspike, na.rm=TRUE)


###################################################
### code chunk number 11: mpyears
###################################################
pfit1 <- pyears(Surv(ptime, pstat) ~ sex, mtemp, scale=12)
round(100* pfit1$event/pfit1$pyears, 1)  # PCM rate per year

temp <- summary(mfit1, rmean="common")  #print the mean survival time
round(temp$table[,1:6], 1)


###################################################
### code chunk number 12: mprev
###################################################
tdata <- expand.grid(mspike=c(.5, 1.5), age=c(6,8), sex=c("F", "M"))
surv1 <- survfit(cfit1, newdata=tdata)  # time to progression curves
surv2 <- survfit(cfit2, newdata=tdata)  # time to death curves


###################################################
### code chunk number 13: mprev2
###################################################
cifun <- function(surv1, surv2) {
    utime <- sort(unique(surv1$time, surv2$time))
    jump1 <- diff(c(0, summary(surv1, times=utime, extend=TRUE)$cumhaz))
    jump2 <- diff(c(0, summary(surv2, times=utime, extend=TRUE)$cumhaz))
    dA  <- diag(3)
    prev  <- matrix(0., nrow= 1+length(utime), ncol=3)
    prev[1,1] <- 1  #initial prevalence at time 0: all are in the left box
    for (i in 1:length(utime)) {
        dA[1,2] <- jump1[i]  #fill in the first row of dA(s)
        dA[1,3] <- jump2[i]
        dA[1,1] <- 1- (jump1[i] + jump2[i])
        prev[i+1,] <- prev[i,] %*% dA
    }
    list(time=c(0, utime), P = prev)
}
# Get curves for the 8 cases, and save them in a matrix.
#  Since they all come from the same pair of Cox models, the time values
#  for all curves will be the same
# The cifun function above is only designed to handle one of the 8 covariate
#  patterns at a time, but survival curves can be subscripted.
temp <- cifun(surv1[1], surv2[1])
coxtime <- temp$time
coxdeath <- coxpcm  <- matrix(0., nrow=length(temp$time), ncol=8)
coxdeath[,1] <- temp$P[,3]
coxpcm[,1]   <- temp$P[,2]  
for (i in 2:8){
    temp <- cifun(surv1[i], surv2[i]) 
    coxdeath[,i] <- temp$P[,3]
    coxpcm[,i]   <- temp$P[,2]
}

# Print out a M/F results at 20 years
indx <- match(20*12, coxtime)
progmat <- matrix(coxpcm[indx,], nrow=4)
dimnames(progmat) <- list(c("a=50/ms=0.5", "a=50/ms=1.5", 
                            "a=80/ms=0.5", "a=80/ms=1.5"),
                          c("female", "male"))
round(100*t(progmat), 1)  #males and females at 20 years


###################################################
### code chunk number 14: mprev3
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mfrow=c(1,2))
matplot(coxtime/12, coxpcm[,c(1,3,5,7)], col=c(1,1,2,2),
        lty=c(1,2,1,2), type='l', lwd=2, ylim=range(coxpcm),
        xlab="Years", ylab="Progression to PCM")
legend(1, .23, c("Female: 60", "Male: 60", "Female: 80", "Male: 80"),
       lty=c(1,1,2,2), col=c(1,2,1,2), lwd=2, bty='n')
matplot(coxtime/12, coxpcm[,c(2,4,6,8)], col=c(1,1,2,2),
        lty=c(1,2,1,2), type='l', lwd=2,
        xlab="Years", ylab="Progression to PCM")


###################################################
### code chunk number 15: finegray
###################################################
require(cmprsk)
temp <- mtemp
temp$fstat <- as.numeric(event)  # 1=censor, 2=pcm, 3=death
temp$msex  <- with(temp, 1* (sex=='M'))
fgfit1 <- with(temp, crr(etime, fstat, cov1= cbind(age, msex,  mspike),
                        failcode=2, cencode=1, variance=TRUE))
fgfit2 <- with(temp, crr(etime, fstat, cov1=cbind(age, msex, mspike),
                         failcode=3, cencode=1, variance=TRUE))
cmat <- rbind("FineGray: PCM" = fgfit1$coef,
              "Cox: PCM" = coef(cfit1),
              "FineGray: death" = fgfit2$coef,
              "Cox: death" = coef(cfit2))
round(cmat,2)


###################################################
### code chunk number 16: compare
###################################################
cox.f <- log(1- progmat)    #log(1-P)
round(cox.f[,1] / cox.f[,2], 2)


###################################################
### code chunk number 17: finegray2
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mfrow=c(1,2))
fdata <- model.matrix(~age + sex + mspike, data=tdata)[,-1] #remove intercept
fpred <- predict(fgfit1, cov1=fdata)
matplot(fpred[,1]/12, fpred[,c(2,4,6,8)], col=c(1,1,2,2), lty=c(1,2,1,2),
        ylim=range(fpred[,-1]),
       type='l', lwd=2, xlab="Years", ylab="FG predicted")
legend(0, .22, c("Female, 60", "Male, 60","Female: 80", "Male, 80"),
       col=c(1,2,1,2), lty=c(1,1,2,2), lwd=2, bty='n')
matplot(fpred[,1]/12, fpred[,c(3,5,7,9)], col=c(1,1,2,2), lty=c(1,2,1,2),
       type='l', lwd=2, xlab="Years", ylab="FG predicted")


