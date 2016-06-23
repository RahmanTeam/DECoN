### R code from vignette source 'adjcurve.Rnw'

###################################################
### code chunk number 1: adjcurve.Rnw:20-24
###################################################
options(continue="  ", width=60)
options(SweaveHooks=list(fig=function() par(mar=c(4.1, 4.1, .3, 1.1))))
pdf.options(pointsize=8) #text in graph about the same as regular text
if (!exists('coxph')) library(survival)


###################################################
### code chunk number 2: adjcurve.Rnw:133-147
###################################################
group3 <- factor(1+ 1*(flchain$flc.grp >7) + 1*(flchain$flc.grp >9),
                      levels=1:3, 
                      labels=c("FLC < 3.38", "3.38 - 4.71", "FLC > 4.71"))
age1 <- cut(flchain$age, c(49,59,69,79, 110))
levels(age1) <- c(paste(c(50,60,70), c(59,69,79), sep='-'), '80+')
temp1 <- table(group3, age1)
temp2 <- round(100* temp1/rowSums(temp1))
pfun <- function(x,y) {
    paste(ifelse(x<1000, "\\phantom{0}", ""), x, " (", 
          ifelse(y<10,   "\\phantom{0}", ""), y,  ") ", sep="")
}
cat(paste(c("FLC low", pfun(temp1[1,], temp2[1,])), collapse=" & "), "\\\\\n")
cat(paste(c("FLC med", pfun(temp1[2,], temp2[2,])), collapse=" & "), "\\\\\n")
cat(paste(c("FLC high", pfun(temp1[3,], temp2[3,])), collapse=" & "), "\n")


###################################################
### code chunk number 3: flc1
###################################################
getOption("SweaveHooks")[["fig"]]()
fdata <- flchain[flchain$futime > 7,]
fdata$group <- factor(1+ 1*(fdata$flc.grp >7) + 1*(fdata$flc.grp >9),
                      levels=1:3, 
                      labels=c("FLC < 3.38", "3.38 - 4.71", "FLC > 4.71"))
fdata$age2 <- cut(fdata$age, c(49,54, 59,64, 69,74,79, 89, 110),
                  labels = c(paste(c(50,55,60,65,70,75,80),
                                 c(54,59,64,69,74,79,89), sep='-'), "90+"))
		      
sfit1 <- survfit(Surv(futime, death) ~ group, fdata)
plot(sfit1, mark.time=F, col=c(1,2,4), lty=1, lwd=2,
     xscale=365.25, xlab="Years from Sample", 
     ylab="Survival")
text(c(11.1, 10.5, 7.5), c(.88, .57, .4),
     c("FLC < 3.38", "3.38 - 4.71", "FLC > 4.71"), col=c(1,2,4))     


###################################################
### code chunk number 4: flc2
###################################################
getOption("SweaveHooks")[["fig"]]()
temp <- with(fdata, table(group, age2, sex))
size <- apply(temp, 2:3, min)

set.seed(1978)
select <- NULL
dd <- dim(temp)
for (i in 1:dd[1]) {
    for (j in 1:dd[2]) {
        for (k in 1:dd[3]) {
            indx <- which(as.numeric(fdata$group)==i &
                          as.numeric(fdata$age2) ==j &
                          as.numeric(fdata$sex) ==k)
            select <- c(select, sample(indx, size[j,k]))
        }
    }
}

data2 <- fdata[select,]
sfit2 <- survfit(Surv(futime, death) ~ group, data2)
plot(sfit2, mark.time=F, col=c(1,2,4), lty=1, lwd=2,
     xscale=365.25, xlab="Years from Sample", 
     ylab="Survival")
lines(sfit1,  mark.time=F, col=c(1,2,4), lty=2, lwd=1,
      xscale=365.25)
legend(2,.4, levels(fdata$group), lty=1, col=c(1,2,4),
               bty='n', lwd=2)


###################################################
### code chunk number 5: adjcurve.Rnw:244-249
###################################################

tab1 <- with(fdata, table(group, age2, sex))
cat("Low&", paste(tab1[1,,1], collapse=" &"), "\\\\\n")
cat("Med&", paste(tab1[2,,1], collapse=" &"), "\\\\\n")
cat("High&", paste(tab1[3,,1], collapse=" &"), "\\\\\n")


###################################################
### code chunk number 6: adjcurve.Rnw:254-257
###################################################
cat("Low&", paste(tab1[1,,2], collapse=" &"), "\\\\\n")
cat("Med&", paste(tab1[2,,2], collapse=" &"), "\\\\\n")
cat("High&", paste(tab1[3,,2], collapse=" &"), "\n")


###################################################
### code chunk number 7: adjcurve.Rnw:281-292
###################################################
tab3 <- with(fdata, table(age2, group))
tab3 <- round(100*scale(tab3, center=F, scale=colSums(tab3)))
tab4 <- with(data2, table(age2, group))
tab4 <- round(100*scale(tab4, center=F, scale=colSums(tab4)))
tab5 <- cbind(tab3[,1], tab4[,1], tab3[,2], tab4[,2], tab3[,3], tab4[,3])
pfun <- function(x) paste(ifelse(x<10, paste("\\phantom{0}", x), x),
                          collapse=" &")
dtemp <- dimnames(tab5)[[1]]
for (j in 1:7) 
    cat(dtemp[j], " &", pfun(tab5[j,]), "\\\\\n")
cat(dtemp[8], " & ",    pfun(tab5[8,]), "\n")


###################################################
### code chunk number 8: adjcurve.Rnw:323-330
###################################################
# I can't seem to put this all into an Sexpr
z1 <- with(fdata,table(age, sex, group))
z2<- apply(z1, 1:2, min)
ztemp <-  3*sum(z2)
z1b <- with(fdata, table(age>64, sex, group))
ztemp2 <- sum(apply(z1b, 1:2, min))



###################################################
### code chunk number 9: adjcurve.Rnw:359-360
###################################################
survdiff(Surv(futime, death) ~ group, data=data2)


###################################################
### code chunk number 10: adjcurve.Rnw:393-399
###################################################
refpop <- uspop2[as.character(50:100),c("female", "male"), "2000"]
pi.us  <- refpop/sum(refpop)
age100 <- factor(ifelse(fdata$age >100, 100, fdata$age), levels=50:100)
tab100 <- with(fdata, table(age100, sex, group))/ nrow(fdata)
us.wt  <- rep(pi.us, 3)/ tab100  #new weights by age,sex, group
range(us.wt)


###################################################
### code chunk number 11: adjcurve.Rnw:410-419
###################################################
temp <- as.numeric(cut(50:100, c(49, 54, 59, 64, 69, 74, 79, 89, 110)+.5))
pi.us<- tapply(refpop, list(temp[row(refpop)], col(refpop)), sum)/sum(refpop)
tab2 <- with(fdata, table(age2, sex, group))/ nrow(fdata)
us.wt <- rep(pi.us, 3)/ tab2
range(us.wt)
index <- with(fdata, cbind(as.numeric(age2), as.numeric(sex), 
                           as.numeric(group)))
fdata$uswt <- us.wt[index]                      
sfit3a <-survfit(Surv(futime, death) ~ group, data=fdata, weight=uswt) 


###################################################
### code chunk number 12: flc3a
###################################################
getOption("SweaveHooks")[["fig"]]()
tab1 <- with(fdata, table(age2, sex))/ nrow(fdata)
matplot(1:8, cbind(pi.us, tab1), pch="fmfm", col=c(2,2,1,1),
        xlab="Age group", ylab="Fraction of population",
        xaxt='n')
axis(1, 1:8, levels(fdata$age2))

tab2 <- with(fdata, table(age2, sex, group))/nrow(fdata)
tab3 <- with(fdata, table(group)) / nrow(fdata)

rwt <- rep(tab1,3)/tab2 
round(rwt[,1,], 1) #show female data
fdata$rwt <- rwt[index]  # add per subject weights to the data set
sfit3 <- survfit(Surv(futime, death) ~ group, data=fdata, weight=rwt)


###################################################
### code chunk number 13: flc3
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(sfit3, mark.time=F, col=c(1,2,4), lty=1, lwd=2,
     xscale=365.25, xlab="Years from Sample", 
     ylab="Survival")
lines(sfit3a,  mark.time=F, col=c(1,2,4), lty=1, lwd=1,
      xscale=365.25)
lines(sfit1,  mark.time=F, col=c(1,2,4), lty=2, lwd=1,
      xscale=365.25)
legend(2,.4, levels(fdata$group), lty=1, col=c(1,2,4),
               bty='n', lwd=2)


###################################################
### code chunk number 14: adjcurve.Rnw:498-507
###################################################
id <- 1:nrow(fdata)
cfit <- coxph(Surv(futime, death) ~ group + cluster(id), data=fdata, 
              weight=rwt)
summary(cfit)$robscore

if (exists("svykm")) { #true if the survey package is loaded
    sdes <- svydesign(id = ~0, weights=~rwt, data=fdata)
    dfit <- svykm(Surv(futime, death) ~ group, design=sdes, se=TRUE)
}


###################################################
### code chunk number 15: ipw
###################################################
options(na.action="na.exclude")
gg <- as.numeric(fdata$group)
lfit1 <- glm(I(gg==1) ~ factor(age2) * sex, data=fdata,
            family="binomial")
lfit2 <- glm(I(gg==2) ~ factor(age2) * sex, data=fdata,
            family="binomial")
lfit3 <- glm(I(gg==3) ~ factor(age2) * sex, data=fdata,
            family="binomial")
temp <- ifelse(gg==1, predict(lfit1, type='response'),
               ifelse(gg==2, predict(lfit2, type='response'),
                      predict(lfit3, type='response')))
all.equal(1/temp, fdata$rwt)


###################################################
### code chunk number 16: flc4
###################################################
getOption("SweaveHooks")[["fig"]]()
lfit1b <-glm(I(gg==1) ~ age + sex, data=fdata,
            family="binomial")
lfit2b <- glm(I(gg==2) ~ age +sex, data=fdata,
            family="binomial")
lfit3b <- glm(I(gg==3) ~ age + sex, data=fdata,
            family="binomial")

# weights for each group using simple logistic
twt <- ifelse(gg==1, 1/predict(lfit1b, type="response"),
              ifelse(gg==2, 1/predict(lfit2b, type="response"),
                                     1/predict(lfit3b, type="response")))
tdata <- data.frame(fdata, lwt=twt)

#grouped plot for the females
temp <- tdata[tdata$sex=='F',]
temp$gg <- as.numeric(temp$group)
c1 <- with(temp[temp$gg==1,], tapply(lwt, age2, sum))
c2 <- with(temp[temp$gg==2,], tapply(lwt, age2, sum))
c3 <- with(temp[temp$gg==3,], tapply(lwt, age2, sum))

xtemp <- outer(1:8, c(-.1, 0, .1), "+") #avoid overplotting
ytemp <- 100* cbind(c1/sum(c1), c2/sum(c2), c3/sum(c3))

matplot(xtemp, ytemp, col=c(1,2,4),
        xlab="Age group", ylab="Weighted frequency (%)", xaxt='n')
ztab <- table(fdata$age2)
points(1:8, 100*ztab/sum(ztab), pch='+', cex=1.5, lty=2)
# Add the unadjusted
temp <- tab2[,1,]
temp <- scale(temp, center=F, scale=colSums(temp))
matlines(1:8, 100*temp, pch='o', col=c(1,2,4), lty=2)
axis(1, 1:8, levels(fdata$age2))


###################################################
### code chunk number 17: adjcurve.Rnw:638-648
###################################################
# compute new weights
wtscale <- table(fdata$group)/ tapply(fdata$rwt, fdata$group, sum)
wt2 <- c(fdata$rwt * wtscale[fdata$group])
c("rescaled cv"= sd(wt2)/mean(wt2), "rwt cv"=sd(fdata$rwt)/mean(fdata$rwt))

cfit2a <- coxph(Surv(futime, death) ~ group + cluster(id),
                data=fdata, weight= rwt)
cfit2b <- coxph(Surv(futime, death) ~ group + cluster(id),
                data=fdata, weight=wt2) 
round(c(cfit2a$rscore, cfit2b$rscore),1)


###################################################
### code chunk number 18: strata
###################################################
allfit <- survfit(Surv(futime/365.25, death) ~ group + 
                               age2 + sex, fdata)
temp <- summary(allfit)$table
temp[1:6, c(1,4)] #abbrev printout to fit page


###################################################
### code chunk number 19: flc5
###################################################
getOption("SweaveHooks")[["fig"]]()
xtime <- seq(0, 14, length=57)  #four points/year for 14 years
smat <- matrix(0, nrow=57, ncol=3) # survival curves
serr <- smat  #matrix of standard errors
pi <- with(fdata, table(age2, sex))/nrow(fdata)  #overall dist
for (i in 1:3) {
    temp <- allfit[1:16 + (i-1)*16] #curves for group i
    for (j in 1:16) {
        stemp <- summary(temp[j], times=xtime, extend=T)
        smat[,i] <- smat[,i] + pi[j]*stemp$surv
        serr[,i] <- serr[,i] + pi[i]*stemp$std.err^2
        }
    }
serr <- sqrt(serr)

matplot(xtime, smat, type='l', lwd=2, col=c(1,2,4), ylim=c(0,1),
        lty=1, xlab="Years from sample", ylab="Survival")
lines(sfit1, mark.time=F, lty=2, col=c(1,2,4), xscale=365.25)


###################################################
### code chunk number 20: adjcurve.Rnw:773-774
###################################################
survdiff(Surv(futime, death) ~ group + strata(age2, sex), fdata)


###################################################
### code chunk number 21: flc8
###################################################
getOption("SweaveHooks")[["fig"]]()
cfit4a <- coxph(Surv(futime, death) ~ age + sex + strata(group),
               data=fdata)
surv4a <- survfit(cfit4a)
plot(surv4a, col=c(1,2,4), mark.time=F, xscale=365.25,
     xlab="Years post sample", ylab="Survival")


###################################################
### code chunk number 22: flc6
###################################################
getOption("SweaveHooks")[["fig"]]()
tab4a <- with(fdata, table(age, sex))
uage <- as.numeric(dimnames(tab4a)[[1]])
tdata <- data.frame(age = uage[row(tab4a)],
                    sex = c("F","M")[col(tab4a)],
                    count= c(tab4a))
tdata3 <- tdata[rep(1:nrow(tdata), 3),]  #three copies
tdata3$group <- factor(rep(1:3, each=nrow(tdata)), 
                       labels=levels(fdata$group))
sfit4a <- survexp(~group, data=tdata3, weight = count, 
                    ratetable=cfit4a) 

plot(sfit4a, mark.time=F, col=c(1,2,4), lty=1, lwd=2,
     xscale=365.25, xlab="Years from Sample", 
     ylab="Survival")
lines(sfit3,  mark.time=F, col=c(1,2,4), lty=2, lwd=1,
      xscale=365.25)
legend(2,.4, c("FLC low", "FLC med", "FLC high"), lty=1, col=c(1,2,4),
               bty='n', lwd=2)


###################################################
### code chunk number 23: adjcurve.Rnw:885-892
###################################################
tfit <- survfit(cfit4a, newdata=tdata, se.fit=FALSE)
curves <- vector('list', 3)
twt <- c(tab4a)/sum(tab4a)
for (i in 1:3) {
    temp <- tfit[i,]
    curves[[i]] <- list(time=temp$time, surv= c(temp$surv %*% twt))
    }


###################################################
### code chunk number 24: flc6b
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mfrow=c(1,2))
cfit4b <- coxph(Surv(futime, death) ~ age*sex + strata(group),
                fdata)
sfit4b <- survexp(~group, data=tdata3, ratetable=cfit4b, weights=count)
plot(sfit4b, fun='event', xscale=365.25,
     xlab="Years from sample", ylab="Deaths")
lines(sfit3, mark.time=FALSE, fun='event', xscale=365.25, lty=2)
lines(sfit4a, fun='event', xscale=365.25, col=2)

temp <- median(fdata$sample.yr)                                
mrate <- survexp.mn[as.character(uage),, as.character(temp)]
crate <- predict(cfit4b, newdata=tdata, reference='sample', type='lp')
crate <- matrix(crate, ncol=2)[,2:1] # mrate has males then females, match it
# crate contains estimated log(hazards) relative to a baseline,
#  and mrate absolute hazards, make both relative to a 70 year old
for (i in 1:2) {
    mrate[,i] <- log(mrate[,i]/ mrate[21,2])
    crate[,i] <- crate[,i] - crate[21,2]
    }
matplot(mrate, crate, col=2:1, type='l')
abline(0, 1, lty=2, col=4)


###################################################
### code chunk number 25: adjcurve.Rnw:963-971
###################################################
getOption("SweaveHooks")[["fig"]]()
obs <- with(fdata, tapply(death, list(age2, sex, group), sum))
pred<- with(fdata, tapply(predict(cfit4b, type='expected'),
              list(age2, sex, group), sum))
excess <- matrix(obs/pred, nrow=8)  #collapse 3 way array to 2
dimnames(excess) <- list(dimnames(obs)[[1]], c("low F", "low M",
                                               "med F", "med M",
                                               "high F", "high M"))
round(excess, 1)


###################################################
### code chunk number 26: adjcurve.Rnw:987-999
###################################################
cfit5a <- coxph(Surv(futime, death) ~ group:age +sex + 
                strata(group), fdata) 
cfit5b <- coxph(Surv(futime, death) ~ group:(age +sex) + 
                strata(group), fdata)
cfit5c <- coxph(Surv(futime, death) ~ group:(age *sex) + 
                strata(group), fdata)

options(show.signif.stars=FALSE) # see footnote
anova(cfit4a, cfit5a, cfit5b, cfit5c)
temp <- coef(cfit5a)
names(temp) <- c("sex", "ageL", "ageM", "ageH")
round(temp,3)


###################################################
### code chunk number 27: flc7
###################################################
getOption("SweaveHooks")[["fig"]]()
pred5a <- with(fdata, tapply(predict(cfit5a, type='expected'),
              list(age2, sex, group), sum))
excess5a <- matrix(obs/pred5a, nrow=8,
                   dimnames=dimnames(excess))
round(excess5a, 1)

sfit5  <- survexp(~group,  data=tdata3, ratetable=cfit5a, weights=count)
plot(sfit3, fun='event', xscale=365.25, mark.time=FALSE, lty=2,
     xlab="Years from sample", ylab="Deaths")
lines(sfit5, fun='event', xscale=365.25)


