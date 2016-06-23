### R code from vignette source 'splines.Rnw'

###################################################
### code chunk number 1: splines.Rnw:20-24
###################################################
options(continue="  ", width=60)
options(SweaveHooks=list(fig=function() par(mar=c(4.1, 4.1, .3, 1.1))))
pdf.options(pointsize=8) #text in graph about the same as regular text
options(contrasts=c("contr.treatment", "contr.poly")) #reset default


###################################################
### code chunk number 2: fit1
###################################################
getOption("SweaveHooks")[["fig"]]()
require(survival)
options(show.signif.stars=FALSE) # display intelligence
fit1 <- coxph(Surv(futime, death) ~ sex + pspline(age, 3), data=flchain)
fit1
termplot(fit1, term=2, se=TRUE, col.term=1, col.se=1,
         ylab="log hazard")


###################################################
### code chunk number 3: fit2
###################################################
agem <- with(flchain, ifelse(sex=="M", age, 60))
agef <- with(flchain, ifelse(sex=="F", age, 60))
fit2 <- coxph(Surv(futime, death) ~ sex + pspline(agef, df=3)
              + pspline(agem, df=3), data=flchain)
anova(fit2, fit1)


###################################################
### code chunk number 4: plot2
###################################################
getOption("SweaveHooks")[["fig"]]()
# predictions
pterm <- termplot(fit2, term=2:3, se=TRUE, plot=FALSE)
# reference
refdata <- data.frame(sex=c('F', 'M'), agef=c(65, 60), agem=c(60,65))
pred.ref <- predict(fit2, newdata=refdata, type="lp")
# females
tempf <- pterm$agef$y + outer(pterm$agef$se, c(0, -1.96, 1.96))
frow <- which(pterm$agef$x == 65)
tempf <- tempf  - tempf[frow,1]  # shift curves
# males
tempm <- pterm$agem$y + outer(pterm$agem$se, c(0, -1.96, 1.96))
mrow  <- which(pterm$agem$x == 65)
tempm <- tempm + diff(pred.ref) - tempm[mrow,1]
# plot
matplot(pterm$agef$x, exp(tempf), log='y', col=1, 
        lty=c(1,2,2), type='l', lwd=c(2,1,1),
        xlab="Age", ylab="Relative risk of death")
matlines(pterm$agem$x, exp(tempm), log='y', 
         col=2, lwd=c(2,1,1), lty=c(1,2,2))
legend(80, 1, c("Female", "Male"), lty=1, lwd=2, col=1:2, bty='n')


