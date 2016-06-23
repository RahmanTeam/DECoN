### R code from vignette source 'timedep.Rnw'

###################################################
### code chunk number 1: preamble
###################################################
options(width=60, continue=" ")
makefig <- function(file, top=1, right=1, left=4) {
    pdf(file, width=9.5, height=7, pointsize=18)
    par(mar=c(4, left, top, right) +.1)
    }
library(survival)


###################################################
### code chunk number 2: timedep.Rnw:75-77 (eval = FALSE)
###################################################
## fit <- coxph(Surv(time1, time2, status) ~ age + creatinine, 
##              data=mydata)


###################################################
### code chunk number 3: timedep.Rnw:155-158
###################################################
cgd[1:10, c("id", "tstart", "tstop", "status", "enum", "treat")]
cfit <- coxph(Surv(tstart, tstop, status) ~ treat + sex + age +
               inherit + cluster(id), data=cgd)


###################################################
### code chunk number 4: timedep.Rnw:228-236
###################################################
load('raheart.rda')
age2 <- tcut(raheart$agechf*365.25, 0:110* 365.25, labels=0:109)
rowid <- 1:nrow(raheart)
pfit <- pyears(Surv(startday, stopday, hospevt) ~ age2 + rowid,
               data=raheart, data.frame=TRUE, scale=1)
print(pfit$offtable)
pdata <- pfit$data
print(pdata[1:6,])


###################################################
### code chunk number 5: timedep.Rnw:265-282
###################################################
index <- as.integer(pdata$rowid)
lagtime <- c(0, pdata$pyears[-nrow(pdata)])
lagtime[1+ which(diff(index)==0)] <- 0 #starts at 0 for each subject
temp <- raheart$startday[index] + lagtime  #start of each new interval
data2 <- data.frame(raheart[index,], 
                    time1= temp,
                    time2= temp + pdata$pyears,
                    event= pdata$event,
                    age2=  1+ as.numeric(pdata$age2) )

afit1 <- coxph(Surv(startday, stopday, hospevt) ~ male + pspline(agechf), 
               data=raheart)
afit2 <- coxph(Surv(time1, time2, event) ~ male + pspline(age2), data2)
#termplot(afit1, terms=2, se=TRUE, xlab="Age at Diagnosis of CHF")
#termplot(afit2, terms=2, se=TRUE, xlab="Current Age")

table(with(raheart, tapply(hospevt, patid, sum)))


###################################################
### code chunk number 6: timedep.Rnw:298-302
###################################################
afit2b <- coxph(Surv(startday, stopday, hospevt) ~ male + tt(agechf),
                data=raheart, 
                tt=function(x, t, ...) pspline(x + t/365.25))
afit2b


###################################################
### code chunk number 7: timedep.Rnw:328-335
###################################################
function(x, t, riskset, weights){ 
    obrien <- function(x) {
        r <- rank(x)
        (r-.5)/(.5+length(r)-r)
    }
    unlist(tapply(x, riskset, obrien))
}


###################################################
### code chunk number 8: timedep.Rnw:345-347
###################################################
function(x, t, riskset, weights) 
    unlist(tapply(x, riskset, rank))


