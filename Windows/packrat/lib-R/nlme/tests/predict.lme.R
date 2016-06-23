## from PR#8905
library(nlme)
fm <- lme(distance ~ poly(age, 3) + Sex, data = Orthodont, random = ~ 1)
# data for predictions
Newdata <- head(Orthodont)
Newdata$Sex <- factor(Newdata$Sex, levels = levels(Orthodont$Sex))
(pr <- predict(fm, Newdata))
stopifnot(all.equal(c(pr), fitted(fm)[1:6]))

## https://stat.ethz.ch/pipermail/r-devel/2013-September/067600.html
## but with a different fix.

m0 <- lme(distance ~ Sex, random = ~1|Subject, data = Orthodont)
Fitted <- predict(m0, level = 0)
Fitted.Newdata <- predict(m0, level = 0, newdata = Orthodont)
stopifnot(sum(abs(Fitted - Fitted.Newdata)) == 0)

Fitted <- predict(m0, level = 1)
Fitted.Newdata <- predict(m0, level = 1, newdata = Orthodont)
sum(abs(Fitted - Fitted.Newdata))
stopifnot(sum(abs(Fitted - Fitted.Newdata)) == 0)

m1 <- lme(distance ~ 1, random = ~1|Subject, data = Orthodont)
Fitted <- predict(m1, level = 0)
Fitted.Newdata <- predict(m1, level = 0, newdata = Orthodont)
stopifnot(sum(abs(Fitted - Fitted.Newdata)) == 0)

Fitted <- predict(m1, level = 1)
Fitted.Newdata <- predict(m1, level = 1, newdata = Orthodont)
stopifnot(sum(abs(Fitted - Fitted.Newdata)) == 0)

m2 <- lme(distance ~ 0, random = ~1|Subject, data = Orthodont)
Fitted <- predict(m2, level = 0)
Fitted.Newdata <- predict(m2, level = 0, newdata = Orthodont)
stopifnot(sum(abs(Fitted - Fitted.Newdata)) == 0)

Fitted <- predict(m2, level = 1)
Fitted.Newdata <- predict(m2, level = 1, newdata = Orthodont)
stopifnot(sum(abs(Fitted - Fitted.Newdata)) == 0)
