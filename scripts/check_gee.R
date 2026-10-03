library(dplyr)
library(geepack)
e <- parse('application.R')
for (x in e) {
  if (is.call(x) && identical(x[[1]], as.name('<-')) &&
      identical(x[[2]], as.name('gee_summary'))) eval(x)
}
set.seed(21)
d <- data.frame(race = rep(1:30, sample(4:12, 30, TRUE)))
d$sex <- sample(c('Female', 'Male'), nrow(d), TRUE)
d$value <- rnorm(nrow(d)) + rep(rnorm(30), table(d$race)) +
  0.4*(d$sex == 'Female')
for (paired in c(TRUE, FALSE)) {
 d$female <- as.numeric(d$sex == 'Female')
 f <- if (paired) value ~ 1 else value ~ female
 full <- geeglm(f, id = race, data = d, family = gaussian(),
   corstr = 'independence', scale.fix = TRUE)
 g <- gee_summary(d$value, d$race, if (paired) NULL else d$sex)
 j <- if (paired) 1 else 2
 stopifnot(abs(g$mean_difference - coef(full)[j]) < 1e-10,
   abs(g$std_error - sqrt(full$geese$vbeta[j,j])) < 1e-10)
 perm <- sample(nrow(d))
 h <- gee_summary(d$value[perm], d$race[perm],
   if (paired) NULL else d$sex[perm])
 stopifnot(abs(h$std_error-g$std_error) < 1e-10)
}
cat('GEE aggregation and row-order checks passed\n')
