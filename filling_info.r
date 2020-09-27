validate <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
    }), warning=function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  list(value=value, warning=warn, error=err)
}

fill_info <- function(inpdt){
# Type 1
  inpdt$tx.sd <- with(inpdt, replmiss(tx.sd, (tx.SEmean * sqrt(tx.n))))
  inpdt$ctl.sd <- with(inpdt, replmiss(ctl.sd, (ctl.SEmean * sqrt(ctl.n))))

  inpdt <- escalc(data=inpdt,measure="SMD",
                  m1i=tx.mean,sd1i=tx.sd,n1i=tx.n,
                  m2i=ctl.mean,sd2i=ctl.sd,n2i=ctl.n,
                  var.names = c("delta", "v"))
# Type 2
  inpdt$g <- with(inpdt,t.IG*sqrt((tx.n + ctl.n)/(tx.n*ctl.n)))
  inpdt$bias.corr <- with(inpdt, (1 - (3 / (4*(tx.n + ctl.n - 2) - 1))))
# Type 3
  inpdt$t.IG <- with(inpdt, replmiss(t.IG,
                            (sign) * qt( (p.twotail/2),(tx.n + ctl.n - 2),
                                         lower.tail = FALSE)))

# Type 4
  inpdt$t.IG <-
    with(inpdt, replmiss(t.IG,
                (sign)*qt(p.t.onetail, (tx.n + ctl.n - 2),lower.tail = FALSE)))

# Type 5
  inpdt$t.IG <- with(inpdt,replmiss(t.IG,(sign)*(sqrt(F.1dfB))))

# Type 6
  inpdt$F.1dfB <-
    with(inpdt, replmiss(F.1dfB,
                qf(p.F1dfB, df1 = 1,df2 = (tx.n + ctl.n - 2),lower.tail = FALSE)))

  inpdt$t.IG <-
    with(inpdt,replmiss(t.IG,(sign)*(sqrt(F.1dfB))))
# Type 7 - 8
  inpdt$F.2dfB <-
    with(inpdt,replmiss(F.2dfB,
               qf(p.F2dfB, df1=2, df2=(tx.n + ctl.n + oth.n - 3),
                    lower.tail = FALSE)))
  inpdt$p.F2dfB
  grandmean <- with(inpdt,
       (((tx.n*tx.mean) + (oth.n*oth.mean) + (ctl.n*ctl.mean)) / (tx.n + oth.n + ctl.n)))

  MSB <- with(inpdt,
       ((tx.n)*((tx.mean - grandmean)^2) + (oth.n)*((oth.mean - grandmean)^2) + (ctl.n)*((ctl.mean - grandmean)^2) ) / 2)

  MSW <- with(inpdt, (MSB / F.2dfB))
# Fill g, delta, v
  inpdt$g <- with(inpdt, replmiss(g, t.IG*sqrt((tx.n + ctl.n)/(tx.n*ctl.n))))
  inpdt$g <- with(inpdt, replmiss(g, ((tx.mean-ctl.mean)/sqrt(MSW))))

  inpdt$bias.corr <- with(inpdt, replmiss(bias.corr, 1 - (3 / (4*(inpdt$tx.n + inpdt$ctl.n - 2) - 1))))
  inpdt$delta <- with(inpdt, replmiss(delta, g*bias.corr))
  inpdt$v <- with(inpdt,replmiss(v,
                         ((tx.n+ctl.n) / (tx.n*ctl.n)+((delta^2)/(2*(tx.n+ctl.n))))))
#
  return(inpdt)
}
