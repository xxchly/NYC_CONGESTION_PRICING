# Companion to weather_sweep__tw79.R: the SAME 2^5 weather-subset sweep run on
# M0 (pooled TWFE, OLS) and M1.1 (interacted TWFE, OLS), so covariate selection
# is validated on the whole model family, not only on M2.
# Records per (model, mode, subset): pooled DiD (M0) or the 15 flow-level
# coefficients (M1.1), plus OLS AIC (ML, valid for fixed-effects comparison).
# Output: results/weather_sweep_m0m1__tw79.csv
suppressMessages({ library(arrow); library(data.table) })
# repo root: 1st command-line argument, else env var CCRE_DID_ROOT, else the working directory
ROOT <- c(commandArgs(trailingOnly=TRUE), Sys.getenv("CCRE_DID_ROOT", "."))[1]
PAN  <- file.path(ROOT, "data/panels")
COV  <- file.path(ROOT, "data/covariates")
OUT  <- file.path(ROOT, "results"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
MODES <- c("yellow","uber","lyft","citibike","subway")
VARS  <- c("temperature_2m","rain","snowfall","wind_gusts_10m","relative_humidity_2m")
zscore <- function(x){ m<-mean(x); s<-sqrt(mean((x-m)^2)); if(s==0) rep(0,length(x)) else (x-m)/s }

rs <- fread(file.path(COV,"rainsnow_cov.csv")); gu <- fread(file.path(COV,"gusts_cov.csv"))
res <- list()
for(mode in MODES){
  md <- mode
  dt <- as.data.table(read_parquet(file.path(PAN, sprintf("ringpair_%s.parquet",mode))))
  dt <- merge(dt, rs[rs[["mode"]]==md, .(ring_flow, year, season, peak, rain, snowfall)],
              by=c("ring_flow","year","season","peak"), all.x=TRUE)
  dt <- merge(dt, gu[gu[["mode"]]==md, .(ring_flow, year, season, peak, wind_gusts_10m)],
              by=c("ring_flow","year","season","peak"), all.x=TRUE)
  dt[, Y := as.numeric(log1p_avg_daily_rides)]
  for(c in VARS) dt[, (paste0("z_",c)) := zscore(get(c))]
  dt[, ysdp := interaction(year, season, dow, peak, drop=TRUE)]
  dt[, PostTreated := as.integer(post & treated)]
  dt[, tcell := factor(ifelse(treated, paste(ring_flow, peak, sep="@"), "ctrl"))]
  for(mask in 0:31){
    sel <- VARS[bitwAnd(mask, 2^(0:4)) > 0]
    rhs <- if(length(sel)) paste(" +", paste(paste0("z_", sel), collapse=" + ")) else ""
    vlab <- if(length(sel)) paste(sel, collapse="+") else "none"
    m0 <- lm(as.formula(paste0("Y ~ PostTreated", rhs, " + ammi + ring_flow + ysdp")), dt)
    res[[length(res)+1]] <- data.table(model="M0", mode=mode, mask=mask, vars=vlab,
      term="pooled", est=coef(m0)[["PostTreated"]], AIC=AIC(m0))
    m1 <- lm(as.formula(paste0("Y ~ PostTreated:tcell", rhs, " + ammi + ring_flow + ysdp")), dt)
    cf <- coef(m1); cf <- cf[grepl("^PostTreated:tcell", names(cf))]
    names(cf) <- sub("PostTreated:tcell", "", names(cf)); cf <- cf[names(cf)!="ctrl"]
    res[[length(res)+1]] <- data.table(model="M1", mode=mode, mask=mask, vars=vlab,
      term=names(cf), est=as.numeric(cf), AIC=AIC(m1))
  }
  cat(mode, "done\n")
}
R <- rbindlist(res)
fwrite(R, file.path(OUT, "weather_sweep_m0m1__tw79.csv"))
cat("DONE:", nrow(R), "rows\n")
