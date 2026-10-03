# Purpose: Weather covariate selection sweep for the paper's Variable Selection
#          experiment. Candidate space = ALL subsets of the five candidate weather
#          variables {temperature, rain, snowfall, wind gusts, humidity}
#          (2^5 = 32 sets, incl. the empty set) x 5 modes = 160 M2 fits.
#          M2 structure is held fixed; only the weather terms change.
# Records: pooled DiD est/se/df/p from the REML fit; AIC/BIC from an ML refit
#          (REML likelihoods are NOT comparable across fixed-effect structures,
#          so information criteria for covariate selection must come from ML).
# Output : results/weather_sweep__tw79.csv
suppressMessages({ library(arrow); library(data.table); library(lme4); library(lmerTest) })
# repo root: 1st command-line argument, else env var CCRE_DID_ROOT, else the working directory
ROOT <- c(commandArgs(trailingOnly=TRUE), Sys.getenv("CCRE_DID_ROOT", "."))[1]
PAN  <- file.path(ROOT, "data/panels")
COV  <- file.path(ROOT, "data/covariates")
OUT  <- file.path(ROOT, "results"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
MODES <- c("yellow","uber","lyft","citibike","subway")
VARS  <- c("temperature_2m","rain","snowfall","wind_gusts_10m","relative_humidity_2m")
PEAK_ORD <- c("am_peak","off_peak","pm_peak")
zscore <- function(x){ m<-mean(x); s<-sqrt(mean((x-m)^2)); if(s==0) rep(0,length(x)) else (x-m)/s }
CTL <- lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))

rs <- fread(file.path(COV,"rainsnow_cov.csv")); gu <- fread(file.path(COV,"gusts_cov.csv"))
res <- list()
for(mode in MODES){
  md <- mode
  dt <- as.data.table(read_parquet(file.path(PAN, sprintf("ringpair_%s.parquet",mode))))
  dt <- merge(dt, rs[rs[["mode"]]==md, .(ring_flow, year, season, peak, rain, snowfall)],
              by=c("ring_flow","year","season","peak"), all.x=TRUE)
  dt <- merge(dt, gu[gu[["mode"]]==md, .(ring_flow, year, season, peak, wind_gusts_10m)],
              by=c("ring_flow","year","season","peak"), all.x=TRUE)
  stopifnot(!anyNA(dt$rain), !anyNA(dt$snowfall), !anyNA(dt$wind_gusts_10m))
  dt[, Y := as.numeric(log1p_avg_daily_rides)]; dt[, treated := as.logical(treated)]
  for(c in VARS) dt[, (paste0("z_",c)) := zscore(get(c))]
  dt[, `:=`(ring_flow=factor(ring_flow), season=factor(season), dow=factor(dow),
            peak=factor(peak,levels=PEAK_ORD))]
  dt[, ysdp := interaction(year, season, dow, peak, drop=TRUE)]
  dt[, PostTreated := as.integer(post & treated)]
  dt[, fp := interaction(ring_flow, peak, drop=TRUE)]
  for(mask in 0:31){
    sel <- VARS[bitwAnd(mask, 2^(0:4)) > 0]
    rhs <- if(length(sel)) paste(" +", paste(paste0("z_", sel), collapse=" + ")) else ""
    form <- as.formula(paste0("Y ~ PostTreated", rhs,
              " + ammi + (1|ring_flow) + (1|ysdp) + (0+PostTreated|fp)"))
    m  <- lmer(form, dt, REML=TRUE,  control=CTL)   # est/p: REML (reporting fit)
    mm <- lmer(form, dt, REML=FALSE, control=CTL)   # AIC/BIC: ML (comparable)
    ct <- coef(summary(m))["PostTreated",]
    res[[length(res)+1]] <- data.table(mode=mode, mask=mask,
      vars=if(length(sel)) paste(sel, collapse="+") else "none", n_vars=length(sel),
      est=ct[["Estimate"]], se=ct[["Std. Error"]], df=ct[["df"]], p=ct[["Pr(>|t|)"]],
      AIC=AIC(mm), BIC=BIC(mm), singular=isSingular(m))
  }
  R <- rbindlist(res)[mode==md]
  cat(sprintf("%-9s done: est range [%.4f, %.4f]  p range [%.4f, %.4f]\n",
      md, min(R$est), max(R$est), min(R$p), max(R$p)))
}
R <- rbindlist(res)
fwrite(R, file.path(OUT, "weather_sweep__tw79.csv"))
full <- R[n_vars==5]
cat("\n=== full set (all 5 candidates) ===\n"); print(full[, .(mode, est=round(est,4), p=round(p,4), AIC=round(AIC,1))])
cat("\nDONE: 160 fits ->", file.path(OUT, "weather_sweep__tw79.csv"), "\n")
