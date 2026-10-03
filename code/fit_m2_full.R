# Purpose: Fit per-mode M2 on the am7-9/pm17-19 window panel and emit BOTH the cell-level (75)
#          and the POOLED (5) DiD coefficients and RDST parallel-trends gates.
#          M2 = (1|ring_flow) + (1|ysdp) + PostTreated (pooled b2) + (0+PostTreated|ring_flow:peak)
#               + 3 z-weather (adopted set: temp + gusts + humidity) + ammi.
#          RDST = residual-difference slope test (see the paper, parallel-trends section).
# Inputs : data/panels/ringpair_<mode>.parquet
#          data/covariates/{rainsnow_cov,gusts_cov}.csv
# Outputs: results/{fixed,did_blup,fitstats,gate_rst_ptr,gate_rst_pooled,varcomp}__tw79.csv
suppressMessages({ library(arrow); library(data.table); library(lme4); library(lmerTest) })
# repo root: 1st command-line argument, else env var CCRE_DID_ROOT, else the working directory
ROOT <- c(commandArgs(trailingOnly=TRUE), Sys.getenv("CCRE_DID_ROOT", "."))[1]
PAN  <- file.path(ROOT, "data/panels")
OUT  <- file.path(ROOT, "results"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
TAG  <- "tw79"
COV  <- file.path(ROOT, "data/covariates")
MODES <- c("yellow","uber","lyft","citibike","subway")
# adopted 3-variable weather set: panel-level temp/humidity + cell-level gusts
WVARS <- c("temperature_2m","wind_gusts_10m","relative_humidity_2m")
PEAK_ORD <- c("am_peak","off_peak","pm_peak")
TR_FLOWS <- c("0->0","0->1","0->2","1->0","2->0"); CTRL <- c("1->1","1->2","2->1","2->2")
SEAS <- c(winter=1,spring=2,summer=3,fall=4); P <- 4L
zscore <- function(x){ m<-mean(x); s<-sqrt(mean((x-m)^2)); if(s==0) rep(0,length(x)) else (x-m)/s }
CTL <- lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
slope_test <- function(form, dat, term){                      # -> c(slope, p) via lmerTest Satterthwaite
  out <- tryCatch(summary(lmer(as.formula(form),dat,REML=TRUE,control=CTL))$coefficients,
                  error=function(e) NULL)
  if(is.null(out) || !(term %in% rownames(out))) return(c(NA,NA)); c(out[term,1], out[term,ncol(out)]) }

rs <- fread(file.path(COV, "rainsnow_cov.csv")); gu <- fread(file.path(COV, "gusts_cov.csv"))
fix<-list(); blup<-list(); stat<-list(); gate<-list(); pool<-list(); vcomp<-list()
for(mode in MODES){
  md <- mode
  dt <- as.data.table(read_parquet(file.path(PAN, sprintf("ringpair_%s.parquet",mode))))
  dt <- merge(dt, rs[rs[["mode"]]==md, .(ring_flow, year, season, peak, rain, snowfall)],
              by=c("ring_flow","year","season","peak"), all.x=TRUE)
  dt <- merge(dt, gu[gu[["mode"]]==md, .(ring_flow, year, season, peak, wind_gusts_10m)],
              by=c("ring_flow","year","season","peak"), all.x=TRUE)
  stopifnot(!anyNA(dt$rain), !anyNA(dt$snowfall), !anyNA(dt$wind_gusts_10m))
  dt[, Y := as.numeric(log1p_avg_daily_rides)]; dt[, treated := as.logical(treated)]
  for(c in WVARS) dt[, (paste0("z_",c)) := zscore(get(c))]
  dt[, rf := as.character(ring_flow)]; dt[, pk := as.character(peak)]
  dt[, `:=`(ring_flow=factor(ring_flow), season=factor(season), dow=factor(dow),
            peak=factor(peak,levels=PEAK_ORD))]
  dt[, ysdp := interaction(year, season, dow, peak, drop=TRUE)]
  dt[, PostTreated := as.integer(post & treated)]
  dt[, fp := interaction(ring_flow, peak, drop=TRUE)]
  zZ <- paste(paste0("z_",WVARS), collapse=" + ")
  form <- paste0("Y ~ PostTreated + (0+PostTreated|fp) + ", zZ,
                 " + (1|ring_flow) + (1|ysdp) + ammi")
  m <- lmer(as.formula(form), dt, REML=TRUE, control=CTL)

  ct <- as.data.frame(coef(summary(m))); ct$term <- rownames(ct); setDT(ct)
  setnames(ct, c("Estimate","Std. Error","df","t value","Pr(>|t|)"),
           c("est","se","df","t","p"), skip_absent=TRUE); ct[, mode := mode]; fix[[mode]] <- ct

  b2 <- fixef(m)[["PostTreated"]]; re <- ranef(m, condVar=TRUE)$fp
  u  <- re[["PostTreated"]]; sd <- sqrt(as.numeric(attr(re,"postVar")[1,1,]))
  bl <- data.table(fp=rownames(re), blup=b2+u, sd=sd)
  bl[, c("ring_flow","peak") := tstrsplit(fp, ".", fixed=TRUE)]
  bl <- bl[ring_flow %in% TR_FLOWS]
  bl[, `:=`(ci_lo=blup-1.96*sd, ci_hi=blup+1.96*sd, mode=mode)]; blup[[mode]] <- bl
  b2se <- as.data.frame(coef(summary(m)))["PostTreated","Std. Error"]
  stat[[mode]] <- data.table(mode=mode, n_obs=nobs(m), logLik=as.numeric(logLik(m)), AIC=AIC(m), BIC=BIC(m),
    sigma=sigma(m), beta2_pooled=b2, b2_lo=b2-1.96*b2se, b2_hi=b2+1.96*b2se,
    converged=is.null(m@optinfo$conv$lme4$code), singular=isSingular(m))
  vc <- as.data.frame(VarCorr(m))
  vcomp[[mode]] <- data.table(mode=mode,
    sigma2_rp    = vc$vcov[vc$grp=="fp"][1],        # flow slope RE variance
    sigma2_r     = vc$vcov[vc$grp=="ring_flow"][1], # ring-pair intercept RE variance
    sigma2_t     = vc$vcov[vc$grp=="ysdp"][1],      # time intercept RE variance
    sigma2_resid = vc$vcov[vc$grp=="Residual"][1])

  # ---- RDST: build the treated-minus-control residual difference series ----
  dt[, eps := residuals(m)]; dt[, tau_seas := as.numeric(SEAS[as.character(season)])]
  ctrl <- dt[rf %in% CTRL, .(cm=mean(eps)), by=.(year, season, dow, pk)]
  poolrows <- list()
  for(f in TR_FLOWS) for(p in PEAK_ORD){
    tr <- dt[rf==f & pk==p, .(year, season, dow, tau_seas, eps)]
    d  <- merge(tr, ctrl[pk==p, .(year, season, dow, cm)], by=c("year","season","dow"))
    d[, dR := eps-cm]; d[, Post := as.integer(year==2025L)]; d[, tau_cal := tau_seas + P*Post]
    d[, yg := interaction(year, season, drop=TRUE)]
    rst <- slope_test("dR ~ tau_cal + Post + (1|yg)", d, "tau_cal")
    gate[[length(gate)+1]] <- data.table(mode=mode, ring_flow=f, peak=p,
      rst_slope=rst[1], rst_p=rst[2], rst_pass=isTRUE(rst[2]>=0.05))
    poolrows[[length(poolrows)+1]] <- d[, .(year, season, dow, tau_seas, dR)]
  }
  # ---- POOLED RDST: average the 15 treated-cell dR at each time point, same slope test ----
  pl <- rbindlist(poolrows)[, .(dR=mean(dR), tau_seas=tau_seas[1]), by=.(year, season, dow)]
  pl[, Post := as.integer(year==2025L)]; pl[, tau_cal := tau_seas + P*Post]
  pl[, yg := interaction(year, season, drop=TRUE)]
  prst <- slope_test("dR ~ tau_cal + Post + (1|yg)", pl, "tau_cal")
  pool[[mode]] <- data.table(mode=mode, beta2_pooled=b2, b2_lo=b2-1.96*b2se, b2_hi=b2+1.96*b2se,
    pooled_slope=prst[1], pooled_p=prst[2], pooled_pass=isTRUE(prst[2]>=0.05))
  cat(sprintf("  %-9s b2=%+.4f pooledRDST p=%.3f pass=%s\n", mode, b2, prst[2], isTRUE(prst[2]>=0.05)))
}
w <- function(x, nm) fwrite(rbindlist(x, fill=TRUE), file.path(OUT, sprintf("%s__%s.csv", nm, TAG)))
w(fix,"fixed"); w(blup,"did_blup"); w(stat,"fitstats"); w(gate,"gate_rst_ptr"); w(pool,"gate_rst_pooled"); w(vcomp,"varcomp")
G <- rbindlist(gate); Pl <- rbindlist(pool)
cat(sprintf("DONE %s: cell RDST %d/%d pass ; pooled RDST %d/5 pass\n",
            TAG, sum(G$rst_pass), nrow(G), sum(Pl$pooled_pass)))
