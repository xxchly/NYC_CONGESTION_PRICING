# Purpose: paper_tw79 Table 3 data. Cell DiD stays the original (fit_m2_full.R, lme4 REML). Emits the
#   4 gate outcomes per mode: (1) Cell RDST pass /15, (2) Cell mean-test pass /15 [pre-period
#   mean test, Zhang, Sang & Wu (2026) baseline: 2024 dR ~ 1 + (1|season), H0 mean=0], (3) Pooled RDST [intensity
#   / traffic-hour weighted], (4) Pooled mean-test [intensity weighted, 2024]. Self-validates that the
#   refit cell DiD matches stored did_blup__tw79 (cell DiD unchanged).
#   Requires results/did_blup__tw79.csv, i.e. run fit_m2_full.R first.
# Output: results/table3__tw79.csv
suppressMessages({ library(arrow); library(data.table); library(lme4); library(lmerTest) })
# repo root: 1st command-line argument, else env var CCRE_DID_ROOT, else the working directory
ROOT <- c(commandArgs(trailingOnly=TRUE), Sys.getenv("CCRE_DID_ROOT", "."))[1]
PAN  <- file.path(ROOT, "data/panels")
OUT  <- file.path(ROOT, "results")
COV  <- file.path(ROOT, "data/covariates")
MODES <- c("yellow","uber","lyft","subway","citibike")
# adopted 3-variable weather set: panel-level temp/humidity + cell-level gusts
WVARS <- c("temperature_2m","wind_gusts_10m","relative_humidity_2m")
PEAK_ORD <- c("am_peak","off_peak","pm_peak"); HRS <- c(am_peak=3, off_peak=18, pm_peak=3)
TR_FLOWS <- c("0->0","0->1","0->2","1->0","2->0"); CTRL <- c("1->1","1->2","2->1","2->2")
SEAS <- c(winter=1,spring=2,summer=3,fall=4); P <- 4L
zscore <- function(x){ m<-mean(x); s<-sqrt(mean((x-m)^2)); if(s==0) rep(0,length(x)) else (x-m)/s }
CTL <- lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
slope_p <- function(d){ o<-tryCatch(summary(lmer(dR~tau_cal+Post+(1|yg),d,REML=TRUE,control=CTL))$coefficients,error=function(e)NULL)
  if(is.null(o)||!("tau_cal"%in%rownames(o)))return(c(NA,NA)); c(o["tau_cal",1],o["tau_cal",ncol(o)]) }
mean_p  <- function(d){ o<-tryCatch(summary(lmer(dR~1+(1|season),d,REML=TRUE,control=CTL))$coefficients,error=function(e)NULL)
  if(is.null(o))return(c(NA,NA)); c(o["(Intercept)",1],o["(Intercept)",ncol(o)]) }

rs <- fread(file.path(COV,"rainsnow_cov.csv")); gu <- fread(file.path(COV,"gusts_cov.csv"))
stored <- fread(file.path(OUT,"did_blup__tw79.csv")); out <- list(); vmax <- 0
for(mode in MODES){ mm <- mode
  dt <- as.data.table(read_parquet(file.path(PAN, sprintf("ringpair_%s.parquet",mode))))
  dt <- merge(dt, rs[rs[["mode"]]==mm, .(ring_flow, year, season, peak, rain, snowfall)],
              by=c("ring_flow","year","season","peak"), all.x=TRUE)
  dt <- merge(dt, gu[gu[["mode"]]==mm, .(ring_flow, year, season, peak, wind_gusts_10m)],
              by=c("ring_flow","year","season","peak"), all.x=TRUE)
  stopifnot(!anyNA(dt$rain), !anyNA(dt$snowfall), !anyNA(dt$wind_gusts_10m))
  dt[, Y := as.numeric(log1p_avg_daily_rides)]; dt[, treated := as.logical(treated)]
  for(c in WVARS) dt[, (paste0("z_",c)) := zscore(get(c))]
  dt[, rf := as.character(ring_flow)]; dt[, pk := as.character(peak)]
  dt[, `:=`(ring_flow=factor(ring_flow), season=factor(season), dow=factor(dow), peak=factor(peak,levels=PEAK_ORD))]
  dt[, ysdp := interaction(year,season,dow,peak,drop=TRUE)]; dt[, PostTreated := as.integer(post & treated)]
  dt[, fp := interaction(ring_flow,peak,drop=TRUE)]
  zZ <- paste(paste0("z_",WVARS), collapse=" + ")
  m <- lmer(as.formula(paste0("Y ~ PostTreated + (0+PostTreated|fp) + ",zZ," + (1|ring_flow) + (1|ysdp) + ammi")),
            dt, REML=TRUE, control=CTL)
  b2 <- fixef(m)[["PostTreated"]]; re <- ranef(m)$fp[["PostTreated"]]
  bl <- data.table(fp=rownames(ranef(m)$fp), blup=b2+re); bl[, c("ring_flow","peak") := tstrsplit(fp,".",fixed=TRUE)]
  bl <- bl[ring_flow %in% TR_FLOWS]
  chk <- merge(bl, stored[mode==mm,.(ring_flow,peak,blup0=blup)], by=c("ring_flow","peak")); vmax <- max(vmax, max(abs(chk$blup-chk$blup0)))
  dt[, eps := residuals(m)]; dt[, tau_seas := as.numeric(SEAS[as.character(season)])]
  ctrl <- dt[rf %in% CTRL, .(cm=mean(eps)), by=.(year,season,dow,pk)]
  cells <- list(); crp <- 0; cmp <- 0
  for(f in TR_FLOWS) for(p in PEAK_ORD){
    tr <- dt[rf==f & pk==p, .(year,season,dow,tau_seas,eps)]
    d  <- merge(tr, ctrl[pk==p,.(year,season,dow,cm)], by=c("year","season","dow"))
    d[, dR := eps-cm]; d[, Post := as.integer(year==2025L)]; d[, tau_cal := tau_seas+P*Post]; d[, yg := interaction(year,season,drop=TRUE)]
    rp <- slope_p(d); mp <- mean_p(d[year==2024L])
    crp <- crp + (isTRUE(rp[2]>=0.05)); cmp <- cmp + (isTRUE(mp[2]>=0.05))
    d[, `:=`(ring_flow=f, peak=p)]; cells[[length(cells)+1]] <- d
  }
  A <- rbindlist(cells)
  rw <- dt[treated==TRUE & year==2024L, .(rides=mean(avg_daily_rides)), by=.(ring_flow=rf, peak=pk)]
  rw[, w := rides/HRS[peak]]; rw[, w := w/sum(w)]
  Aw <- merge(A, rw[,.(ring_flow,peak,w)], by=c("ring_flow","peak"))
  pl_all <- Aw[, .(dR=weighted.mean(dR,w), tau_seas=tau_seas[1]), by=.(year,season,dow)]
  pl_all[, Post := as.integer(year==2025L)]; pl_all[, tau_cal := tau_seas+P*Post]; pl_all[, yg := interaction(year,season,drop=TRUE)]
  prdst <- slope_p(pl_all); pmean <- mean_p(pl_all[year==2024L])
  out[[mode]] <- data.table(mode=mode,
    pooled_rdst_p=prdst[2], pooled_rdst_pass=isTRUE(prdst[2]>=0.05),
    cell_rdst_pass=crp,
    pooled_mean_p=pmean[2], pooled_mean_pass=isTRUE(pmean[2]>=0.05),
    cell_mean_pass=cmp)
  cat(sprintf("  %-9s RDST[pool p=%.3f %s | cell %d/15]  MEAN[pool p=%.3f %s | cell %d/15]\n",
    mode, prdst[2], ifelse(prdst[2]>=0.05,"pass","FAIL"), crp, pmean[2], ifelse(pmean[2]>=0.05,"pass","FAIL"), cmp))
}
cat(sprintf("cell-DiD refit vs stored: max|diff| = %.2e\n", vmax))
O <- rbindlist(out); fwrite(O, file.path(OUT,"table3__tw79.csv"))
cat(sprintf("DONE: pooled RDST %d/5, cell RDST %d/75, pooled mean %d/5, cell mean %d/75\n",
    sum(O$pooled_rdst_pass), sum(O$cell_rdst_pass), sum(O$pooled_mean_pass), sum(O$cell_mean_pass)))
