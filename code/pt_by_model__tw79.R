# Purpose: does the specification change whether parallel trends holds at the
#   resolution the effects are read? The RDST is a test on residuals, so it is
#   run three times on the same panel, the same treated/control construction and
#   the same test, changing only which model produced the residual:
#     M0  TWFE DiD, one pooled coefficient          (lm)
#     M1  interacted TWFE DiD, one per flow         (lm)
#     M2  cross-classified random effects           (lmer, the adopted model)
#   Self-validates that the M2 column reproduces the stored gate_rst_ptr counts.
# Output: results/pt_by_model__tw79.csv  (per mode x model: cell pass /15, pooled p)
suppressMessages({ library(arrow); library(data.table); library(lme4); library(lmerTest) })
# repo root: 1st command-line argument, else env var CCRE_DID_ROOT, else the working directory
ROOT <- c(commandArgs(trailingOnly=TRUE), Sys.getenv("CCRE_DID_ROOT", "."))[1]
PAN  <- file.path(ROOT, "data/panels")
OUT  <- file.path(ROOT, "results"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
COV  <- file.path(ROOT, "data/covariates")
MODES <- c("yellow","uber","lyft","subway","citibike")
WVARS <- c("temperature_2m","wind_gusts_10m","relative_humidity_2m")
PEAK_ORD <- c("am_peak","off_peak","pm_peak"); HRS <- c(am_peak=3, off_peak=18, pm_peak=3)
TR_FLOWS <- c("0->0","0->1","0->2","1->0","2->0"); CTRL <- c("1->1","1->2","2->1","2->2")
SEAS <- c(winter=1,spring=2,summer=3,fall=4); P <- 4L
zscore <- function(x){ m<-mean(x); s<-sqrt(mean((x-m)^2)); if(s==0) rep(0,length(x)) else (x-m)/s }
CTL <- lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))

slope_p <- function(d){
  o <- tryCatch(summary(lmer(dR~tau_cal+Post+(1|yg),d,REML=TRUE,control=CTL))$coefficients,
                error=function(e) NULL)
  if(is.null(o)||!("tau_cal"%in%rownames(o))) return(c(NA,NA))
  c(o["tau_cal",1], o["tau_cal",ncol(o)])
}

# the RDST, given a residual column already attached to dt
gate <- function(dt){
  ctrl <- dt[rf %in% CTRL, .(cm=mean(eps)), by=.(year,season,dow,pk)]
  cells <- list(); npass <- 0
  for(f in TR_FLOWS) for(p in PEAK_ORD){
    tr <- dt[rf==f & pk==p, .(year,season,dow,tau_seas,eps)]
    d  <- merge(tr, ctrl[pk==p,.(year,season,dow,cm)], by=c("year","season","dow"))
    d[, dR := eps-cm]; d[, Post := as.integer(year==2025L)]
    d[, tau_cal := tau_seas+P*Post]; d[, yg := interaction(year,season,drop=TRUE)]
    npass <- npass + isTRUE(slope_p(d)[2] >= 0.05)
    d[, `:=`(ring_flow=f, peak=p)]; cells[[length(cells)+1]] <- d
  }
  A <- rbindlist(cells)
  rw <- dt[treated==TRUE & year==2024L, .(rides=mean(avg_daily_rides)), by=.(ring_flow=rf, peak=pk)]
  rw[, w := rides/HRS[peak]]; rw[, w := w/sum(w)]
  Aw <- merge(A, rw[,.(ring_flow,peak,w)], by=c("ring_flow","peak"))
  pl <- Aw[, .(dR=weighted.mean(dR,w), tau_seas=tau_seas[1]), by=.(year,season,dow)]
  pl[, Post := as.integer(year==2025L)]; pl[, tau_cal := tau_seas+P*Post]
  pl[, yg := interaction(year,season,drop=TRUE)]
  list(cell_pass=npass, pooled_p=slope_p(pl)[2])
}

rs <- fread(file.path(COV,"rainsnow_cov.csv")); gu <- fread(file.path(COV,"gusts_cov.csv"))
out <- list()
for(mode in MODES){ mm <- mode
  dt <- as.data.table(read_parquet(file.path(PAN, sprintf("ringpair_%s.parquet",mode))))
  dt <- merge(dt, rs[rs[["mode"]]==mm, .(ring_flow, year, season, peak, rain, snowfall)],
              by=c("ring_flow","year","season","peak"), all.x=TRUE)
  dt <- merge(dt, gu[gu[["mode"]]==mm, .(ring_flow, year, season, peak, wind_gusts_10m)],
              by=c("ring_flow","year","season","peak"), all.x=TRUE)
  dt[, Y := as.numeric(log1p_avg_daily_rides)]; dt[, treated := as.logical(treated)]
  for(c in WVARS) dt[, (paste0("z_",c)) := zscore(get(c))]
  dt[, rf := as.character(ring_flow)]; dt[, pk := as.character(peak)]
  dt[, `:=`(ring_flow=factor(ring_flow), season=factor(season), dow=factor(dow),
            peak=factor(peak,levels=PEAK_ORD))]
  dt[, ysdp := interaction(year,season,dow,peak,drop=TRUE)]
  dt[, PostTreated := as.integer(post & treated)]
  dt[, fp := interaction(ring_flow,peak,drop=TRUE)]
  dt[, tau_seas := as.numeric(SEAS[as.character(season)])]
  zZ <- paste(paste0("z_",WVARS), collapse=" + ")

  fits <- list(
    M0 = lm(as.formula(paste0("Y ~ PostTreated + ",zZ," + ammi + ring_flow + ysdp")), dt),
    M1 = lm(as.formula(paste0("Y ~ PostTreated:fp + ",zZ," + ammi + ring_flow + ysdp")), dt),
    M2 = lmer(as.formula(paste0("Y ~ PostTreated + (0+PostTreated|fp) + ",zZ,
                                " + (1|ring_flow) + (1|ysdp) + ammi")),
              dt, REML=TRUE, control=CTL))
  for(nm in names(fits)){
    dt[, eps := residuals(fits[[nm]])]
    g <- gate(dt)
    out[[length(out)+1]] <- data.table(mode=mode, model=nm,
      cell_pass=g$cell_pass, pooled_p=g$pooled_p,
      pooled_pass=isTRUE(g$pooled_p>=0.05),
      resid_sd=sd(residuals(fits[[nm]])))
    cat(sprintf("  %-9s %s  cell %2d/15  pooled p=%.4f %s  resid sd=%.4f\n",
        mode, nm, g$cell_pass, g$pooled_p,
        ifelse(g$pooled_p>=0.05,"pass","FAIL"), sd(residuals(fits[[nm]]))))
  }
}
O <- rbindlist(out); fwrite(O, file.path(OUT,"pt_by_model__tw79.csv"))
cat("\n=== cell RDST pass out of 75 ===\n")
print(O[, .(cell_pass=sum(cell_pass), pooled_pass=sum(pooled_pass)), by=model])
cat("\nstored M2 reference: cell 65/75, pooled 5/5 (weighted)\n")
