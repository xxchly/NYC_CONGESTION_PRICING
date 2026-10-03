# E3 strips + E2 dow grid, on the M2 side:
#  (a) AMMI variant comparison on M2: fit M2 with each AMMI variant and save
#      the 15 flow-level BLUPs per variant in DiD % (strip for Fig 3 right).
#  (b) Grid decomposition on M2 with day-of-week: random DiD slope on the FULL
#      ring x season3 x peak x dow grid (225 cells), save the 225 BLUPs.
# Outputs: results/ammi_necessity_strip_m2__tw79.csv,
#          results/decomp_m2_blups_dow__tw79.csv
suppressMessages({ library(arrow); library(data.table); library(lme4) })
# repo root: 1st command-line argument, else env var CCRE_DID_ROOT, else the working directory
ROOT <- c(commandArgs(trailingOnly=TRUE), Sys.getenv("CCRE_DID_ROOT", "."))[1]
PAN  <- file.path(ROOT, "data/panels")
COV  <- file.path(ROOT, "data/covariates")
OUT  <- file.path(ROOT, "results"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
MODES <- c("yellow","uber","lyft","citibike","subway")
WVARS <- c("temperature_2m","wind_gusts_10m","relative_humidity_2m")
PEAK_ORD <- c("am_peak","off_peak","pm_peak")
TR_FLOWS <- c("0->0","0->1","0->2","1->0","2->0")
zscore <- function(x){ m<-mean(x); s<-sqrt(mean((x-m)^2)); if(s==0) rep(0,length(x)) else (x-m)/s }
CTL <- lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))

rank1_ammi <- function(pre, dt, rowcol, colcols){
  kr_p <- as.character(pre[[rowcol]])
  kc_p <- do.call(paste, c(lapply(colcols, function(cc) as.character(pre[[cc]])), sep="|"))
  kr_d <- as.character(dt[[rowcol]])
  kc_d <- do.call(paste, c(lapply(colcols, function(cc) as.character(dt[[cc]])), sep="|"))
  cm <- tapply(pre$Y, list(kr_p, kc_p), mean)
  cm[is.na(cm)] <- mean(cm, na.rm=TRUE)
  Mc <- cm - rowMeans(cm) %o% rep(1, ncol(cm)) - rep(1, nrow(cm)) %o% colMeans(cm) + mean(cm)
  sv <- svd(Mc); g <- sv$u[,1]*sqrt(sv$d[1]); h <- sv$v[,1]*sqrt(sv$d[1])
  if (h[which.max(abs(h))] < 0){ g <- -g; h <- -h }
  gmap <- setNames(g, rownames(cm)); hmap <- setNames(h, colnames(cm))
  as.numeric(gmap[kr_d] * hmap[kc_d])
}

strip <- list(); dec <- list()
for(mode in MODES){
  md <- mode
  dt <- as.data.table(read_parquet(file.path(PAN, sprintf("ringpair_%s.parquet",mode))))
  gu <- fread(file.path(COV,"gusts_cov.csv"))
  dt <- merge(dt, gu[gu[["mode"]]==md, .(ring_flow, year, season, peak, wind_gusts_10m)],
              by=c("ring_flow","year","season","peak"), all.x=TRUE)
  dt[, Y := as.numeric(log1p_avg_daily_rides)]
  for(c in WVARS) dt[, (paste0("z_",c)) := zscore(get(c))]
  dt[, season3 := ifelse(season %in% c("spring","fall"), "rest", as.character(season))]
  dt[, `:=`(ring_flow=factor(ring_flow), season=factor(season), dow=factor(dow),
            peak=factor(peak,levels=PEAK_ORD))]
  dt[, ysdp := interaction(year, season, dow, peak, drop=TRUE)]
  dt[, PostTreated := as.integer(post & treated)]
  dt[, fp := interaction(ring_flow, peak, drop=TRUE)]
  pre <- dt[year==2024L]
  zZ <- paste(paste0("z_",WVARS), collapse=" + ")

  # ---- (a) AMMI variants on M2: save the 15 flow-level BLUPs per variant ----
  VAR <- list(none=NULL,
              ring_season=list("ring_flow", "season3"),
              ring_dow=list("ring_flow", "dow"),
              peak_season=list("peak", "season3"),
              ring_peak=list("ring_flow", "peak"))
  for(v in names(VAR)){
    if(is.null(VAR[[v]])) dt[, av := 0] else
      dt[, av := rank1_ammi(pre, dt, VAR[[v]][[1]], VAR[[v]][2])]
    aterm <- if(v=="none") "" else " + av"
    m <- lmer(as.formula(paste0("Y ~ PostTreated + (0+PostTreated|fp) + ", zZ,
              aterm, " + (1|ring_flow) + (1|ysdp)")), dt, REML=TRUE, control=CTL)
    b2 <- fixef(m)[["PostTreated"]]; re <- ranef(m)$fp
    bl <- data.table(fp=rownames(re), blup=b2+re[["PostTreated"]])
    bl[, c("ring_flow","peak") := tstrsplit(fp, ".", fixed=TRUE)]
    bl <- bl[ring_flow %in% TR_FLOWS]
    strip[[length(strip)+1]] <- bl[, .(mode=mode, variant=v, ring_flow, peak,
                                       pct=100*(exp(blup)-1))]
  }

  # ---- (b) full-grid with dow (225 cells) random slope on M2 ----
  dt[, gcell := factor(ifelse(dt$treated,
      paste(as.character(ring_flow), season3, as.character(peak),
            as.character(dow), sep="@"), "ctrl"))]
  dt[, av := rank1_ammi(pre, dt, "ring_flow", list("season3","peak","dow"))]  # composite nuisance
  m <- lmer(as.formula(paste0("Y ~ PostTreated + (0+PostTreated|gcell) + ", zZ,
            " + av + (1|ring_flow) + (1|ysdp)")), dt, REML=TRUE, control=CTL)
  b2 <- fixef(m)[["PostTreated"]]; re <- ranef(m)$gcell
  bl <- data.table(cell=rownames(re), blup=b2+re[["PostTreated"]])
  bl <- bl[cell!="ctrl"]
  bl[, c("ring_flow","season3","peak","dow") := tstrsplit(cell, "@", fixed=TRUE)]
  dec[[length(dec)+1]] <- bl[, .(mode=mode, ring_flow, season3, peak, dow, blup)]
  cat(mode, "done\n")
}
fwrite(rbindlist(strip), file.path(OUT, "ammi_necessity_strip_m2__tw79.csv"))
fwrite(rbindlist(dec), file.path(OUT, "decomp_m2_blups_dow__tw79.csv"))
cat("DONE\n")
