# Purpose: Table "M2 reproduces M0 and M1.1.0" on the am7-9/pm17-19 window panel.
#          Per mode, two-sided Welch-Satterthwaite t-test that each of M2's 16 estimates (pooled
#          beta1 + 15 ring x peak flow effects) equals its M0 / M1.1.0 counterpart. Uses the panel's
#          ring x peak AMMI column (the matched AMMI, M1.1.0). df: Satterthwaite (lmerTest) on the
#          M2 side, residual df on the lm side, combined by Welch-Satterthwaite.
# Inputs : data/panels/ringpair_<mode>.parquet, data/covariates/{rainsnow_cov,gusts_cov}.csv
# Outputs: results/integration_pvals__tw79.csv
suppressMessages({library(arrow); library(data.table); library(lme4); library(lmerTest)})
# repo root: 1st command-line argument, else env var CCRE_DID_ROOT, else the working directory
ROOT <- c(commandArgs(trailingOnly=TRUE), Sys.getenv("CCRE_DID_ROOT", "."))[1]
PNL  <- file.path(ROOT, "data/panels")
OUT  <- file.path(ROOT, "results"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
COV  <- file.path(ROOT, "data/covariates")
MODES <- c("yellow","uber","lyft","citibike","subway")
# adopted 3-variable weather set: panel-level temp/humidity + cell-level gusts
WVARS <- c("temperature_2m","wind_gusts_10m","relative_humidity_2m")
FLOWS <- c("0->0","0->1","0->2","1->0","2->0"); PEAKS <- c("am_peak","off_peak","pm_peak")
CELLORD <- as.vector(t(outer(FLOWS, PEAKS, paste, sep="@")))          # flow-major: 15 cells fixed order
zc <- function(x){m<-mean(x); s<-sqrt(mean((x-m)^2)); if(s==0) rep(0,length(x)) else (x-m)/s}
CTL <- lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=3e5))
# Welch-Satterthwaite t-test that two independent estimates are equal:
# df = (sa^2+sb^2)^2 / (sa^4/dfa + sb^4/dfb)  (paper convention: Satterthwaite df)
pt_ws <- function(d, sa, sb, dfa, dfb){
  v <- sa^2 + sb^2
  df <- v^2 / (sa^4/dfa + sb^4/dfb)
  2*pt(-abs(d)/sqrt(v), df)
}
rs <- fread(file.path(COV,"rainsnow_cov.csv")); gu <- fread(file.path(COV,"gusts_cov.csv"))
rows <- list()
for(mode in MODES){
  md <- mode
  dt <- as.data.table(read_parquet(file.path(PNL, sprintf("ringpair_%s.parquet", mode))))
  dt <- merge(dt, rs[rs[["mode"]]==md, .(ring_flow, year, season, peak, rain, snowfall)],
              by=c("ring_flow","year","season","peak"), all.x=TRUE)
  dt <- merge(dt, gu[gu[["mode"]]==md, .(ring_flow, year, season, peak, wind_gusts_10m)],
              by=c("ring_flow","year","season","peak"), all.x=TRUE)
  stopifnot(!anyNA(dt$rain), !anyNA(dt$snowfall), !anyNA(dt$wind_gusts_10m))
  dt[, Y:=as.numeric(log1p_avg_daily_rides)]
  dt[, PostTreated:=as.integer(post & treated)]; dt[, ysdp:=interaction(year,season,dow,peak,drop=TRUE)]
  for(w in WVARS) dt[,(paste0("z_",w)):=zc(get(w))]                    # ammi = panel column (ring x peak)
  dt[, tcell:=factor(ifelse(treated, paste(ring_flow,peak,sep="@"),"ctrl"))]
  zZ <- paste(paste0("z_",WVARS),collapse=" + ")
  m0 <- lm(as.formula(paste0("Y ~ PostTreated + ",zZ," + ammi + ring_flow + ysdp")), dt)
  b0<-coef(m0)[["PostTreated"]]; se0<-summary(m0)$coefficients["PostTreated","Std. Error"]
  m1 <- lm(as.formula(paste0("Y ~ PostTreated:tcell + ",zZ," + ammi + ring_flow + ysdp")), dt)
  s1<-summary(m1)$coefficients; s1<-s1[grepl("^PostTreated:tcell",rownames(s1)),,drop=FALSE]
  rownames(s1)<-sub("PostTreated:tcell","",rownames(s1)); s1<-s1[rownames(s1)!="ctrl",,drop=FALSE]
  b1<-s1[,"Estimate"]; se1<-s1[,"Std. Error"]
  m2 <- lmer(as.formula(paste0("Y ~ PostTreated + ",zZ," + ammi + (1|ring_flow)+(1|ysdp)+(0+PostTreated|tcell)")),
             dt, REML=TRUE, control=CTL)
  b2p<-fixef(m2)[["PostTreated"]]; se2p<-sqrt(vcov(m2)["PostTreated","PostTreated"])
  df2 <- summary(m2)$coefficients["PostTreated","df"]        # Satterthwaite df (lmerTest)
  df0 <- df.residual(m0); df1 <- df.residual(m1)
  re<-ranef(m2,condVar=TRUE)$tcell; pv<-attr(re,"postVar")[1,1,]
  cell2<-b2p+re[["PostTreated"]]; names(cell2)<-rownames(re); se2<-sqrt(pv+se2p^2); names(se2)<-rownames(re)
  p_pool <- pt_ws(b2p-b0, se2p, se0, df2, df0)
  p_cells <- sapply(CELLORD, function(k)
    pt_ws(cell2[k]-b1[k], se2[k], se1[k], df2, df1))
  allp <- c(p_pool, p_cells); pct <- round(mean(allp>=0.05)*100,0)
  r <- as.list(round(c(p_pool, p_cells),3)); names(r) <- c("pooled", CELLORD)
  rows[[mode]] <- data.table(mode=mode, as.data.table(r), pct_notsig=pct)
  cat(sprintf("  %-9s min p=%.3f  %%n.s.=%d\n", mode, min(allp), pct))
}
R <- rbindlist(rows); fwrite(R, file.path(OUT,"integration_pvals__tw79.csv"))
cat(sprintf("DONE: overall min p = %.3f ; modes with 100%% n.s. = %d/5\n",
            min(sapply(rows, function(x) min(unlist(x[, -c("mode","pct_notsig"), with=FALSE])))),
            sum(sapply(rows, function(x) x$pct_notsig==100))))
