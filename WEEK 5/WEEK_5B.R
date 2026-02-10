# ============================================================
# Exercise 5.7 (Second): NEON–AmeriFlux co-location & FC/H comparison
# - Find co-located sites with metScanR (fallback to known pair if needed)
# - Download AmeriFlux BASE & NEON dp04 (DP4.00200.001)
# - Align timestamps (AMF local standard time) and compare FC/H vs 1:1
# - Print RMSE, Bias, r; draw scatter plots
# ============================================================

# -------------------------
# 0) Packages & setup
# -------------------------
base_pkgs <- c("dplyr","ggplot2","lubridate","stringr")
for (p in base_pkgs) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
lapply(base_pkgs, library, character.only = TRUE)

if (!requireNamespace("amerifluxr", quietly=TRUE)) install.packages("amerifluxr")
library(amerifluxr)   # amf_site_info(), amf_download_base(), amf_read_base

if (!requireNamespace("neonUtilities", quietly=TRUE)) install.packages("neonUtilities")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
if (!requireNamespace("rhdf5", quietly=TRUE)) BiocManager::install("rhdf5")
library(neonUtilities) # zipsByProduct(), stackEddy()

# Try installing metScanR from GitHub ZIP (avoids some API 404 issues)
get_metScanR <- function() {
  if (requireNamespace("metScanR", quietly=TRUE)) return(TRUE)
  if (!requireNamespace("remotes", quietly=TRUE)) install.packages("remotes")
  zip_url <- "https://github.com/jaroberti/metScanR/archive/refs/heads/master.zip"
  tf <- tempfile(fileext = ".zip")
  ok <- try({
    utils::download.file(zip_url, tf, mode="wb", quiet=TRUE)
    install.packages(tf, repos=NULL, type="source")
  }, silent=TRUE)
  requireNamespace("metScanR", quietly=TRUE)
}
has_metScanR <- isTRUE(get_metScanR())
if (has_metScanR) library(metScanR)

# -------------------------
# 1) Pick region & find co-located pair (metScanR → fallback)
# -------------------------
# Example search near Harvard Forest (you can change lat/lon)
lat0 <- 42.53
lon0 <- -72.18

if (has_metScanR) {
  ms <- metScanR::siteFinder(Lat = lat0, Lon = lon0)       # Nearby stations across networks
  nearby <- ms$site.data %>%
    arrange(distance_km) %>%
    filter(network %in% c("NEON","AMERIFLUX"))
  message("Nearest NEON & AmeriFlux by metScanR (top 10):")
  print(head(nearby, 10))
  
  neon_pick <- nearby %>% filter(network=="NEON") %>% slice(1)
  amf_pick  <- nearby %>% filter(network=="AMERIFLUX") %>% slice(1)
  
  neon_site_raw <- as.character(neon_pick$station_id)
  amf_site      <- as.character(amf_pick$station_id)
  
  # Extract 4-letter NEON code if station_id includes extra text
  neon_site <- stringr::str_match(neon_site_raw, "([A-Z]{4})$")[,2]
  if (is.na(neon_site)) neon_site <- neon_site_raw
} else {
  message("metScanR not available; using known co-located pair: US-xHA (AmeriFlux) ↔ HARV (NEON).")
  amf_site  <- "US-xHA"
  neon_site <- "HARV"
}

message(sprintf("Chosen pair → AmeriFlux: %s ; NEON: %s", amf_site, neon_site))

# -------------------------
# 2) Define overlapping time window
# -------------------------
# Choose a recent ~6-month window; adjust as needed
t1 <- Sys.Date() - 1
t0 <- t1 - 180

# -------------------------
# 3) AmeriFlux BASE download & subset (FC/H)
# -------------------------
# Supply your AmeriFlux credentials
user  <- "Ogonna"
email <- "Ogonnaeli@yahoo.com"
if (!nzchar(user) || !nzchar(email)) stop("Please set your AmeriFlux user/email.")

amf_files <- amf_download_base(
  user_id = user, user_email = email,
  site_id = amf_site,
  data_product = "BASE-BADM",
  data_policy = "CCBY4.0",
  agree_policy = TRUE,
  intended_use = "education",
  intended_use_text = "Exercise 5.7 NEON–AmeriFlux FC/H comparison",
  out_dir = tempdir(), verbose = TRUE
)
amf_df <- amf_read_base(amf_files)

# BASE timestamps are local standard time (no DST). Keep FC & H in [t0, t1]
amf_sub <- amf_df %>%
  mutate(tstart_local = lubridate::ymd_hm(as.character(TIMESTAMP_START))) %>%
  filter(tstart_local >= t0, tstart_local <= t1) %>%
  select(tstart_local, FC_amf = FC, H_amf = H)

# -------------------------
# 4) NEON dp04 (DP4.00200.001) download & extract
# -------------------------
save_dir <- file.path(tempdir(), paste0("neon_dp04_", neon_site))
dir.create(save_dir, showWarnings = FALSE)

start_month <- format(t0, "%Y-%m")
end_month   <- format(t1, "%Y-%m")

# Try released files; if none, retry with include.provisional=TRUE
ok <- TRUE
tryCatch({
  zipsByProduct(
    dpID="DP4.00200.001", site=neon_site,
    startdate=start_month, enddate=end_month,
    package="basic", savepath=save_dir,
    check.size=FALSE, include.provisional=FALSE
  )
}, error=function(e){ message("Released dp04 not available; retrying with include.provisional=TRUE ..."); ok <<- FALSE })

if (!ok) {
  zipsByProduct(
    dpID="DP4.00200.001", site=neon_site,
    startdate=start_month, enddate=end_month,
    package="basic", savepath=save_dir,
    check.size=FALSE, include.provisional=TRUE
  )
}

flux_list <- stackEddy(file.path(save_dir, "filesToStack00200"), level="dp04", metadata=TRUE)
neon_dp04 <- flux_list[[neon_site]]

# -------------------------
# 5) Convert NEON UTC → AmeriFlux local standard time & select FC/H
# -------------------------
# Get UTC offset from AmeriFlux site metadata (if available), else fallback (e.g., EST = -5)
amf_sites <- amf_site_info()
utc_offset <- suppressWarnings(as.numeric(amf_sites$UTC_OFFSET[amf_sites$SITE_ID == amf_site]))
if (is.na(utc_offset)) {
  message("UTC offset missing; defaulting to -5 (EST). Adjust if needed.")
  utc_offset <- -5
}

# Map NEON dp04 turbulent fluxes to AMF variables: FC (CO2), H (sensible heat)
neon_sub <- neon_dp04 %>%
  transmute(
    # NEON times are UTC → convert to local STANDARD time with fixed offset (no DST)
    tstart_local = ymd_hm(format(timeBgn + lubridate::hours(utc_offset), "%Y%m%d%H%M")),
    FC_neon = `data.fluxCo2.turb.flux`,  # µmol CO2 m-2 s-1
    H_neon  = `data.fluxTemp.turb.flux`  # W m-2
  ) %>%
  filter(tstart_local >= t0, tstart_local <= t1)

# -------------------------
# 6) Join by time & compare
# -------------------------
cmp <- inner_join(amf_sub, neon_sub, by="tstart_local")

rmse <- function(a,b) sqrt(mean((a-b)^2, na.rm=TRUE))
stats <- tibble::tibble(
  metric = c("RMSE_FC","Bias_FC","r_FC","RMSE_H","Bias_H","r_H"),
  value  = c(
    rmse(cmp$FC_amf, cmp$FC_neon),
    mean(cmp$FC_neon - cmp$FC_amf, na.rm=TRUE),
    suppressWarnings(cor(cmp$FC_amf, cmp$FC_neon, use="pairwise.complete.obs")),
    rmse(cmp$H_amf, cmp$H_neon),
    mean(cmp$H_neon - cmp$H_amf, na.rm=TRUE),
    suppressWarnings(cor(cmp$H_amf, cmp$H_neon, use="pairwise.complete.obs"))
  )
)
print(stats)

# -------------------------
# 7) Plots: FC and H scatter vs 1:1
# -------------------------
p_fc <- ggplot(cmp, aes(FC_amf, FC_neon)) +
  geom_point(alpha=0.4, color="#1b9e77") +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="gray40") +
  labs(x="AmeriFlux FC (µmol m⁻² s⁻¹)", y="NEON dp04 FC (µmol m⁻² s⁻¹)",
       title=paste0(amf_site, " — FC: AmeriFlux vs NEON (", t0, " to ", t1, ")"))

p_h <- ggplot(cmp, aes(H_amf, H_neon)) +
  geom_point(alpha=0.4, color="#d95f02") +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="gray40") +
  labs(x="AmeriFlux H (W m⁻²)", y="NEON dp04 H (W m⁻²)",
       title=paste0(amf_site, " — H: AmeriFlux vs NEON (", t0, " to ", t1, ")"))

print(p_fc); print(p_h)

# Optionally save outputs
# ggsave(filename=file.path(getwd(),"fc_scatter_USxHA_vs_NEON.png"), plot=p_fc, width=7, height=5, dpi=150)
# ggsave(filename=file.path(getwd(),"h_scatter_USxHA_vs_NEON.png"),  plot=p_h, width=7, height=5, dpi=150)

message("Done. If alignment looks off, verify UTC offset and the overlap window, or widen the time window.")
