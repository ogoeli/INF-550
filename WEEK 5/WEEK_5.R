install.packages("neonUtilities")
# install packages if needed
install.packages(c("amerifluxr", "neonUtilities", "tidyverse", "lubridate"))
if(!require(devtools)){install.packages("devtools")}
devtools::install_github("chuhousen/amerifluxr")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("rhdf5", quietly = TRUE)) BiocManager::install("rhdf5")


# load neonUtilities
library(neonUtilities)
library("amerifluxr")
library(neonUtilities)
library(tidyverse)
library(lubridate)

#NEON TOKEN SOURCE FILE
source("C:\\Users\\Ope4\\OneDrive - Northern Arizona University\\Desktop\\ACADEMIC_SEMSTER\\SPRING 2026\\INF_550\\neon_token_source.R")

foliar <- loadByProduct(dpID="DP1.10026.001", site="all", 
                        package="expanded", check.size=F,
                        token=NEON_TOKEN,
                        savepath = "C:/Users/Ope4/NEON_DATA")

# See AmeriFlux sites; NEON-in-AmeriFlux sites start with 'US-x'
amf_sites <- amf_site_info()

#	US-xHA, NEON Harvard Forest (HARV), USA,MA,DBF,2017,NA,https://ameriflux.lbl.gov/sites/siteinfo/US-xHA
# US-xHA started in 2017 and ended in 2024
# SELECT ALL SITES IN USA
neon_amf <- dplyr::filter(amf_sites, stringr::str_detect(SITE_ID, "^US-x"))
neon_amf %>% dplyr::select(SITE_ID, SITE_NAME, IGBP)


# WEEK 5A
# Site choices for this run
amf_site  <- "US-xHA"   # AmeriFlux ID for Harvard Forest
neon_site <- "HARV"     # NEON code mapped to US-xHA  (official mapping)           [1](https://www.neonscience.org/impact/observatory-blog/gap-filled-and-partitioned-neon-data-products-available-part-ameriflux)
site_tz   <- "America/New_York"
utc_offset_hours <- -5  # Eastern Standard Time offset for AmeriFlux "local standard time" (no DST).  [3](https://ameriflux.lbl.gov/data/uploading-half-hourly-hourly-data/)

# =======================
# 1) Download AmeriFlux BASE for US-xHA
# =======================
user  <- "Ogonna"
email <- "Ogonnaeli@yahoo.com"
password <- "Miracle2009+++"

base_paths <- amf_download_base(
  user_id = user, user_email = email,
  site_id = amf_site, data_product = "BASE-BADM",
  data_policy = "CCBY4.0", agree_policy = TRUE,
  intended_use = "education",
  intended_use_text = "Extending US-xHA with NEON dp04 for class exercise",
  out_dir = tempdir(), verbose = TRUE
) # API and function docs: amerifluxr  [6](https://search.r-project.org/CRAN/refmans/amerifluxr/html/amf_download_base.html)

amf_base <- amf_read_base(base_paths) %>%          # parses BASE file to data frame  [7](https://www.rdocumentation.org/packages/amerifluxr/versions/1.0.0)
  select(TIMESTAMP_START, TIMESTAMP_END, FC, H, LE)

# Convert AmeriFlux local standard times -> UTC to know where to extend.
# (UTC = local_standard + |utc_offset_hours| hours; for EST, +5h)
amf_last_end_utc <- ymd_hm(as.character(max(amf_base$TIMESTAMP_END)), tz = "UTC") +
  hours(+abs(utc_offset_hours))

# =======================
# 2) Download NEON dp04 HDF5 *after* AmeriFlux coverage
# =======================
start_month <- format(amf_last_end_utc + days(1), "%Y-%m")
end_month   <- format(Sys.Date(), "%Y-%m")

save_dir <- file.path(tempdir(), "neon_dp04_usxha")
dir.create(save_dir, showWarnings = FALSE)

zipsByProduct(                                     # NEON Data API bundle downloader  [2](https://www.neonscience.org/resources/learning-hub/tutorials/eddy-data-intro)
  dpID = "DP4.00200.001", site = neon_site,
  startdate = start_month, enddate = end_month,
  package = "basic", savepath = save_dir, check.size = FALSE
)

flux_list <- stackEddy(
  filepath = file.path(save_dir, "filesToStack00200"),
  level    = "dp04", metadata = TRUE
) # returns list: [[HARV]], variables, objDesc         [2](https://www.neonscience.org/resources/learning-hub/tutorials/eddy-data-intro)

neon_dp04 <- flux_list[[neon_site]]

# ---- 3) Standardize AmeriFlux timestamps to character (12-digit) ----
amf_base <- amf_base %>%
  mutate(
    TIMESTAMP_START = sprintf("%012.0f", as.numeric(TIMESTAMP_START)),
    TIMESTAMP_END   = sprintf("%012.0f", as.numeric(TIMESTAMP_END))
  )
# (AmeriFlux uses YYYYMMDDHHMM strings in local standard time)  #  [1](https://ameriflux.lbl.gov/data/uploading-half-hourly-hourly-data/)

# ---- 4) Rebuild NEON → AmeriFlux timestamps with *fixed* UTC offset (no DST) ----
# NEON dp04 times are UTC; AmeriFlux requires local standard time (EST = UTC - 5h).
utc_offset_hours <- -5  # US-xHA (Harvard Forest) standard offset
neon_after <- neon_dp04 %>%
  filter(timeEnd > amf_last_end_utc) %>%         # this 'amf_last_end_utc' you already created
  transmute(
    TIMESTAMP_START = format(timeBgn + hours(utc_offset_hours), "%Y%m%d%H%M"),
    TIMESTAMP_END   = format(timeEnd + hours(utc_offset_hours), "%Y%m%d%H%M"),
    FC = `data.fluxCo2.turb.flux`,               # µmol CO2 m-2 s-1 (AmeriFlux FC)
    H  = `data.fluxTemp.turb.flux`,              # W m-2 (AmeriFlux H)
    LE = `data.fluxH2o.turb.flux`                # W m-2 (AmeriFlux LE)
  ) %>%
  arrange(TIMESTAMP_START)
# (NEON dp04 variable mapping → FC/H/LE; units align with AmeriFlux)  #  [2](https://www.neonscience.org/resources/learning-hub/tutorials/eddy-data-intro)[3](https://ameriflux.lbl.gov/data/aboutdata/data-variables/)

# ---- 5) Bind rows now that both sides are character timestamps ----
extended <- bind_rows(
  amf_base %>% select(TIMESTAMP_START, TIMESTAMP_END, FC, H, LE),
  neon_after
) %>% arrange(TIMESTAMP_START)

# ---- 6) Plot (parse timestamps for display only) ----
# These timestamps are in *local standard time*; parse as such for plotting.
site_tz <- "America/New_York"
p_fc <- ggplot(extended, aes(x = ymd_hm(TIMESTAMP_START, tz = site_tz), y = FC)) +
  geom_line(color="#1b9e77", linewidth=0.3) +
  labs(x=NULL, y="FC (µmol m⁻² s⁻¹)", title=paste0("US-xHA — Extended FC time series"))

p_h <- ggplot(extended, aes(x = ymd_hm(TIMESTAMP_START, tz = site_tz), y = H)) +
  geom_line(color="#d95f02", linewidth=0.3) +
  labs(x=NULL, y="H (W m⁻²)", title=paste0("US-xHA — Extended H time series"))

p_fc; p_h
