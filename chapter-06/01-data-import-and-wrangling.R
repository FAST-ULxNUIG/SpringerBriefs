# Packages: ---------------------------------------------------------------
library(readr) # CRAN v1.3.1
library(fda) # CRAN v5.1.9
library(tidyverse) # CRAN v1.3.0

# Data Import: ------------------------------------------------------------

# From:
# https://springernature.figshare.com/collections/GaitRec_A_large-scale_ground_reaction_force_dataset_of_healthy_and_impaired_gait/4788012/1

# Citation:
############
# Horsak, Brian; Slijepcevic, Djordje; Raberger, Anna-Maria; Schwab, Caterine; Worisch, Marianne; Zeppelzauer, Matthias (2020).
# GaitRec: A large-scale ground reaction force dataset of healthy and impaired gait. figshare. Collection. 
# https://doi.org/10.6084/m9.figshare.c.4788012.v1
############

# Download the following Ground Reaction Forces Data:
# Force (F)
# vertical, anterior-posterior and medio-lateral force components (V, AP, ML)
# Legs: Both Left (left) and Right (right)
# Processing: Smoothed and Normalized (PRO)
# This is a big data set, it might take a while:

# Vertical Force:
GRF_F_V_PRO_left <- read_csv("https://springernature.figshare.com/ndownloader/files/22063191")
GRF_F_V_PRO_right <- read_csv("https://springernature.figshare.com/ndownloader/files/22063119")

# Medio-lateral Force:
GRF_F_ML_PRO_left <- read_csv("https://springernature.figshare.com/ndownloader/files/22063113")
GRF_F_ML_PRO_right <- read_csv("https://springernature.figshare.com/ndownloader/files/22063086")

# Anterior-Posterior Force:
GRF_F_AP_PRO_left <- read_csv("https://springernature.figshare.com/ndownloader/files/22063185")
GRF_F_AP_PRO_right <- read_csv("https://springernature.figshare.com/ndownloader/files/22063101")

# And associated MetaDeta, e.g., session and subject information:
GRF_metadata <- read_csv("https://ndownloader.figshare.com/files/22062960")
#.. you might have to wait.. it Will download!




# Data Wrangling: ---------------------------------------------------------
# Do this manually, possibly could do it not manually.
GRF_F_V_PRO_left$side <- "left"
GRF_F_V_PRO_right$side <- "right"
GRF_F_ML_PRO_left$side <- "left"
GRF_F_ML_PRO_right$side <- "right"
GRF_F_AP_PRO_left$side <- "left"
GRF_F_AP_PRO_right$side <- "right"

GRF_F_V_PRO_left$component <- "vertical"
GRF_F_V_PRO_right$component <- "vertical"
GRF_F_ML_PRO_left$component <- "medio_lateral"
GRF_F_ML_PRO_right$component <- "medio_lateral"
GRF_F_AP_PRO_left$component <- "anterior_posterior"
GRF_F_AP_PRO_right$component <- "anterior_posterior"

# renaming the grf columns:
GRF_F_V_PRO_left <- rename_at(.tbl = GRF_F_V_PRO_left,
          .vars = vars(paste0("F_V_PRO_", 1:101)), 
          .funs = ~ paste0("time_", 0:100))

GRF_F_V_PRO_right <- rename_at(.tbl = GRF_F_V_PRO_right,
                              .vars = vars(paste0("F_V_PRO_", 1:101)), 
                              .funs = ~ paste0("time_", 0:100))

GRF_F_ML_PRO_left <- rename_at(.tbl = GRF_F_ML_PRO_left,
                              .vars = vars(paste0("F_ML_PRO_", 1:101)), 
                              .funs = ~ paste0("time_", 0:100))

GRF_F_ML_PRO_right <- rename_at(.tbl = GRF_F_ML_PRO_right,
                               .vars = vars(paste0("F_ML_PRO_", 1:101)), 
                               .funs = ~ paste0("time_", 0:100))

GRF_F_AP_PRO_left <- rename_at(.tbl = GRF_F_AP_PRO_left,
                              .vars = vars(paste0("F_AP_PRO_", 1:101)), 
                              .funs = ~ paste0("time_", 0:100))

GRF_F_AP_PRO_right <- rename_at(.tbl = GRF_F_AP_PRO_right,
                               .vars = vars(paste0("F_AP_PRO_", 1:101)), 
                               .funs = ~ paste0("time_", 0:100))


stopifnot(all.equal(names(GRF_F_V_PRO_left), names(GRF_F_V_PRO_right)))
stopifnot(all.equal(names(GRF_F_V_PRO_left), names(GRF_F_ML_PRO_left)))
stopifnot(all.equal(names(GRF_F_V_PRO_left), names(GRF_F_ML_PRO_right)))
stopifnot(all.equal(names(GRF_F_V_PRO_left), names(GRF_F_AP_PRO_left)))
stopifnot(all.equal(names(GRF_F_V_PRO_left), names(GRF_F_AP_PRO_right)))

# create full dataset for analysis:

GRF_dataset_PRO <- rbind(
  GRF_F_V_PRO_left,
  GRF_F_V_PRO_right,
  GRF_F_AP_PRO_left,
  GRF_F_AP_PRO_right,
  GRF_F_ML_PRO_left,
  GRF_F_ML_PRO_right
)


inner_join(x = GRF_dataset_PRO, y = GRF_metadata)

GRF_dataset_PRO_meta <- merge.data.frame(
  x = GRF_dataset_PRO, 
  y = GRF_metadata,
  by = c("SUBJECT_ID", "SESSION_ID"))

stopifnot(nrow(GRF_dataset_PRO) == nrow(GRF_dataset_PRO_meta))

GRF_dataset_PRO_meta <- relocate(GRF_dataset_PRO_meta,
                                 time_0:time_100,
                                 .after = last_col())

# For now work only with balanced, controlled subset of the data:
GRF_dataset_PRO_meta <- GRF_dataset_PRO_meta %>%
  filter(SHOD_CONDITION == 1 & SPEED == 2 & TRAIN_BALANCED == 1)

saveRDS(object = GRF_dataset_PRO_meta, 
        file = here::here(here::here(
          "chapter-06",
          "data", 
          "GRF_dataset_PRO_meta.rds")))
