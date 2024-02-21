# Libraries --------------------------------------------------------------------
library(here)

# Data processing --------------------------------------------------------------

source(here("scripts", "processing_scripts", "01-import_spring_data.R"))
source(here("scripts", "processing_scripts", "02-import_fall_data.R"))
source(here("scripts", "processing_scripts", "03-get_plant_status.R"))
source(here("scripts", "processing_scripts", "04-calc_props.R"))

# Data analysis ----------------------------------------------------------------

source(here("scripts", "analysis_scripts", "01-glm.R"))
source(here("scripts", "analysis_scripts", "02-model_props.R"))
source(here("scripts", "analysis_scripts", "03-create_figure_1_map.R"))
source(here("scripts", "analysis_scripts", "04-glm_reco_reviewer.R"))
source(here("scripts", "analysis_scripts", "04-rda_reco_reviewer2.R"))
source(here("scripts", "functions.R"))
