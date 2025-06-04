## MONTHLY VEGETATION + SHEEP-OFFTAKE MODEL (Armstrong et al. adapted)
###Modified to intergrate with st-sim
### adding in young and mature trees
### Jo Furtado
#June 2025

#set-up
library(dplyr)
library(terra)
library(rsyncrosim) #package for working with syncrosim
 ##random change

#Heather_grass <- ssimLibrary("C:/Local data/Modelling/Heather_grasssimplev3.ssim")
#Heather_grass_grazing <- rsyncrosim::project(Heather_grass, project = "Phase change, proportions, external program")
#Heather_scenario <- rsyncrosim::scenario(Heather_grass_grazing, scenario ="New R script runs")
#inputRaster2 <- datasheetSpatRaster(Heather_scenario, datasheet="stsim_OutputSpatialState", timestep=1)
#resultsscenario <- rsyncrosim::run(Heather_scenario)
#resultsscenario <- scenario(Heather_grass_grazing, scenario = 41)
#datasheet(resultsscenario, summary=TRUE)
#datasheet(resultsscenario, name = 'stsim_OutputStratum')
#datasheet(resultsscenario, name = 'stsim_StateAttributeValue')
#datasheet(resultsscenario, name = 'stsim_OutputSpatialStateAttribute')
#datasheet(resultsscenario, name = 'stsim_OutputSpatialState')
#datasheet(resultsscenario, name = 'stsim_OutputStock')

# Set working directory from SyncroSim environment
e <- ssimEnvironment()
workingDir <- e$TempDirectory

# ensure dirs exist then move into workingDir so every relative path below
# resolves there; keeps temp and archive separate
if (!dir.exists(workingDir)) dir.create(workingDir, recursive = TRUE)
setwd(workingDir)

## Set up environment variables
timestep <- e$BeforeTimestep #retrieves value of current timestep from e
prevTimestep <- timestep - 1
iteration <- e$BeforeIteration

it_tag <- sprintf("it%03d", iteration)
ts_tag <- sprintf("ts%03d", timestep)

#it_tag <- 1
#ts_tag <- 1

Spatial_state <- "stsim_OutputSpatialState" #habitat raster from stsm NB
#SpatialNPP <- 'stsim_OutputSpatialStateAttribute' #NPP and biomass
Grazing_levels <- 'stsim_TransitionSpatialMultiplier'#name of grazing rasters to be passed back to syncrosim
#Standing_biomass <- 'stsim_FlowSpatialMultiplier'#name of biomass rasters to be passed back to st-sim, flowGroupID= 'Grazing biomass removed'
#name of transition groups to hand back to syncrosim
TransitionGroupid_nograzing <- "No grazing"
TransitionGroupid_vlightgrazing <- "Very light grazing [Type]"
TransitionGroupid_lightgrazing <- "Light grazing [Type]"
TransitionGroupid_lightmodgrazing <- "Light-moderate grazing [Type]"
TransitionGroupid_moderategrazing <- "Moderate grazing [Type]"
TransitionGroupid_modheavygrazing <- "Moderate-heavy grazing [Type]"
TransitionGroupid_heavygrazing <- "Heavy grazing [Type]"
TransitionGroupid_heavyvheavygrazing <- "Heavy to very heavy grazing [Type]"
TransitionGroupid_vheavygrazing <- "Very heavy grazing [Type]"
#FlowGroupID <- "Grazing biomass removed"

## Set path to store modified TIF file once analysis is complete
tempFileName <- stringr::str_c("inputRaster.it",
                               iteration, ".ts",
                               timestep, ".tif")

tempFilePath <- file.path(workingDir, tempFileName)

## Load Library, Result Scenario, and Parent Scenario
myLibrary <- ssimLibrary()
myScenario <- scenario()
myScenarioParentId <- parentId(myScenario) 
myParentScenario <- scenario(myLibrary, scenario=myScenarioParentId)

## Load habitat raster output from previous timestep and Result Scenario
inputRaster <- rsyncrosim::datasheetSpatRaster(myScenario,
                                               datasheet=Spatial_state,
                                               iteration=iteration,
                                               timestep=prevTimestep)

writeRaster(inputRaster, "C:/Local data/Modelling/input_raster.tif", overwrite=TRUE)

#read in files - NB these can/should all be read in from St-SIM later on, reading straight in now whilst developing script
#inputRaster <- terra::rast("C:/Local data/Modelling/Pratice landscapes/Heather wood grass/heatherwoodgrass_habitats_simp.tif") #habitat IDs
management_rast <- terra::rast("C:/Local data/Modelling/Pratice landscapes/Heather wood grass/management_model.tif") #management rasters: for hefts and tree-planting exclosures
ecozone_rast <- terra::rast("C:/Local data/Modelling/Pratice landscapes/Heather wood grass/ecozone.tif") #tells you ecozones
NPP <- rast("C:/Local data/Modelling/Pratice landscapes/Heather wood grass/NPP_hwg2.tif")
# #this tells me annual DM production according to habitat and ecozone. NB: IN g/10x10m2 cell. Based on figures from Armstrong et al (1997a)
biomass  <- rast("C:/Local data/Modelling/Pratice landscapes/Heather wood grass/Initialbiomass_variablehwg.tif") #read in existing standing biomass for a cell: NPP plus existing biomass

## Defining some constants used below - can be changed/added to depending on the stocking calendar used, when more habitats are added in etc.
days_in_month  <- c(31,28,31,30,31,30,31,31,30,31,30,31)
graze_months   <- 5:9                    # May–Sept currently only months grazed in model
utilisation_cap<- 0.60                  # 60 % of standing pool can be removed
digest_grass   <- 0.60                  # U4 grass digestibility (fraction) - from Armstrong et al 
digest_heather <- 0.35                  # heather digestibility (fraction) - from Armstrong et al
W_kg           <- 40                    # average Herwick ewe mass - change if grazing with Swaledales too (significantly heavier)
W075           <- W_kg^0.75             # metabolic weight for average Herdwick
ID_grass       <- function(D) 166.6*D - 43.6      # Armstrong eqn for daily dry matter intake from grass per ewe (g kg^-0.75 d^-1)
phys_cap_g_day <- ID_grass(digest_grass) * W075   # Total potential daily intake from grass, again from Armstrong et al  (g DM day^-1)
ddm_req_g_day  <- phys_cap_g_day * digest_grass   # Daily digestible DM requirement for maintenance of ewe (g/day)

# raster‑safe minima helpers
rmin2 <- function(a, b)            terra::app(c(a, b), fun = min)
rmin3 <- function(a, b, c)         terra::app(c(a, b, c), fun = min)

# Division that works for raster / scalar & raster / raster
safe_div <- function(num, den){
  if(inherits(den, "SpatRaster") && inherits(num, "SpatRaster"))
    return(terra::ifel(den == 0, Inf, num / den))
  if(inherits(num, "SpatRaster") && is.numeric(den)){
    if(den == 0) return(num * Inf)
    return(num / den)
  }
  if(is.numeric(num) && inherits(den, "SpatRaster")){
    num <- terra::rast(num, extent(den), nlyrs = 1)
    return(terra::ifel(den == 0, Inf, num / den))
  }
  if(den == 0) return(Inf)
  num / den
}

#-------finding monthly and annual NPP production for vegetation at different altitudes, splitting into monthly layers--------

#defining monthly fractions of annual NPP for each habitat type, as set out in Armstrong et al
frac_heather <- c(0,0,0,0, 0.081,0.325,0.532,0.062, 0,0,0,0)
frac_u4grass <- c(0.013,0.015,0.037,0.067,0.115,0.135,
                  0.288,0.148,0.068,0.060,0.034,0.020)

#masks raster for each habitat
mask_u4  <- inputRaster == 3      # TRUE where habitat == 3
mask_heath<- inputRaster == 14     # TRUE where habitat == 14
mask_ywood <- inputRaster == 12
mask_mwood <- inputRaster == 10

mask_heathlike <- mask_heath | mask_ywood | mask_mwood

#loop over 12 months, where each habitat type is masked, multiplied by correct habitat fraction, then rasters are written to disk/saved
monthly_rasters <- vector("list", 12)

for (m in 1:12) {
  
  ## A raster of fractions for this month, assembled cell-wise: heather and grass. 
  ## Cells that were neither heather nor U4 are 0; keeps zeros for U4 & heather cells but still masks non-habitat
  frac_rast <-  mask_heath * frac_heather[m] +
  mask_u4   * frac_u4grass[m] +
  mask_ywood   * frac_heather[m] +
  mask_mwood   * (frac_heather[m] * 0.8)  
  
  frac_rast <- ifel(mask_u4 | mask_heathlike, frac_rast, NA)
  
  ## Convert annual NPP (g/cell) to this month’s growth
  npp_month <- NPP * frac_rast
  names(npp_month) <- sprintf("NPP_%02d_%s_%s", m, it_tag, ts_tag)
  
  ## Save 
  writeRaster(npp_month,
              filename = paste0("NPP_", sprintf("%02d_%s_%s", m, it_tag, ts_tag), ".tif"),
              overwrite = TRUE)
  monthly_rasters[[m]] <- npp_month
}

## combine the 12 layers into one multi-band file
npp_monthly_stack <- rast(monthly_rasters) 
npp_stack_file  <- sprintf("NPP_monthly_stack_%s_%s.tif", it_tag, ts_tag)
writeRaster(npp_monthly_stack, npp_stack_file, overwrite = TRUE)

## Convert list → list of SpatRasters for fast indexed access in the loop coming next!
npp_rast <- lapply(sprintf("NPP_%02d_%s_%s.tif", 1:12, it_tag, ts_tag), rast)

#-------Sheep density raster

#Reclassify habitats according to grazing preference (according to distribution of preferred forage, grasses)
grazing_preference <- classify(inputRaster, rcl = matrix(c(
  3, 10,  # grassland
  14,  3, #heathland, assuming all mature
  12, 3, #young woodland assumed to be the same as heathland
  10, 2 # 
), ncol=2, byrow=TRUE))

## Local sum of attractiveness (within 500 m radius, 5 cells)
focal_weights <- terra::focalMat(grazing_preference, 50, type = "circle")
local_pressure <- focal(grazing_preference, w = focal_weights, fun = sum, na.policy="omit", na.rm=TRUE)

# Assign overall density per area: 1 ewe/ha = 0.01 per 100m² cell, 0.6 ewes/ha = 0.006 per cell
density_rcl <- matrix(c(
  1, 0.01,   # Area A
  2, 0.006   # Area B
), ncol = 2, byrow = TRUE)

density_per_cell <- classify(management_rast, rcl = density_rcl)

# Initialise empty raster for sheep density, ensuring it is same size as grazing preference, set values to NA
sheep_density <- grazing_preference
sheep_density[] <- NA

#Loop that builds a raster of final sheep density that scales by local attraction to match zone totals, e.g., sheep are spread out according to forage preference but also makes sure that total stocking density remains as set
#NB: zones refer to management areas here, e.g., hefts, areas of afforestation
zones <- c(1, 2)

for (z in zones) {
  # Mask to this zone
  zone_mask <- management_rast == z #binary mask to select only cells belonging to zone z
  
  # Local pressure and stocking density rasters in this zone
  zone_pressure <- mask(local_pressure, zone_mask, maskvalues = FALSE) #the summed attractiveness of the surrounding 50m radius
  zone_density  <- mask(density_per_cell, zone_mask, maskvalues = FALSE) #the intended per cell stocking rate for each area
  
  # Compute expected total sheep in this zone based on stocking density
  zone_density_vals <- values(zone_density, mat = FALSE) #extract all non-NA values
  expected_total_sheep <- sum(zone_density_vals, na.rm = TRUE) #sum to check totoal sheep in area, should be as expected
  
  # Multiply local pressure by per-cell stocking density (not normalized yet)
  zone_sheep_density <- zone_pressure * zone_density #this multiplies grazing preference per cell by stocking density
  
  # Compute actual (raw) total from the above
  actual_total_sheep <- global(zone_sheep_density, fun = "sum", na.rm = TRUE)[["sum"]] #sums ratser values to see total sheep numbers from this
  
  # Scale back to match expected total sheep
  scale_factor <- expected_total_sheep / actual_total_sheep #scales if there is a mismatch between actual_total_sheep and expected_total sheep
  zone_sheep_density <- zone_sheep_density * scale_factor
  
  # Combine with main raster
  sheep_density <- ifel(is.na(sheep_density), zone_sheep_density, sheep_density) #adds calculted sheep densities per cell for zone Z to empty raster
}

writeRaster(sheep_density, "sheep_density.tif", overwrite=TRUE)

## ---------- OFF-TAKE + MONTHLY GRAZING LOOP -------------------------------
## Three live pools (0,1,2 months) which sheep can graze before growth becomes “standing”.      
## Eat fresh grass first, then *forced* 5 % of demand must come from heather, then older live pools, finally capped standing pool.
## All demand / capacity is cell-wise, scaled by sheep_density.    
## --------------------------------------------------------------------------
## digestibility map (0.60 on grass, 0.35 on heather)
digest_rast <- mask_u4 * digest_grass + 
  mask_heathlike * digest_heather

# live-growth buffers (0-, 1-, 2-month): basically keeping track of NPP stocks available for grazing, assuming only three month's worth of NPP are 'live'
live0 <- biomass * 0 # this-month growth
live1 <- biomass * 0  # 1-month old
live2 <- biomass * 0 # 2-month old

util_list <- vector("list", length(graze_months))
names(util_list) <- sprintf("%02d", graze_months)

offtake_log <- data.frame(month = integer(),
                          ddm_demand_kg = numeric(),
                          ddm_eaten_kg  = numeric(),
                          dm_eaten_grass_kg = numeric(),
                          dm_eaten_heat_kg  = numeric(),
                          mean_util_pct = numeric())

for (m in graze_months) {
  
  days <- days_in_month[m]
  
  ## 1.  add this-month NPP into the live0 'stock'
  npp_month <- npp_rast[[m]]
  live0     <- npp_month
  
  ## snapshot of *all* live growth before grazing (for utilisation map)
  live_before <- live0 + live1 + live2
  
  ## 2. per-cell demand / capacity
  ddm_need <- sheep_density * ddm_req_g_day  * days #total dry matter requirement for maintence of flock
  phys_cap <- sheep_density * phys_cap_g_day * days #most that can be eaten (physical cap)
  ddm_need[is.na(ddm_need)] <- 0
  phys_cap[is.na(phys_cap)] <- 0
  
  ## helper to graze any pool
  eat_pool <- function(pool, cap_lim = 1) {
    avail <- if (cap_lim < 1) pool * cap_lim else pool
    take  <- rmin3(avail,
                   safe_div(ddm_need, digest_rast),
                   phys_cap)
    ddm_need <<- ddm_need - take * digest_rast
    phys_cap <<- phys_cap - take
    take
  }
  
  ## 3. grass from live0 -----------------------------------------------------
  take0 <- eat_pool(live0) ; live0 <- live0 - take0
  
  ## 4. forced 5 % heather share --------------------------------------------
  forced_frac <- 0.05
  need_DDM_H  <- ddm_need * forced_frac
  need_DM_H   <- safe_div(need_DDM_H, digest_heather)
  
  takeH_forced <- rmin3(live0 * mask_heathlike, need_DM_H, phys_cap)
  
  live0    <- live0 - takeH_forced
  ddm_need <- ddm_need - takeH_forced * digest_heather
  phys_cap <- phys_cap - takeH_forced
  
  consumed_new <- take0 + takeH_forced          # running total
  
  ## 5. older live pools -----------------------------------------------------
  take1 <- eat_pool(live1) ; live1 <- live1 - take1
  take2 <- eat_pool(live2) ; live2 <- live2 - take2
  consumed_new <- consumed_new + take1 + take2  # final live tally
  
  ## 6. standing pool (60 % cap) --------------------------------------------
  stand_cap <- biomass * utilisation_cap
  takeS <- rmin3(stand_cap,
                 safe_div(ddm_need, digest_rast),
                 phys_cap)
  biomass <- biomass - takeS
  consumed_stand <- takeS
  
  ## 7. ageing: 2-month shoots drop into standing ---------------------------
  biomass <- biomass + live2
  live2 <- live1
  live1 <- live0
  live0 <- live0 * 0
  
  ## 8. utilisation of live growth ------------------------------------------
  util_rast <- consumed_new / live_before
  util_rast[live_before == 0] <- NA
  
  ## 9. write rasters --------------------------------------------------------
  terra::writeRaster(biomass,
                     sprintf("biomass_after_%02d_%s_%s.tif", m, it_tag, ts_tag),
                     overwrite = TRUE)
  terra::writeRaster(util_rast,
                     sprintf("utilisation_%02d_%s_%s.tif", m, it_tag, ts_tag),
                     overwrite = TRUE)
  util_list[[sprintf("%02d", m)]] <- util_rast
  
  ## 10. bookkeeping ---------------------------------------------------------
  dm_grass <- terra::global((consumed_new + consumed_stand) * mask_u4,
                            "sum", na.rm = TRUE)[1,1]
  dm_heat  <- terra::global((consumed_new + consumed_stand) * mask_heathlike,
                            "sum", na.rm = TRUE)[1,1]
  
  offtake_log <- rbind(offtake_log,
                       data.frame(month = m,
                                  ddm_demand_kg = terra::global(sheep_density *
                                                                  ddm_req_g_day * days,
                                                                "sum", na.rm = TRUE)[1,1] / 1000,
                                  ddm_eaten_kg  = terra::global((consumed_new + consumed_stand) *
                                                                  digest_rast, "sum",
                                                                na.rm = TRUE)[1,1] / 1000,
                                  dm_eaten_grass_kg = dm_grass / 1000,
                                  dm_eaten_heat_kg  = dm_heat  / 1000,
                                  mean_util_pct = terra::global(util_rast,
                                                                "mean", na.rm = TRUE)[1,1] * 100))
}

## 11. stack utilisation layers once ----------------------------------------
util_stack <- terra::rast(util_list)
util_stack_file <- sprintf("utilisation_stack_%s_%s.tif", it_tag, ts_tag)
terra::writeRaster(util_stack, util_stack_file, overwrite = TRUE)

print(offtake_log)

###############################################################################
##  ANNUAL GRAZING PRESSURE  -----------------------------------------------
##  1. collapse monthly utilisation stack  →  single annual proportion
##  2. re-classify into 9 grazing-level bands
##  3. write one binary raster per band and register them with SyncroSim
###############################################################################


## --- folders --------------------------------------------------------------
#dir.create(workingDir, showWarnings = FALSE, recursive = TRUE)
#setwd(workingDir)                                         # all files land here

#archiveDir <- "C:/Local data/Modelling/grazing_binaries"  # ↩ your browse-able folder
#ir.create(archiveDir, showWarnings = FALSE, recursive = TRUE)

## --- labels & transition-group IDs ---------------------------------------
labels <- c("no_grazing", "very_light_grazing", "light_grazing",
            "light_moderate_grazing", "moderate_grazing",
            "moderate_heavy_grazing", "heavy_grazing",
            "heavy_very_heavy_grazing", "very_heavy_grazing")

transition_group_ids <- c(
  TransitionGroupid_nograzing,
  TransitionGroupid_vlightgrazing,
  TransitionGroupid_lightgrazing,
  TransitionGroupid_lightmodgrazing,
  TransitionGroupid_moderategrazing,
  TransitionGroupid_modheavygrazing,
  TransitionGroupid_heavygrazing,
  TransitionGroupid_heavyvheavygrazing,
  TransitionGroupid_vheavygrazing
)



## --------------------------------------------------------------------------
## 1.  annual utilisation (May–Sep summed)
## --------------------------------------------------------------------------
annual_util <- app(util_stack, fun = sum, na.rm = TRUE)
annual_util <- terra::ifel(annual_util > 1, 1, annual_util)  #cap so it cannot exceed 1
names(annual_util) <- "annual_util"

## --------------------------------------------------------------------------
## 2.  classify into grazing bands
## --------------------------------------------------------------------------
breaks  <- c(0, .10, .20, .30, .40, .50, .60, .70, .80, 1)
rclmat  <- cbind(breaks[-length(breaks)], breaks[-1], 1:9)

grazing_cat <- classify(annual_util, rclmat,
                        include.lowest = TRUE, right = FALSE)
names(grazing_cat) <- "grazing_class"

cat("layers in util_stack:", nlyr(util_stack), "\n")
print(terra::global(util_stack, c("min","max"), na.rm=TRUE))

cat_file <- sprintf("grazing_class_%s_%s.tif", it_tag, ts_tag)
writeRaster(grazing_cat, cat_file,
            datatype = "INT1U", overwrite = TRUE)

# Quick sanity checklist ----------------------------------------------------

rng <- terra::global(annual_util, "range", na.rm = TRUE)[1, ] 
cat("annual util range:", rng$min, "to", rng$max, "\n")
cat("grazing_cat frequencies:\n")
print(terra::freq(grazing_cat))

util_file <- sprintf("annual_util_%s_%s.tif", it_tag, ts_tag)
writeRaster(annual_util, util_file,
            datatype = "FLT4S", overwrite = TRUE)

if (all(is.na(annual_util[]))) stop("annual_util is all NA – check earlier steps")
## --------------------------------------------------------------------------
## 3.  make 9 binary rasters  ➜  append rows to SyncroSim
## --------------------------------------------------------------------------
#workingDir <- "C:/Local data/Modelling/Heather_grasssimplev3.ssim.temp"
  
newRows <- vector("list", 9)            # collect rows, then rbind()

for (i in seq_along(labels)) {
  
  fname <- sprintf("grazing_%d_%s_%s_%s.tif",
                   i, labels[i], it_tag, ts_tag)
  
  dest <- file.path(workingDir, fname)
  
  bin <- as.int(grazing_cat == i)
  names(bin) <- labels[i]
  
  writeRaster(bin, dest, datatype = "INT1U", overwrite = TRUE)
  
  newRows[[i]] <- data.frame(
    Iteration          = iteration,
    Timestep           = timestep,
    TransitionGroupId  = transition_group_ids[i],
    MultiplierFileName = dest,
    stringsAsFactors   = FALSE
  )
}

## --------------------------------------------------------------------------
## 4.  append (don’t overwrite) the nine new rows
## --------------------------------------------------------------------------
saveDatasheet(
  myParentScenario,
  do.call(rbind, newRows),
  name       = Grazing_levels,
  append     = TRUE,        # <— crucial: keep previous timesteps
  breakpoint = TRUE         # memory-friendly for long runs
)