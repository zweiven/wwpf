# R code to run all extra simulations in the Supplementary Information

# Load packages
library(here)
library(data.table)
library(sf)

# Load the code to solve
source(here("code", "SolveBulldozersProbabilistic.R"))
source(here("code", "SolveBulldozersHeuristic.R"))

############
## Set up problem parameters
##############

# Core ecoregion information
ecoregion.data = fread(here("data", "alldata.csv"))  
# Mammal-based information on complementarity/species overlap across ecoregions
profile.data = fread(here("data", "mammal_based_profiles.csv"))
profile.abundance = profile.data[,filt.abundances]               
profiles = as.matrix(profile.data[,-ncol(profile.data), with=F])

# Baseline fixed parameters
kNumPeriods = 50
kBudgetPerPeriod = 1E9
kSARExponent = 0.2
kCostType = "GDPperm2_partial"
kDefSource = "hansen2018gross"
kDefType = "dozer_rate_2000_2018"
kFAType = "forest_area_2018"
kPFAType = "protected_forest_area_2018"
kCostEl = -6
# Endemic profiles with restoration value at 80% of pristine forest
# (see Newbold et al. 2015); 
kRestMultiplier = 0.8 

# For heuristic solutions
kBudgetIncr = 100000


# Calculate area which *could* be reforested and its implications
# for the SAR multiplier
ecoregion.data[,forestable.area.init:=forest_area_2018+18*get(kDefType)]
ecoregion.data[,ae:=exp(log(plant_spcs) - kSARExponent*log(forestable.area.init))]

# Drop regions with missing or clearly bad cost data (i.e. cost should be positive), and
base.land.setup = ecoregion.data[!is.na(GDPperm2) & GDPperm2>0,
                            .(eco_code = eco_code, 
                              area.forestable.init = forestable.area.init,
                              area.forested.init = get(kFAType),
                              area.reserved.init = get(kPFAType),
                              species = plant_spcs,
                              ce = get(kCostType),
                              def.rate = get(kDefType),
                              ae=ae)] 
base.land.setup[,ce:=ce/((area.forestable.init-area.reserved.init)^(1/kCostEl))]

# Set up alternative species profiles under assumption of endemism
endemic.profiles = diag(nrow(base.land.setup))
endemic.abundances = base.land.setup$species


##########
# Alternate z values (SAR exponents)
altz.land.setup = copy(base.land.setup)
altz = 0.1
altz.land.setup[,ae:=exp(log(species) - altz*log(area.forestable.init))]

SolveBulldozersProbabilistic(altz.land.setup,
                             num.periods = kNumPeriods,
                             budget.per.period = kBudgetPerPeriod, 
                             profiles = profiles,
                             profile.abundance = profile.abundance,
                             cost.el = kCostEl,
                             restoration.value = kRestMultiplier,
                             budget.conv.tol = 100,
                             conv.tol = 5E-8,
                             inner.conv.tol = 1E-9,
                             realloc.tol = 0.1, 
                             budget.nudge = 0.1,
                             out.file=here("results","altz_0.1.Rdata"))

altz = 0.3
altz.land.setup[,ae:=exp(log(species) - altz*log(area.forestable.init))]

SolveBulldozersProbabilistic(altz.land.setup,
                             num.periods = kNumPeriods,
                             budget.per.period = kBudgetPerPeriod, 
                             profiles = profiles,
                             profile.abundance = profile.abundance,
                             cost.el = kCostEl,
                             restoration.value = kRestMultiplier,
                             budget.conv.tol = 100,
                             conv.tol = 5E-8,
                             inner.conv.tol = 1E-9,
                             realloc.tol = 0.1, 
                             budget.nudge = 0.1,
                             out.file=here("results","altz_0.3.Rdata"))

##########
# Geographically varying SAR exponents (z values)
##########

# z drawn from Storch, Keil, & Jetz (2012)

load(here("data", "ecoregions_simplified.Rdata"))
sf::sf_use_s2(FALSE)
ecoregions.simple = st_buffer(ecoregions.simple, dist=0) # clean up self-intersections
regions = read_sf(here("data", "worldregions", "region.shp"))
ecoregion.region.overlap = st_intersection(regions, ecoregions.simple)
sf::sf_use_s2(TRUE)
ecoregion.region.overlap$overlap.area = as.numeric(st_area(ecoregion.region.overlap))
ecoregion.region.overlap.dt = data.table(ecoregion.region.overlap)
max.overlap.dt = ecoregion.region.overlap.dt[, .(region = REGION[overlap.area==max(overlap.area)], 
                                                 overlap.area = overlap.area[overlap.area==max(overlap.area)]),
                                             by=.(ECO_CODE)]
# Group into regions at scale used in data source for z-values
max.overlap.dt[,lgregion:=region]

max.overlap.dt[region %in% c("Southern Africa",
                             "Eastern Africa",
                             "Middle Africa",
                             "Northern Africa",
                             "Western Africa"),
               lgregion:="Africa"]
max.overlap.dt[region %in% c("Northern Europe",
                             "Eastern Europe",
                             "Southern Europe",
                             "Western Europe",
                             "European Russia",
                             "Southern Asia",
                             "Eastern Asia",
                             "Central Asia",
                             "Southeastern Asia",
                             "Western Asia",
                             "Asiatic Russia"),
               lgregion:="Eurasia"]
max.overlap.dt[region %in% c("Northern America"), 
               lgregion:="N. America"]
max.overlap.dt[region %in% c("South America"), 
               lgregion:="S. America"]
max.overlap.dt[region %in% c("Australia/New Zealand"), 
               lgregion:="Australia"]

# The following regions aren't really well represented in the
# study with z-values. We'll try to group sensibly with regions
# where we have data OR group together & use placeholder value
max.overlap.dt[region %in% c("Central America"), 
               lgregion:="S. America"]

max.overlap.dt[region %in% c("Melanesia",
                             "Micronesia",
                             "Polynesia",
                             "Caribbean"), 
               lgregion:="Islands"]


max.overlap.dt[,z:=0.216] # For any regions not overwritten below, default is mean of values from other regions
max.overlap.dt[lgregion=="Eurasia", z:=0.26]
max.overlap.dt[lgregion=="Africa", z:=0.23]
max.overlap.dt[lgregion=="N. America", z:=0.21]
max.overlap.dt[lgregion=="S. America", z:=0.17]
max.overlap.dt[lgregion=="Australia", z:=0.21]

storch.land.setup = copy(base.land.setup)
storch.land.setup = merge(storch.land.setup, max.overlap.dt[,.(eco_code=ECO_CODE, ze=z)], by=c("eco_code"), all.x=T)
# one ecoregion missing from polygon overlap. It's a small island ecoregion in oceania (OC0204).
# So for consistency with other islands in the region, set to z=0.216 like others without data
storch.land.setup[is.na(ze), ze:=0.216]
storch.land.setup[,ae:=exp(log(species) - ze*log(area.forestable.init))]
SolveBulldozersProbabilistic(storch.land.setup,
                             num.periods = kNumPeriods,
                             budget.per.period = kBudgetPerPeriod, 
                             profiles = profiles,
                             profile.abundance = profile.abundance,
                             cost.el = kCostEl,
                             restoration.value = kRestMultiplier,
                             budget.conv.tol = 100,
                             conv.tol = 5E-8,
                             inner.conv.tol = 1E-9,
                             realloc.tol = 0.1, 
                             budget.nudge = 0.1,
                             out.file=here("results","altz_storch.Rdata"))



# Alternate ecoregion-specific z values drawn from Table 1 of Kier et al. (2005)
# 
kier.land.setup = copy(base.land.setup)
kier.land.setup = merge(kier.land.setup, max.overlap.dt[,.(eco_code=ECO_CODE, region=region, lgregion=lgregion)], by=c("eco_code"), all.x=T)
kier.land.setup = merge(kier.land.setup, ecoregion.data[,.(eco_code, wwf_mhtnam)], by="eco_code")

kier.land.setup[,ze:=as.numeric(NA)]
kier.land.setup[wwf_mhtnam=="Tropical and Subtropical Moist Broadleaf Forests"
           & region %in% c("Eastern Asia","Melanesia","Micronesia","Polynesia","Southeastern Asia","Southern Asia"),
           ze:=0.26]
kier.land.setup[wwf_mhtnam=="Tropical and Subtropical Moist Broadleaf Forests"
           & region %in% c("Central America", "Caribbean"),
           ze:=0.33]
kier.land.setup[wwf_mhtnam=="Tropical and Subtropical Moist Broadleaf Forests"
           & region %in% c("South America"),
           ze:=0.32]
kier.land.setup[wwf_mhtnam=="Tropical and Subtropical Moist Broadleaf Forests"
           & region %in% c("Australia/New Zealand","Eastern Africa","Middle Africa","Southern Africa", "Western Africa"),
           ze:=0.24]
kier.land.setup[wwf_mhtnam=="Tropical and Subtropical Dry Broadleaf Forests",
           ze:=0.21]
kier.land.setup[wwf_mhtnam=="Tropical and Subtropical Coniferous Forests",
           ze:=0.19]
kier.land.setup[wwf_mhtnam=="Temperate Broadleaf and Mixed Forests",
           ze:=0.17]
kier.land.setup[wwf_mhtnam=="Temperate Conifer Forests",
           ze:=0.14]
kier.land.setup[wwf_mhtnam=="Boreal Forests/Taiga",
           ze:=0.16]
kier.land.setup[wwf_mhtnam=="Mediterranean Forests, Woodlands and Scrub",
           ze:=0.2]

kier.land.setup[,ae:=exp(log(species) - ze*log(area.forestable.init))]
SolveBulldozersProbabilistic(kier.land.setup,
                             num.periods = kNumPeriods,
                             budget.per.period = kBudgetPerPeriod, 
                             profiles = profiles,
                             profile.abundance = profile.abundance,
                             cost.el = kCostEl,
                             restoration.value = kRestMultiplier,
                             budget.conv.tol = 100,
                             conv.tol = 5E-8,
                             inner.conv.tol = 1E-9,
                             realloc.tol = 0.1, 
                             budget.nudge = 0.1,
                             out.file=here("results","altz_kier.Rdata"))





#################
# Alternate land cost
altcost.land.setup = ecoregion.data[!is.na(GDPperm2) & GDPperm2>0,
                                 .(eco_code = eco_code, 
                                   area.forestable.init = forestable.area.init,
                                   area.forested.init = get(kFAType),
                                   area.reserved.init = get(kPFAType),
                                   species = plant_spcs,
                                   ce = GDPperm2,
                                   def.rate = get(kDefType),
                                   ae=ae)] 
altcost.land.setup[,ce:=ce/((area.forestable.init-area.reserved.init)^(1/kCostEl))]

SolveBulldozersProbabilistic(altcost.land.setup,
                             num.periods = kNumPeriods,
                             budget.per.period = kBudgetPerPeriod, 
                             profiles = profiles,
                             profile.abundance = profile.abundance,
                             cost.el = kCostEl,
                             restoration.value = kRestMultiplier,
                             budget.conv.tol = 100,
                             conv.tol = 5E-8,
                             inner.conv.tol = 1E-9,
                             realloc.tol = 0.1, 
                             budget.nudge = 0.1,
                             out.file=here("results","altcost.Rdata"))


#################
# Alternate deforestation rates

# Based on net forest loss from Hansen et al. (2013)
# note we only use that for the rate; we don't alter forested/forestable area
hansen.net.land.setup = copy(base.land.setup)
hansen.net.land.setup[,def.rate:=NULL]
hansen.net.land.setup = merge(hansen.net.land.setup, ecoregion.data[!is.na(GDPperm2) & GDPperm2>0, .(eco_code, def.rate=dozer_rate_2000_2012_net)], by="eco_code")
SolveBulldozersProbabilistic(hansen.net.land.setup,
                             num.periods = kNumPeriods,
                             budget.per.period = kBudgetPerPeriod, 
                             profiles = profiles,
                             profile.abundance = profile.abundance,
                             cost.el = kCostEl,
                             restoration.value = kRestMultiplier,
                             budget.conv.tol = 100,
                             conv.tol = 5E-8,
                             inner.conv.tol = 1E-9,
                             realloc.tol = 0.1, 
                             budget.nudge = 0.1,
                             out.file=here("results","altdef_hansennet.Rdata"))

# Based on Copernicus data
copernicus.land.setup = copy(base.land.setup)
copernicus.land.setup[,def.rate:=NULL]
copernicus.land.setup = merge(copernicus.land.setup, ecoregion.data[!is.na(GDPperm2) & GDPperm2>0, .(eco_code, def.rate=dozer_rate_2009_2015)], by="eco_code")
SolveBulldozersProbabilistic(copernicus.land.setup,
                             num.periods = kNumPeriods,
                             budget.per.period = kBudgetPerPeriod, 
                             profiles = profiles,
                             profile.abundance = profile.abundance,
                             cost.el = kCostEl,
                             restoration.value = kRestMultiplier,
                             budget.conv.tol = 100,
                             conv.tol = 5E-8,
                             inner.conv.tol = 1E-9,
                             realloc.tol = 0.1, 
                             budget.nudge = 0.1,
                             out.file=here("results","altdef_copernicus.Rdata"))



# Using SSP projections for alternate assumptions about
# - rates of deforestation
# - costs
# - budget

# For forest cover, we'll use a fairly pessimistic projection at a global scale.
# (SSP3-Baseline)
# We do this
# - so that there are actually threats everywhere -- if there is large scale
#   reforestation embedded in the SSPs, then our problem formulation isn't as
#   relevant.
# - because the regional cover projections are not fine-grained enough to give
#   useful projections per ecoregion. Regional projections also bounce around 
#   in ways that are non-monotonic and occasionally quite sporadic, raising 
#   questions about both credibility and requiring further assumptions about
#   what happens to species in those regions.
# - so that we focus solely on positive deforestation rates. The point of our
#   approach is to globally optimize reforestation. Growing forest cover in the
#   SSPs themselves implies reforestation is happening, which could deviate from
#   optimization. We'd like (to the extent possible) to use the SSPs purely for
#   projections of cost and threat.
ssp.forestcover = fread(here("data", "iamc_db_forest.csv"))

# to use these, we'll calculate annual rates per decade, 
# normalize to the 2005-2010 rate (which falls within the time span of the historical rates),
# and then scale historical rates by those normalized factors. Since our 
# simulation starts at 2018 forest levels, we'll consider years 1-2 to follow 
# 2010-2020 rates, then shift to 2020-2030 rates, etc.

hist.rate = ssp.forestcover[Scenario=="SSP3-Baseline" & Region=="World",(`2005`-`2010`)/5]
rel.future.rates = diff(as.numeric(ssp.forestcover[Scenario=="SSP3-Baseline" & Region=="World",7:16, with=F]))/-10/hist.rate

ssp.land.setup = copy(base.land.setup)

ssp.land.setup$def.rate = lapply(ssp.land.setup$def.rate, FUN=function(hist.def.rate) {
  c(rep(hist.def.rate*rel.future.rates[1], times=2),  #2018-2019
    rep(hist.def.rate*rel.future.rates[2], times=10), #2020-2029
    rep(hist.def.rate*rel.future.rates[3], times=10), #2030-2039
    rep(hist.def.rate*rel.future.rates[4], times=10), #2040-2049
    rep(hist.def.rate*rel.future.rates[5], times=10), #2050-2059
    rep(hist.def.rate*rel.future.rates[5], times=8)   #2060-2067
  )
})

# For costs, we'll use the GDP projections and GDP per capita projections to scale
# land cost multipliers for conservation. We will use regional projections here, together
# with a mapping from ecoregions to the 5 SSP regions.
# For consistency with the deforestation data, we'll again use
# SSP3-Baseline.
ssp.gdp = fread(here("data","iamc_db_gdp.csv"))
ssp.pop = fread(here("data", "iamc_db_pop.csv"))
region.crosswalk = fread(here("data","ecoregions_to_sspregions.csv"))
# here, since we're interested in levels rather than rates of change, we'll
# calculate gdppc and linearly interpolate (gdp per capita) levels provided in
# different years under the SSPs
target.yrs = 2018:2067
filt.gdp.cap = as.matrix(ssp.gdp[Scenario=="SSP3-Baseline" & Region!="World",7:16]/ssp.pop[Scenario=="SSP3-Baseline" & Region!="World",7:16])

decades = seq(from=2010, to=2100, by=10)
gdppc.proj = matrix(NA, nrow=nrow(filt.gdp.cap), ncol=50)
for(yr in target.yrs) {
  dec = as.integer(cut(yr, decades, right=F))
  gdppc.proj[,yr-2017]=filt.gdp.cap[,dec]+(yr-decades[dec])/10*(filt.gdp.cap[,dec+1]-filt.gdp.cap[,dec])
}
gdppc.proj = gdppc.proj/gdppc.proj[,1]
row.names(gdppc.proj) = substr(ssp.gdp[Scenario=="SSP3-Baseline" & Region!="World",Region], start=5, stop=10)

# Now, we need to align ecoregions and regions, and use these multipliers
# to scale and create costs. Note the land cost multipliers (ce) have already
# been computed in base.land.setup (which was the starting point for ssp.land.setup),
# so the code here just rescales those based on ratios of future gdppc to current 
# (ratios computed above). As a result, the ce values should be the right ones to use.
ssp.land.setup$ce = mapply(ssp.land.setup$eco_code, 
                           ssp.land.setup$ce, 
                           SIMPLIFY=F,
                           FUN=function(eco_code, hist.cost) {
                             hist.cost * gdppc.proj[rownames(gdppc.proj)==region.crosswalk[ECO_CODE==eco_code, ssp.region],]
                           })

# Third, we'll consider the budget to scale along with rising land costs.
# We'll take projections of global GDP per capita, and assume that the annual
# budget scales up and down with that.
global.gdppc.proj = as.numeric(rep(NA, times=50))
global.gdp.cap = as.numeric(ssp.gdp[Scenario=="SSP3-Baseline" & Region=="World",7:16]/ssp.pop[Scenario=="SSP3-Baseline" & Region=="World",7:16])
for(yr in target.yrs) {
  dec = as.integer(cut(yr, decades, right=F))
  global.gdppc.proj[yr-2017]=global.gdp.cap[dec]+(yr-decades[dec])/10*(global.gdp.cap[dec+1]-global.gdp.cap[dec])
}
global.gdppc.proj = global.gdppc.proj/global.gdppc.proj[1]

# Run with SSP-derived projections of deforestation and costs

SolveBulldozersProbabilistic(ssp.land.setup,
                             num.periods = kNumPeriods,
                             budget.per.period = kBudgetPerPeriod*global.gdppc.proj, 
                             profiles = profiles,
                             profile.abundance = profile.abundance,
                             cost.el = kCostEl,
                             restoration.value = kRestMultiplier,
                             budget.conv.tol = 100,
                             conv.tol = 5E-8,
                             inner.conv.tol = 1E-9,
                             realloc.tol = 0.1, 
                             budget.nudge = 0.1,
                             out.file=here("results","altdef_ssp3.Rdata"))

# Run again, but at half SSP3 deforestation rates
ssp.halfdef.land.setup = copy(ssp.land.setup)
ssp.halfdef.land.setup[,def.rate:=lapply(ssp.halfdef.land.setup$def.rate, FUN=function(ecoregion.def.rate) {
  ecoregion.def.rate/2
})]
SolveBulldozersProbabilistic(ssp.halfdef.land.setup,
                             num.periods = kNumPeriods,
                             budget.per.period = kBudgetPerPeriod*global.gdppc.proj, 
                             profiles = profiles,
                             profile.abundance = profile.abundance,
                             cost.el = kCostEl,
                             restoration.value = kRestMultiplier,
                             budget.conv.tol = 100,
                             conv.tol = 5E-8,
                             inner.conv.tol = 1E-9,
                             realloc.tol = 0.1, 
                             budget.nudge = 0.1,
                             out.file=here("results","altdef_ssp3_halved.Rdata"))

# Based on halving current rates based on Hansen et al. (2013) gross forest loss
lower.def.setup = copy(land.setup)
lower.def.setup[,def.rate:=def.rate/2]

SolveBulldozersProbabilistic(lower.def.setup,
                             num.periods = kNumPeriods,
                             budget.per.period = kBudgetPerPeriod, 
                             profiles = profiles,
                             profile.abundance = profile.abundance,
                             cost.el = kCostEl,
                             restoration.value = kRestMultiplier,
                             budget.conv.tol = 100,
                             conv.tol = 5E-8,
                             inner.conv.tol = 1E-9,
                             realloc.tol = 0.1, 
                             budget.nudge = 0.1,
                             out.file=here("results","half_def_rate.Rdata"))



#################
# Alternate planning horizons
SolveBulldozersProbabilistic(base.land.setup,
                             num.periods = 100,
                             budget.per.period = kBudgetPerPeriod, 
                             profiles = profiles,
                             profile.abundance = profile.abundance,
                             cost.el = kCostEl,
                             restoration.value = kRestMultiplier,
                             budget.conv.tol = 100,
                             conv.tol = 5E-8,
                             inner.conv.tol = 1E-9,
                             realloc.tol = 0.1, 
                             budget.nudge = 0.1,
                             out.file=here("results","althorizon_100.Rdata"))

SolveBulldozersProbabilistic(base.land.setup,
                             num.periods = 200,
                             budget.per.period = kBudgetPerPeriod, 
                             profiles = profiles,
                             profile.abundance = profile.abundance,
                             cost.el = kCostEl,
                             restoration.value = kRestMultiplier,
                             budget.conv.tol = 100,
                             conv.tol = 5E-8,
                             inner.conv.tol = 1E-9,
                             realloc.tol = 0.1, 
                             budget.nudge = 0.1,
                             out.file=here("results","althorizon_200.Rdata"))


##################################################
# Models missing one of three key model features: reforestation, isoelastic costs, and complementarity
##################################################

# Restoration and endogenous costs (no species complementarity/assume endemism)
SolveBulldozersProbabilistic(base.land.setup,
                             num.periods = kNumPeriods,
                             budget.per.period = kBudgetPerPeriod, 
                             profiles = endemic.profiles,
                             profile.abundance = endemic.abundances,
                             cost.el = kCostEl,
                             restoration.value = kRestMultiplier,
                             budget.conv.tol = 100,
                             conv.tol = 1E-8,
                             inner.conv.tol = 1E-8,
                             realloc.tol = 1, 
                             budget.nudge = 0.1,
                             out.file=here("results","nocomplementarity.Rdata"))

# Endogenous costs and complementarity (but no restoration)
SolveBulldozersProbabilistic(base.land.setup,
                             num.periods = kNumPeriods,
                             budget.per.period = kBudgetPerPeriod, 
                             profiles = profiles,
                             profile.abundance = profile.abundance,
                             cost.el=kCostEl,
                             restoration.value = 0,
                             budget.conv.tol = 100,
                             conv.tol = 1E-8,
                             inner.conv.tol = 1E-8,
                             realloc.tol = 1, 
                             budget.nudge = 0.1,
                             out.file=here("results","norestoration.Rdata"))

# Restoration and complementarity (but exogenous costs)
fixedcost.land.setup = copy(base.land.setup)
fixedcost.land.setup[,ce:=GDPperm2_partial]

SolveBulldozersProbabilistic(fixedcost.land.setup,
                             num.periods = kNumPeriods,
                             budget.per.period = kBudgetPerPeriod, 
                             profiles = profiles,
                             profile.abundance = profile.abundance,
                             cost.el=NA,
                             restoration.value = kRestMultiplier,
                             budget.conv.tol = 100,
                             conv.tol = 1E-8,
                             inner.conv.tol = 1E-8,
                             realloc.tol = 1, 
                             budget.nudge = 0.1,
                             out.file=here("results","norisingcosts.Rdata"))


# Model missing all three key model features: reforestation, isoelastic costs, and complementarity
SolveBulldozersProbabilistic(fixedcost.land.setup,
                             num.periods = kNumPeriods,
                             budget.per.period = kBudgetPerPeriod, 
                             profiles = endemic.profiles,
                             profile.abundance = endemic.profile.abundance,
                             cost.el = NA,
                             restoration.value = 0,
                             budget.conv.tol = 100,
                             conv.tol = 1E-8,
                             inner.conv.tol = 1E-8,
                             realloc.tol = 1, 
                             budget.nudge = 0.1,
                             out.file=here("results","noextensions.Rdata"))
