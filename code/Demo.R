# Code to demonstrate running of optimization code
# This code uses the data underlying our main analysis, but filters to a subset for 
# the purposes of the demo (there is no separate demo dataset in the repository)

# Load packages
library(here)
library(data.table)

# Load the code to solve
source(here("code", "SolveBulldozersProbabilistic.R"))

############
## Set up problem parameters
##############

# Core ecoregion information.
ecoregion.data = fread(here("data", "alldata.csv"))  
ecoregion.data = ecoregion.data[!is.na(GDPperm2) & GDPperm2>0,]
# Mammal-based information on complementarity/species overlap across ecoregions
profile.data = fread(here("data", "mammal_based_profiles.csv"))
profile.abundance = profile.data[,filt.abundances]               
profiles = as.matrix(profile.data[,-ncol(profile.data), with=F])

# For the purposes of the demo, we'll just filter to a few ecoregions
# and corresponding profiles
kNumDemoRegions = 20
ecoregion.data = ecoregion.data[1:kNumDemoRegions,]
profiles = profiles[,1:kNumDemoRegions]
profiles.to.keep = which(rowsums(profiles)>0)
profiles = profiles[profiles.to.keep,]

profile.abundance = profile.abundance[profiles.to.keep]

# Optimization parameters: few periods and low budget so that optimization completes quickly
kNumPeriods = 5
kBudgetPerPeriod = 1E4
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

# Calculate area which *could* be reforested and its implications
# for the SAR multiplier
ecoregion.data[,forestable.area.init:=forest_area_2018+18*get(kDefType)]
ecoregion.data[,ae:=exp(log(plant_spcs) - kSARExponent*log(forestable.area.init))]

# Drop regions with missing or clearly bad cost data (i.e. cost should be positive), and
land.setup = ecoregion.data[!is.na(GDPperm2) & GDPperm2>0,
                            .(eco_code = eco_code, 
                              area.forestable.init = forestable.area.init,
                              area.forested.init = get(kFAType),
                              area.reserved.init = get(kPFAType),
                              species = plant_spcs,
                              ce = get(kCostType),
                              def.rate = get(kDefType),
                              ae=ae)] 
land.setup[,ce:=ce/((area.forestable.init-area.reserved.init)^(1/kCostEl))]


##########
# Main model involving complementarity, endogenous costs, and restoration
SolveBulldozersProbabilistic(land.setup,
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
                             out.file=here("results","demo.Rdata"))

# Load and inspect spend allocation.
# Budget matrix has rows of ecoregions and columns of periods.
load(here("results","demo.Rdata"))
global.soln$budget