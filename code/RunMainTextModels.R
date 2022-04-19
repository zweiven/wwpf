# Code for running optimization models underlying
# results in main text.

# Load packages
library(here)
library(data.table)

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

# Budget increment for heuristic solutions
kBudgetIncr = 100000

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

# Set up alternative species profiles under assumption of endemism
endemic.profiles = diag(nrow(land.setup))
endemic.abundances = land.setup$species


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
                             out.file="results/final/main.Rdata")


# Alternate budget models
for(b in c(1E7,5E7,1E8,2.5E8,5E8,7.5E8,1E9,1E10,2.5E10)) {
#for(b in c(2.5E8,5E8,7.5E8,1E9,1E10,2.5E10)) {
  SolveBulldozersProbabilistic(land.setup,
                               num.periods = kNumPeriods,
                               budget.per.period = b, 
                               profiles = profiles,
                               profile.abundance = profile.abundance,
                               cost.el = kCostEl,
                               restoration.value = kRestMultiplier,
                               budget.conv.tol = 100,
                               conv.tol = 1E-8,
                               inner.conv.tol = 1E-8,
                               realloc.tol = 1, 
                               budget.nudge = 0.1,
                               out.file=paste0("results/final/altbudget_",b,".Rdata"))
}

# What if we had all budget up front? We don't discount for consistency with 
# how we treat costs.
single.per.budget = sum(rep(kBudgetPerPeriod, kNumPeriods))
SolveBulldozersProbabilistic(land.setup,
                             num.periods = 1,
                             budget.per.period = single.per.budget, 
                             profiles = profiles,
                             profile.abundance = profile.abundance,
                             cost.el = kCostEl,
                             restoration.value = kRestMultiplier,
                             budget.conv.tol = 100,
                             conv.tol = 1E-8,
                             inner.conv.tol = 1E-8,
                             realloc.tol = 1, 
                             budget.nudge = 0.1,
                             out.file=paste0("results/final/upfrontbudget.Rdata"))


# Model variants dropping one of three features:
# - Species complementarity
# - Endogenous land costs
# - Restoration

# Restoration and endogenous costs (no species complementarity/assume endemism)
SolveBulldozersProbabilistic(land.setup,
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
                             out.file="results/final/endemism.Rdata")

# Restoration and complementarity (but exogenous costs)

SolveBulldozersProbabilistic(land.setup,
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
                             out.file="results/final/fixedmc.Rdata")

# Endogenous costs and complementarity (but no restoration)

SolveBulldozersProbabilistic(land.setup,
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
                             out.file="results/final/norestoration.Rdata")


############################
# Heuristic solutions
############################

SolveBulldozersProblemHeuristic(land.setup,
                                 num.periods = kNumPeriods,
                                 budget.per.period = kBudgetPerPeriod, 
                                 profiles = profiles,
                                 profile.abundance = profile.abundance,
                                 cost.el = kCostEl,
                                 restoration.value = kRestMultiplier,
                                 heur.type="bcratio",
                                 budget.incr = kBudgetIncr,
                                 out.file="results/final/heur_bcratio.Rdata")


SolveBulldozersProblemHeuristic(land.setup,
                                 num.periods = kNumPeriods,
                                 budget.per.period = kBudgetPerPeriod, 
                                 profiles = profiles,
                                 profile.abundance = profile.abundance,
                                 cost.el = kCostEl,
                                 restoration.value = kRestMultiplier,
                                 heur.type="mb",
                                 budget.incr = kBudgetIncr,
                                 out.file="results/final/heur_mb.Rdata")


SolveBulldozersProblemHeuristic(land.setup,
                                 num.periods = kNumPeriods,
                                 budget.per.period = kBudgetPerPeriod, 
                                 profiles = profiles,
                                 profile.abundance = profile.abundance,
                                 cost.el = kCostEl,
                                 restoration.value = kRestMultiplier,
                                 heur.type="hotspot",
                                 budget.incr = kBudgetIncr,
                                 out.file="results/final/heur_hotspot.Rdata")


SolveBulldozersProblemHeuristic(land.setup,
                                 num.periods = kNumPeriods,
                                 budget.per.period = kBudgetPerPeriod, 
                                 profiles = profiles,
                                 profile.abundance = profile.abundance,
                                 cost.el = kCostEl,
                                 restoration.value = kRestMultiplier,
                                 heur.type="threat",
                                 budget.incr = kBudgetIncr,
                                 out.file="results/final/heur_threat.Rdata")

SolveBulldozersProblemHeuristic(land.setup,
                                num.periods = kNumPeriods,
                                budget.per.period = kBudgetPerPeriod, 
                                profiles = profiles,
                                profile.abundance = profile.abundance,
                                cost.el = kCostEl,
                                restoration.value = kRestMultiplier,
                                heur.type="cost",
                                budget.incr = kBudgetIncr,
                                out.file="results/final/heur_cost.Rdata")

SolveBulldozersProblemHeuristic(land.setup,
                                num.periods = kNumPeriods,
                                budget.per.period = kBudgetPerPeriod, 
                                profiles = profiles,
                                profile.abundance = profile.abundance,
                                cost.el = kCostEl,
                                restoration.value = kRestMultiplier,
                                heur.type="percprot",
                                budget.incr = kBudgetIncr,
                                out.file="results/final/heur_percprot.Rdata")
