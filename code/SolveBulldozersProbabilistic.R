# Set of functions that solve the general
# problem of dynamically optimizing spend
# across a set of regions to maximize expected
# number of species protection.

##############
# Dependencies
##############
library(Rfast) # for colprods

##############
# Constants
##############
# Tolerance for detecting when land has run out
EPSILON = 1E-3   
# Number of budget reallocation iterations to temporarily ignore an ecoregion
# if it is simply acting as a pass-through for budget transfers from
# one region to another. This occurs if a region is on a steep part
# of the ROI curve and can oscillate between highest and lowest ROI.
# This is just to speed things up; we will always return to the ecoregion.
NIGNOREITER = 10  

# First, define a helper function that characterizes the implications of a particular budget allocation
# Treat all problem parameters as (global) constants -- make this a function of the current candidate budget choice only.
# Gives back a list with the following elements given the candidate budget allocation:
# - cons.end.periods: the end period (race over period) per region
# - area.avail.consendper: the area still available during the race over period per region
# - end.areas.conserved: areas conserved at race end per region 
# - rois: the marginal ROI at the end of the race per region 
# Note the addition of profiles and profile.abundance. 
# - profiles is a matrix in which rows represent types of species, columns represent ecoregions, 
#   and entries indicate presence absence of that species type in that ecoregion.
# - profile.abundances is a vector with length equal to the number of rows in profile.abundance, indicating
#   how many species of that type we think there are.


FindRaceEndConditionsProbabilistic = function(cur.budget, 
                                              land.setup, 
                                              cum.devel.mat, 
                                              profiles, 
                                              profile.abundance, 
                                              cost.el=NA, 
                                              restoration.value=0, 
                                              roi.period) {
  num.regions = nrow(land.setup)
  num.periods = nrow(cum.devel.mat)
  anti.profiles = 1 - profiles
  fixed.mc = is.na(cost.el)
  
  
  
  with(land.setup,
       {
    # Find end period. This is the first period in which, given rates of deforestation and protection,
    # if both sides did not collide, the available land at the end of that period would be non-positive.
    
    # Cumulative land bought per current allocation at end of each period
    # - calculate cumulative spend as of period t, then compute land purchased using that vector
    cum.spend = apply(cur.budget, 1, cumsum)
    initial.forested.mat = matrix(area.forested.init, ncol=num.regions, nrow=num.periods, byrow=T)
    initial.forestable.mat = matrix(area.forestable.init, ncol=num.regions, nrow=num.periods, byrow=T)
    initial.protected.mat = matrix(area.reserved.init, ncol=num.regions, nrow=num.periods, byrow=T)
    if(class(ce)=="list") {
      ce.mat = matrix(unlist(ce), ncol=num.regions, nrow=num.periods, byrow=F)
    } else {
      ce.mat = matrix(ce, ncol=num.regions, nrow=num.periods, byrow=T)
    }
    
    # calculate hypothetical cumulative *new* area reserved by end of each period
    # if not constrained by bulldozers
    if(fixed.mc) {
      res.totals =  apply(cur.budget/t(ce.mat), 1, cumsum) 
    } else {
      if(cost.el!=-1) {
        ep.oneplusep = cost.el/(1+cost.el)
        if(class(ce)=="list") {
          # (potentially) time-varying cost multiplier per ecoregion
          res.totals = initial.forestable.mat - ((initial.forestable.mat-initial.protected.mat)^(1/ep.oneplusep) - apply(cur.budget/(t(ce.mat)*ep.oneplusep), 1, cumsum))^ep.oneplusep - initial.protected.mat
        } else {
          # time-invariant cost multiplier per ecoregion
          res.totals = initial.forestable.mat - ((initial.forestable.mat-initial.protected.mat)^(1/ep.oneplusep) - cum.spend/(ce.mat*ep.oneplusep))^ep.oneplusep - initial.protected.mat
        }
      } else {
        stop("cost elasticity of -1 not implemented")
      }
      # If we end up with NAs/NaNs, it's because there's more than enough spend to conserve everything, and the 
      # marginal and total cost functions aren't well defined once all forestable land has been protected.
      # we'll define both to be infinite beyond that point. Equivalently, we set res.totals to all available forestable land in those cases.
      res.totals[is.na(res.totals)] = (initial.forestable.mat-initial.protected.mat)[is.na(res.totals)]
      
      # We could end up with tiny negative numbers due to apparent numerical precision issues. If cum.spend is zero 
      # for a period and region, raising (initial.forestable.mat-initial.protected.mat) to 1/ep.oneplusep then to ep.oneplusep
      # doesn't result in an entry that exactly equals initial.forestable.mat-initial.protected.mat (but it should).
      # zero those out.
      res.totals[res.totals<0] = 0
    }

    # deal with R automatically retyping as numeric for single period case
    if(num.periods==1) {
      res.totals = matrix(res.totals, ncol=num.regions, nrow=num.periods)
    }
    rem.forested.land = forested.area.available.init - t(res.totals + cum.devel.mat) # How much available land is remaining at end of each period given current plan
    
    # Find out when the 'race' is over in each region. Earlier of:
    # a) when conservation meets development: pristine forest runs out by the end of that period
    # b) the end of the planning horizon.
    # we use epsilon as a threshold for detecting when land has run out since numerical precision
    # issues may result in failure to detect that event.
    meeting.periods = apply(rem.forested.land, 1, FUN=function(rem.vec) { min(which(rem.vec<=EPSILON)) })
    cons.end.periods = pmin(meeting.periods, num.periods)
    
    # Next, we want land conserved in each region when the 'race' is over in each region.
    # Because the budget may be more than enough to conserve all remaining land in the 
    # ending period for a particular region, we'll calculate land conserved at the start of
    # that period, calculate how much is conserved in the 'race over' period subject to 
    # feasibility constraints, and add both to the land initially reserved in a region.
    
    # Compute land newly reserved at start of conservation ending period 
    area.res.consendper.start = sapply(1:num.regions, function(r) {
      # find available land at start of period in which race ends.
      race.over.per = cons.end.periods[r]
      if(race.over.per > 1) {
        area.newly.res = res.totals[race.over.per-1, r]
      } else {
        area.newly.res = 0
      }
      return(area.reserved.init[r]+area.newly.res)
    })
        
    # Now find actual land conserved during final conservation period in each region.
    # To do so, we'll need to know constraints: how much is available in the final 
    # period during which conservation takes place
    area.avail.consendper = sapply(1:num.regions, function(r) {
      # find available land at start of period in which race ends.
      cons.end.per = cons.end.periods[r]
      if(cons.end.per > 1) {
        area.avail.consendper = rem.forested.land[r, cons.end.per-1]
      } else {
        area.avail.consendper = forested.area.available.init[r]
      }
      return(area.avail.consendper)
    })
    area.newly.conserved = 
      sapply(1:num.regions, function(r) {  # land added during race
        
        # find available land at start of period in which race ends.
        race.over.per = cons.end.periods[r]
        

        # Need to properly handle case in which race ends in first period
        if(race.over.per==1) {
          return(min(area.avail.consendper[r], res.totals[race.over.per,r]))
        } else {
          # In final period, the amount conserved is the minimum of that implied by the budget 
          # and the amount still available. Where the amount implied by the budget is larger, we
          # consider that to be wasted budget rather than an infeasible allocation.
          area.conserved.lastper = min(area.avail.consendper[r], res.totals[race.over.per,r]-res.totals[race.over.per-1,r])
          return(res.totals[race.over.per-1,r] + area.conserved.lastper)
        }
      })
    res.totals.at.consend = area.reserved.init + area.newly.conserved # add in land initially conserved
    # If restoration is allowed, we need to know conditions for when the restoration phase ends.
    # It will be the earlier of 
    # - when all forested land is either conserved or restored
    # - the end of the planning horizon
    # Restoration begins when conservation stops.
    # As specified, note restoration costs the same as protection, but is just less effective per 
    # the restoration.value multiplier. This means we can use protected land totals computed earlier, but just
    # need to know how much is conserved and how much is restored. 
    if(restoration.value > 0) {
      area.restored = res.totals-matrix(area.newly.conserved, nrow=num.periods, ncol=num.regions, byrow=T)
      area.restored[area.restored<0] = 0
      # for restoration, there is no development to consider, so we only need to constrain
      # restoration totals by the total amount of forested land initially available
      # (the sum of conservation and restoration can't exceed initially forested land)
      area.avail.for.rest = matrix(area.forestable.init - res.totals.at.consend,
                                   ncol=num.regions,
                                   nrow=num.periods, 
                                   byrow=T)
      
      # can't restore more than is available for restoration
      area.restored[area.restored>area.avail.for.rest] = area.avail.for.rest[area.restored>area.avail.for.rest]
      
      # Find out when restoration 'finishes' (if it does): when
      # area available for restoration runs out. 
      # We allow for epsilon error in detecting end periods as otherwise numerical precision 
      # issues may lead to a detection failure.
      rest.end.periods = sapply(1:num.regions, FUN=function(r) {
        rest.end = min(which((res.totals.at.consend[r]+area.restored[,r])>=(area.forestable.init[r]-EPSILON)))
      })
      # If we don't actually run out of land to restore in a region,
      # there will be no periods in which the area restored is at least as
      # large as the area available for restoration. The result of the
      # sapply call above will be Inf for that region.
      # In that case, the 'end' period for restoration is the end of the planning horizon,
      # and pmin will do the trick.
      # 
      rest.end.periods = pmin(rest.end.periods, num.periods)
      
      # Find out ending area restored per region
      rest.totals.at.restend = sapply(1:num.regions, FUN=function(r) {
        rest.end = area.restored[rest.end.periods[r],r]
      })
      # now we have a region x period tally of restoration totals (area.restored)
      # and a vector of periods during which restoration finishes (rest.end.periods)
      # and a vector of final restoration areas
    } else {
      rest.end.periods = cons.end.periods
      rest.totals.at.restend = rep(0, num.regions)
    }
    
    # Calculate ROI per region in the specified period, i.e., the marginal species benefit for a dollar spent
    # in each ecoregion in the specified period.
    # To do so in our probabilistic framework, we will compute
    # the expected increase in aggregate species protected across all ecoregions if we marginally
    # increase protection in each ecoregion. 
    # That increase is the sum over types of species in a given ecoregion of the increase in probability
    # that marginal land protection in the ecoregion will newly protect the species AND it is not protected
    # elsewhere.
    # 
    # We need to deal with two discontinuities in the ROI: 
    # 1. All forested land is conserved or developed
    # -- Here the restoration.value multiplier kicks in, which drops the marginal benefit
    # 2. All forestable land is conserved or reforested.
    # -- Here the marginal benefit drops to zero because the probability curve flattens at 1.
    # 
    # Converging on an optimal budget will work better if we treat these discontinuities specially,
    # so that we can get as close to the equimarginal principle defining the simple interior solution.
    # Otherwise we can get to oscillating budget movements where we shift from a region with, e.g., zero 
    # ROI, which then has the highest ROI, we shift budget back so it has zero, etc.
    # 
    # One option is to look at left and right ROIs (approaching from below and above), and then when
    # identifying source regions, we examine the left ROI (what would happen if we removed budget)
    # while we use the right ROI for targets (what if we added $).
 
    
    
    # Calculate probability of protection for each species type in each ecoregion
    # We assume that probability is constant across species in an ecoregion.
    # If a candidate budget allocation would result in more than the initially
    # foreseted area being protected, we cap the probability at one.
    if(restoration.value > 0) {
      p.protect = ((res.totals.at.consend + restoration.value*rest.totals.at.restend)/area.forestable.init)^ze
    } else {
      p.protect = (res.totals.at.consend/area.forestable.init)^ze
    }
    p.protect[p.protect>1]=1
    
    
    if(restoration.value > 0) {
      # marginal benefits depend upon whether the marginal unit of area affected by the marginal 
      # unit of budget in the specified period is conserved or restored
      marg.p.protect.right = marg.p.protect.left = ze*((res.totals.at.consend + restoration.value * rest.totals.at.restend)^(ze-1))/(area.forestable.init^ze)
      # Need to adjust marginal probability of protection if the marginal dollar would 
      # go toward restoration. That happens if, by the end of the roi period of interest,
      # remaining forested land that is not protected and not developed is at or below zero.
      roi.period.rest = (rem.forested.land[,roi.period]<0)
      marg.p.protect.left[roi.period.rest] = marg.p.protect.left[roi.period.rest]*restoration.value      
      roi.period.rest = (rem.forested.land[,roi.period]<=0)
      marg.p.protect.right[roi.period.rest] = marg.p.protect.right[roi.period.rest]*restoration.value
      
    } else {
      marg.p.protect.right=marg.p.protect.left = ze*(res.totals.at.consend^(ze-1))/(area.forestable.init^ze)
    }

    # And build up terms indicating probability of each species type being unprotected in each ecoregion.
    # this will be a matrix with a column per species type, a row per ecoregion,
    # and the entry representing the probability the species of that type is unprotected
    # in that ecoregion.
    marg.p.terms.left = profiles * matrix(marg.p.protect.left, ncol=nrow(land.setup), nrow=length(profile.abundance), byrow=T)
    marg.p.terms.right = profiles * matrix(marg.p.protect.right, ncol=nrow(land.setup), nrow=length(profile.abundance), byrow=T)
    # We could get an NaN if a profile entry is zero for a region and the marginal probability of protection is Infinite
    # because there's zero land protected (neither initially nor via the current budget). In those cases 0*Inf is NaN, but
    # we should treat the marginal probability terms as zeros since those species' probability of protection can't be 
    # influenced by protection in that region.
    marg.p.terms.right[is.nan(marg.p.terms.right)] = 0
    marg.p.terms.left[is.nan(marg.p.terms.left)] = 0
    
    unprot.terms = apply(profiles, 1, function(sp.type.pa.vec) { (1-p.protect*sp.type.pa.vec)})

    # now for each ecoregion x species type, we want the probability that
    # additional protected area would newly protect that species type. 
    # We then sum that across species in an ecoregion to get MB, and divide by MC to get ROI.
    
    # we can take full products of unprotected terms across all ecoregions once since that's expensive
    # then in the region-specific loop, we divide by the unprotected term from that focal region to remove its effect.
    fullprods = colprods(unprot.terms)
    mbs.right = sapply(1:num.regions, function(i) {
      sum(profile.abundance*marg.p.terms.right[,i]*(fullprods/unprot.terms[i,]))
    })
    mbs.left = sapply(1:num.regions, function(i) {
      sum(profile.abundance*marg.p.terms.left[,i]*(fullprods/unprot.terms[i,]))
    })
    
    # We'll now compute marginal costs given allocations up to and including roi.period
    if(fixed.mc) {
      period.mcs.right = period.mcs.left = period.mcs = ce.mat[roi.period,]

    } else {
      # Need to calculate total reserved land (accounting for constraints of deforestation and forestable land)
      # in each area for comparison with forestable land, which will give us our MC
      roi.period.newly.res.totals = res.totals[roi.period,] 
      if(restoration.value==0) {
        roi.period.newly.res.totals[cons.end.periods==roi.period] =  area.newly.conserved[cons.end.periods==roi.period]
      } else {
        roi.period.newly.res.totals[rest.end.periods==roi.period] = (area.newly.conserved+rest.totals.at.restend)[rest.end.periods==roi.period]
      }
      # MC is a function of forestable land that is not yet protected via either conservation or reforestation
      if(class(ce)=="list") {
        # (potentially) time-varying cost multiplier per ecoregion
        period.mcs = ce.mat[roi.period,]*(area.forestable.init-(roi.period.newly.res.totals+area.reserved.init))^(1/cost.el)
      } else {
        # fixed cost multiplier per ecoregion
        period.mcs = ce*(area.forestable.init-(roi.period.newly.res.totals+area.reserved.init))^(1/cost.el)
      }
      # note again this is not well defined if all available land has been conserved via protection or reforestation, in which case 
      # MC should be infinite and no additional land can be conserved (we correct res.totals earlier).
      period.mcs.right = period.mcs.left = period.mcs
      period.mcs.right[res.totals[roi.period,]>=(area.forestable.init-area.reserved.init)] = Inf
      period.mcs.left[res.totals[roi.period,]>=(area.forestable.init-area.reserved.init)] = Inf
    }
    rois.right = mbs.right/period.mcs.right
    rois.left = mbs.left/period.mcs.left
    
    # In the edge case in which both ze and res.totals.at.race.over are zero, this will be NaN.
    # If no land is conserved AND conserving is useless (such that ze=0), then the ROI must 
    # be zero since spending can't increase forested area conserved
    rois.right[is.nan(rois.right)] = 0
    rois.left[is.nan(rois.left)] = 0
    
    
    
    # give back info about ending conditions, ROI, etc.
    if(restoration.value > 0) {
      race.end.periods = rest.end.periods
    } else {
      race.end.periods = cons.end.periods
    }
    
    
    # Calculate two things about the final period during which 
    # either conservation or restoration takes place. 
    # 1. area reserved or restored at start of that period
    area.resrest.raceendper.start = sapply(1:num.regions, function(r) {
      # find available land at start of period in which race ends.
      race.over.per = race.end.periods[r]
      if(race.over.per > 1) {
        area.newly.resrest = res.totals[race.over.per-1, r]
      } else {
        area.newly.resrest = 0
      }
      return(area.reserved.init[r]+area.newly.resrest)
    })
    
    # 2. Area available at the start of that period (end of preceding period).
    
    # note what we consider "available" in the race end period
    # will differ based on whether or not restoration is allowed.
    # - If yes, it's remaining initially forestable land that is not yet conserved.
    # - If no, it's remaining initially forested land that is neither conserved nor developed.
    area.avail.raceendper = sapply(1:num.regions, function(r) {
      # find available land at start of period in which race ends.
      race.end.per = race.end.periods[r]
      if(race.end.per > 1) {
        if(restoration.value > 0) {
          area.avail.raceendper = (forestable.area.available.init[r] - res.totals[race.end.per-1,r]) 
        } else {
          area.avail.raceendper = rem.forested.land[r, race.end.per-1]
        }
      } else {
        if(restoration.value > 0) {
          area.avail.raceendper = forestable.area.available.init[r]
        } else {
          area.avail.raceendper = forested.area.available.init[r]
        }
      }
      return(area.avail.raceendper)
    })
    
    
    return(list(cons.end.periods = cons.end.periods,
                rest.end.periods = rest.end.periods,
                race.end.periods = race.end.periods,

                area.res.consendper.start = area.res.consendper.start,
                area.resrest.raceendper.start = area.resrest.raceendper.start,
                
                area.avail.consendper = area.avail.consendper,
                area.avail.raceendper = area.avail.raceendper,
                
                end.areas.conserved = res.totals.at.consend,
                end.areas.restored = rest.totals.at.restend,
                rois.right = rois.right,
                rois.left = rois.left)) 
    
  })
}


# Function to actually solve the race as a function
# of the specified setup. Most arguments are vectors, with one entry
# per region
# Arguments:
# land.setup: data.table with information on ecoregions (one row per ecoregion)
# num.periods: number of periods in optimization (50 in base case in paper)
# budget.per.period: max budget that can be spent per period
# profiles: matrix (rows: profiles, columns: ecoregions) of presence/absence data for species profiles
# profile.abundance: vector of abundances of different species profiles
# restoration.value: multiplier indicating relative value of restored forest vs. intact forest for species protection
# cost.el: cost elasticity of land (as a proxy for more general protection costs)
# budget.conv.tol: tolerance for when budget shifts are sufficiently small to consider the algorithm to have converged
# conv.tol: tolerance for when ROI differences are sufficiently small to consider the algorithm to have converged
# inner.conv.tol: tolerance for when ROI differences are sufficiently small for inner search for budget reallocation to be complete (avoid overshooting on individual reallocations)
# realloc.tol: tolerance for when budget shifts are sufficiently small for inner search for budget reallocation to be complete (avoid overshooting on individual reallocations)
# budget.nudge: small amount by which we shift budget either direction when computing ROI for a budget increase or decrease to deal with potential discontinuities in ROI and numerical imprecision
# max.fwd.pass: limit on forward passes through the time horizon before exiting out
# out.file: name of filename to save results to. If null, results simply returned to calling code.
SolveBulldozersProbabilistic =  function(land.setup,
                                         num.periods,
                                         budget.per.period,
                                         profiles, 
                                         profile.abundance,
                                         cost.el=-999999999,
                                         restoration.value = 0,
                                         budget.conv.tol = 100,
                                         conv.tol = 1E-8,
                                         inner.conv.tol = 1E-8,
                                         realloc.tol=1E-12,
                                         budget.nudge = 0.1,
                                         max.fwd.pass=100,
                                         out.file=NULL) {


  fixed.mc = is.na(cost.el)
  land.setup[,forestable.area.available.init:= area.forestable.init-area.reserved.init] # Initially Available Area  
  land.setup[,forested.area.available.init:= area.forested.init-area.reserved.init] # Initially Available Area
  land.setup[,ze:=log(species/ae)/log(area.forestable.init)] # species area curve exponent
  
  # create some locals from the setup
  lapply(colnames(land.setup), FUN=function(cname) { 
    assign(x=cname, value=land.setup[,get(cname)], pos=1) 
  })
  
  # Compute a few convenience quantities based on problem setup
  num.regions = nrow(land.setup)  
  
  # Assemble matrix of cumulative development (if unchecked).
  # This will be derived from either a single fixed rate per ecoregion OR
  # a vector of rates per ecoregion, in which case def.rate is a list (not unlike cost)
  if(class(def.rate)=="list") {
    cum.devel.mat = apply(matrix(unlist(def.rate), nrow=num.regions, ncol=num.periods, byrow=TRUE),
                          1,
                          cumsum)
  } else {
    cum.devel.mat = apply(matrix(def.rate, nrow=num.regions, ncol=num.periods, byrow=FALSE), 
                          1, 
                          cumsum)
  }
  # convert to matrix for special case of 1 period problem since otherwise gets 
  # retyped as numeric            
  if(num.periods == 1) {
    cum.devel.mat = matrix(cum.devel.mat, nrow=num.periods, ncol=num.regions)
  }


  # Move budget to maximize objective, subject to feasibility
  # Guiding principle: 
  # - iteratively move budget from low to high marginal ROI 
  
  
  # Initial allocation equal across regions in each period
  cur.budget = c(1/num.regions) * matrix(budget.per.period, nrow = num.regions, ncol=num.periods, byrow=T) #Set initial proportion of budget allocation
  
  # Iniitalize other numerical tracking variables
  roi.diff = conv.tol*10
  num.pass.adjustments = 1
  last.race.over.time =ncol(cur.budget)
  fwd.pass = 1
  total.budget.chg = 100000
  # Loop until we've converged at a solution. 
  # Start with the initial budget guess, then repeatedly
  # 1. Identify the ending conditions, including ROI, that follow from that budget
  # 2. Sweep forward through time, reallocating budget from lower ROI to higher ROI areas.
  # 3. Check for convergence
  
  # Details of reallocation and convergence
  # Consider ROIs if we add a dollar vs remove a dollar:
  # - If we add a dollar to the budget for any region, its (right) ROI should be lower than the highest
  #   (left) ROI from any region where we could potentially take a dollar away.
  # - If a region is fully protected, the right ROI will be zero, which will be below the left ROI for any potential donor region.
  # - If a region is fully conserved but reforestation is possible, the right ROI will be the lower reforestation ROI, while
  #   its left ROI will be the conservation ROI
  # - If a region has no spend, its right ROI and left ROI should both be quite low. The right ROI will satisfy the condition
  #   for not wanting more spend in that region, while the low left ROI will be irrelevant since it won't be a potential donor region.
  #
  # For the reallocation procedure, we'll use this logic. Once we identify the lowest left (source) ROI and highest right (target) ROI,
  # we'll start with a candidate allocation of ALL of the source region's ROI. If the target right ROI is still higher than the source
  # left ROI after that, do the full allocation. We'll then do a binary search to find the smallest 
  # reallocation we can that makes the target's right ROI smaller than the source's left ROI. The idea is that when we consider
  # the target region as a source, we want its left ROI to be higher than other potential target region ROIs, and this choice makes
  # that most likely. We'll of course have to iterate, but it favors convergence.
  
  # Go until we make no adjustments. We're checking in each period for the ROI-based convergence criterion. If we make any reallocations,
  # it failed, so keep going. If we make none, then the ROI-based convergence holds in all periods.
  while(num.pass.adjustments>0 && total.budget.chg > budget.conv.tol & fwd.pass<=max.fwd.pass) {
    last.source = last.target = -1
    recent.target = rep(F, num.regions)
    recent.source = rep(F, num.regions)
    t = 1
    # set up copy of budget for convergence checks. We'll look at the difference
    # across the entire time sweep.
    last.budget = cur.budget
    num.pass.adjustments = 0
    total.budget.chg = 0
    target.ignore.ctr = rep(0, nrow(cur.budget))
    while(t <= last.race.over.time) {
      print(t)
      # find race over times and ROIs for current allocation.
      # we'll actually do this twice: once considering a small increase in budget everywhere (for right ROIs), 
      # and once with a small decrease in budget everywhere (for left ROIs). This should help deal with 
      # numerical precision issues, especially around discontinuities.
      right.budget = left.budget = cur.budget
      right.budget[,t] = right.budget[,t] + budget.nudge
      left.budget[,t] = pmax(left.budget[,t] - budget.nudge,0)

      race.over.info = FindRaceEndConditionsProbabilistic(cur.budget, land.setup, cum.devel.mat, profiles, profile.abundance, cost.el, restoration.value, t)
      
      race.over.info.right = FindRaceEndConditionsProbabilistic(right.budget, land.setup, cum.devel.mat, profiles, profile.abundance, cost.el, restoration.value, t)
      race.over.info.left = FindRaceEndConditionsProbabilistic(left.budget, land.setup, cum.devel.mat, profiles, profile.abundance, cost.el, restoration.value, t)
    
      last.race.over.time = max(race.over.info$race.end.periods)
      
      rois.right = race.over.info.right$rois.right
      rois.left = race.over.info.left$rois.left
      period.budgets = cur.budget[,t]
      if(class(ce)=="list") {
        period.ces = sapply(ce, function(rce) { rce[t]})
      } else {
        period.ces = ce
      }
      
      # - calculate how much budget would be required to purchase remaining land in last period. To do this
      #   we need to know starting and ending conservation in that period, not just the difference
      # - subtract from actual budget
      if(fixed.mc) {
        cost.of.rem.land = race.over.info$area.avail.raceendper*period.ces
      } else {
        # isoelastic cost. Note that if we allow restoration, the cost of available land
        # (including land that could be restored) will have the second term be
        # exactly zero. To avoid any numerical precision issues, we can just handle this
        # case explicitly.
        if(cost.el!=-1) {
          ep.oneplusep = cost.el/(1+cost.el)
          if(restoration.value>0) {
            cost.of.rem.land = ep.oneplusep*period.ces*((area.forestable.init-race.over.info$area.resrest.raceendper.start)^(1/ep.oneplusep)) 
          } else {
            cost.of.rem.land = ep.oneplusep*period.ces*((area.forestable.init-race.over.info$area.resrest.raceendper.start)^(1/ep.oneplusep) - (area.forestable.init-(race.over.info$area.resrest.raceendper.start+race.over.info$area.avail.raceendper))^(1/ep.oneplusep)) 
          }
        } else {
          stop("Cost elasticity of -1 not implemented")
        }
      }

      # Find out what extra budget we have: beyond whats needed to purchase land remaining in the 
      # final period (which may not be this period -- we include that condition below)
      excess.budget = period.budgets - cost.of.rem.land
      
      # The main ROI code doesn't set ROIs to zero when the race is over (because in earlier periods
      # we want to consider it to still be positive due to discreteness of time periods).
      # Set it to zero for regions in which this is the final period and we have more than enough
      # $ to conserve everything.
      rois.right[race.over.info.right$race.end.periods==t & excess.budget>=0] = 0
      rois.left[race.over.info.left$race.end.periods==t & excess.budget>0] = 0
      rois.right[race.over.info.right$race.end.periods<t] = 0
      rois.left[race.over.info.left$race.end.periods<t] = 0
        
      
      # Identify source region. Use left ROIs, and only consider regions that
      # have budget to give.
      roi.feas = rois.left
      roi.feas[period.budgets==0] = NA
      source = which.min(roi.feas)
      source.roi = roi.feas[source]
 
      # Identify target region. Use right ROIs, and only consider regions where
      # the race is not already over as of this period given current budgets.
      
      roi.feas = rois.right
      roi.feas[race.over.info.right$race.end.periods<t] = NA
      roi.feas[race.over.info.right$race.end.periods==t & cost.of.rem.land <= period.budgets] = NA
      
      # We'll also temporarily ignore any region that acted as an intermediary, where we transferred
      # budget in and immediately back out. That might happen if ROI changes much more rapidly in an area (e.g.,
      # almost out of area to reforest where MC rises sharply) than others, so that swapping ROI ranking happens
      # without much change to the source ROI, so that the target becomes the new source.
      roi.feas[target.ignore.ctr>0] = NA
      target = which.max(roi.feas)
      target.roi = roi.feas[target]
      
      # there exists a beneficial reallocation in this period and the ROIs aren't already sufficiently close
      # to be considered converged, reallocate.
      if(length(source)>0 && length(target)>0 && source!=target && source.roi < target.roi && (target.roi-source.roi)>conv.tol) {
        num.pass.adjustments = num.pass.adjustments + 1
        print(paste0("Finding reallocation for ", source, "->", target, ": ", source.roi, ", ",target.roi))
        # start by considering reallocating ALL of source regions budget.
        realloc.amt = cur.budget[source, t]
 
        if((excess.budget[target]<0) && ((-excess.budget[target])< realloc.amt)) {
          realloc.amt = -excess.budget[target]
        }
        
        # Now we'll start a binary search between this initial realloc.amt and zero 
        # for the smallest reallocation that makes the target's right ROI smaller than
        # the source's left ROI.
        last.realloc.delta = init.realloc.amt =  smallest.flip.realloc = realloc.amt
        
        while(T) {
          trial.budget = cur.budget
          trial.budget[source, t] = cur.budget[source, t] - realloc.amt
          trial.budget[target, t] = cur.budget[target, t] + realloc.amt
          right.trial.budget = left.trial.budget = trial.budget
          right.trial.budget[,t] = right.trial.budget[,t] + budget.nudge
          left.trial.budget[,t] = pmax(left.trial.budget[,t] - budget.nudge,0)
          

          trial.race.over.info =  FindRaceEndConditionsProbabilistic(trial.budget, land.setup, cum.devel.mat, profiles, profile.abundance, cost.el, restoration.value, t)
          trial.race.over.info.left =  FindRaceEndConditionsProbabilistic(left.trial.budget, land.setup, cum.devel.mat, profiles, profile.abundance, cost.el, restoration.value, t)
          trial.race.over.info.right =  FindRaceEndConditionsProbabilistic(right.trial.budget, land.setup, cum.devel.mat, profiles, profile.abundance, cost.el, restoration.value, t)
          
          trial.rois.left = trial.race.over.info.left$rois.left
          trial.rois.right = trial.race.over.info.right$rois.right
          
          trial.period.budgets = trial.budget[,t]
          
          # - calculate how much budget would be required to purchase remaining land in last period. To do this
          #   we need to know starting and ending conservation in that period, not just the difference
          # - subtract from actual budget
          if(fixed.mc) {
            trial.cost.of.rem.land = trial.race.over.info$area.avail.raceendper*period.ces
          } else {
            if(cost.el!=-1) {
              ep.oneplusep = cost.el/(1+cost.el)
              trial.cost.of.rem.land = ep.oneplusep*period.ces*((area.forestable.init-trial.race.over.info$area.resrest.raceendper.start)^(1/ep.oneplusep) - (area.forestable.init-(trial.race.over.info$area.resrest.raceendper.start+trial.race.over.info$area.avail.raceendper))^(1/ep.oneplusep)) 
            } else {
              stop("Cost elasticity of -1 not implemented")
            }
          }
          
          # Find out what extra budget we have: beyond whats needed to purchase remaining land in the 
          # final period (which may not be this period -- we include that condition below)
          trial.excess.budget = trial.period.budgets - trial.cost.of.rem.land
          
          # The main ROI code doesn't set ROIs to zero when the race is over (because in earlier periods
          # we want to consider it to still be positive due to discreteness of time periods).
          # Set it to zero for regions in which this is the final period and we have more than enough
          # $ to conserve everything.
          trial.rois.right[trial.race.over.info.right$race.end.periods==t & trial.excess.budget>=0] = 0
          trial.rois.left[trial.race.over.info.left$race.end.periods==t & trial.excess.budget>0] = 0
          trial.rois.right[trial.race.over.info.right$race.end.periods<t] = 0
          trial.rois.left[trial.race.over.info.left$race.end.periods<t] = 0
          
          trial.roi.target = trial.rois.right[target]
          trial.roi.source = trial.rois.left[source]
          
          if((trial.roi.target>trial.roi.source)) {
            if(realloc.amt == init.realloc.amt) {
              # liquidating budget of source leaves target with higher ROI still. Make reallocation
              break
            } 
            if((last.realloc.delta<=realloc.tol)) {
              # We've reached the reallocation tolerance, but the current reallocation doesn't bring target ROI below source.
              # revert to smallest reallocation we considered that did make target ROI smaller than source ROI.
              realloc.amt = smallest.flip.realloc
              break
            } else {
              # we must be here because we decreased the reallocation, and that made the target ROI higher than the source ROI.
              # Go back the other way.
              realloc.amt = realloc.amt + last.realloc.delta/2
              last.realloc.delta = last.realloc.delta/2
            }
          }  
          if(trial.roi.target<=trial.roi.source) { # we got target ROI down below source ROI
            if(realloc.amt < smallest.flip.realloc) {
              smallest.flip.realloc = realloc.amt
            }
            if((last.realloc.delta<=realloc.tol)) {
              # binary search converged on a reallocation. Even if we don't meet criteria for ROI convergence
              # because of numerical precision & ROI discontinuities, stop.& reallocate
              break
            }
            # What we're really after: shift that makes the target ROI 
            # lower than the source ROI AFTER the transfer
            # close to the source ROI AFTER the transfer
            # higher than the source ROI BEFORE the transfer (so we get contraction of the ROI span)
            # AND
            # source ROI AFTER the transfer is no higher than target ROI BEFORE the transfer.
            # This suggests a contraction of the ROIs among regions having spend, which should
            # push us toward convergence.
            if(((trial.roi.source-trial.roi.target) <= inner.conv.tol*realloc.amt) && 
                 trial.roi.target >= source.roi &&
                 trial.roi.source <= target.roi
               ) {
              # Break and reallocate.
              break
            } else {
              # We flipped the ranking of the ROIs, but they're not very close. Reallocate less.
              realloc.amt = realloc.amt - last.realloc.delta/2
              last.realloc.delta = last.realloc.delta/2
            }
          }          
        }
        
        print(paste("fp: ", fwd.pass, "; t=", t, ": realloc", source, target, realloc.amt, source.roi, target.roi))
        cur.budget[source, t] = cur.budget[source, t] - realloc.amt
        cur.budget[target, t] = cur.budget[target, t] + realloc.amt
        recent.target[target] = T
        recent.source[source] = T
        total.budget.chg = total.budget.chg + realloc.amt
        
        target.ignore.ctr = pmax(target.ignore.ctr-1, 0)
        if(source == last.target) {
          target.ignore.ctr[source] = NIGNOREITER
        }
        
        
        last.source = source 
        last.target = target
      
      } else {
        # If this test fails, no sufficiently large, profitable reallocation possible in this period. 
        # Go to the next time step
        t = t+1
        recent.target = rep(F, num.regions)
        recent.source = rep(F, num.regions)
      }
      
    }
    fwd.pass = fwd.pass + 1
  }
  
  # Return a list with
  # 1. Solution budget matrix: rows are regions and columns are periods
  # 2. Race end conditions
  prelim.res = FindRaceEndConditionsProbabilistic(cur.budget, land.setup, cum.devel.mat, profiles, profile.abundance, cost.el, restoration.value, 1)
  final.rois = sapply(1:num.regions, FUN=function(r) {
    r.res = FindRaceEndConditionsProbabilistic(cur.budget, land.setup, cum.devel.mat, profiles, profile.abundance, cost.el, restoration.value, prelim.res$race.end.periods[r])
    return(c(r.res$rois.left[r], r.res$rois.right[r]))
  })
  prelim.res$rois.left = final.rois[1,]
  prelim.res$rois.right = final.rois[2,]
  global.soln = list(budget=cur.budget,
                     race.end=prelim.res,
                     land.setup = land.setup,
                     profiles = profiles,
                     profile.abundance = profile.abundance,
                     cost.el = cost.el,
                     restoration.value = restoration.value)
  if(is.null(out.file)) {
    return(global.soln)
  } else {
    save(global.soln, file=out.file)
  }
}