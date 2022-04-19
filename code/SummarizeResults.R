GetFinalProtectionStats = function(soln) {
  cur.budget = soln$budget 
  num.regions = nrow(cur.budget)
  num.periods = ncol(cur.budget)
  land.setup = soln$land.setup
  
  # Assemble matrix of cumulative development (if unchecked).
  # This will be derived from either a single fixed rate per ecoregion OR
  # a vector of rates per ecoregion, in which case def.rate is a list (not unlike cost)
  if(class(land.setup$def.rate)=="list") {
    cum.devel.mat = apply(matrix(unlist(land.setup$def.rate), nrow=num.regions, ncol=num.periods, byrow=TRUE),
                          1,
                          cumsum)
  } else {
    cum.devel.mat = apply(matrix(land.setup$def.rate, nrow=num.regions, ncol=num.periods, byrow=FALSE), 
                          1, 
                          cumsum)
  }

  # convert to matrix for special case of 1 period problem since otherwise gets 
  # retyped as numeric            
  if(num.periods == 1) {
    cum.devel.mat = matrix(cum.devel.mat, nrow=num.periods, ncol=num.regions)
  }
  profiles = soln$profiles
  profile.abundance = soln$profile.abundance 
  num.profiles = length(profile.abundance)
  cost.el = soln$cost.el
  restoration.value = soln$restoration.value
  anti.profiles = 1 - profiles
  fixed.mc = is.na(cost.el)
  land.setup[,area.restored.init:=0]
  with(land.setup,
       {
         if(class(ce)=="list") {
           ce.mat = matrix(unlist(ce), ncol=num.regions, nrow=num.periods, byrow=F)
         } else {
           ce.mat = matrix(ce, ncol=num.regions, nrow=num.periods, byrow=T)
         }
         
         # Calculate and store initial mb, mc, mb:mc ratio
         if(fixed.mc) {
           land.setup[,initial.mc := ce.mat[1,]]
         } else {
           # MC is a function of forestable land that is not yet protected via either conservation or reforestation
           land.setup[, initial.mc := ce.mat[1,]*(area.forestable.init-area.reserved.init)^(1/cost.el)]
           # Note this is not well defined if all available land has been conserved via protection or reforestation, in which case 
           # MC should be infinite and no additional land can be conserved (we correct res.totals earlier).
           # Shouldn't happen, but want to handle edge case in case.
           land.setup[area.forestable.init<area.reserved.init, initial.mc:=Inf] 
         }
         
         p.protect = ((area.reserved.init + restoration.value*area.restored.init)/area.forestable.init)^ze
         p.protect[p.protect>1]=1
         
         exp.spp.unprotected=(1-p.protect)*species
         
         
         # Marginal probability of protection will depend upon whether current purchase
         # levels have put us into a restoration case. For initialization, that is simply
         # any area with <= 0 available forested land. For iteration below, we'll have to calculate
         # whether we're in a restoration case based on budget allocation.
         marg.p.protect = ze*((area.reserved.init + restoration.value * area.restored.init)^(ze-1))/(area.forestable.init^ze)
         in.restoration = forested.area.available.init<=0
         marg.p.protect[in.restoration] = marg.p.protect[in.restoration]*restoration.value
         
         marg.p.terms = profiles * matrix(marg.p.protect, ncol=nrow(land.setup), nrow=length(profile.abundance), byrow=T)
         marg.p.terms[is.nan(marg.p.terms)] = 0
         unprot.terms = apply(profiles, 1, function(sp.type.pa.vec) { (1-p.protect*sp.type.pa.vec)})
         fullprods = colprods(unprot.terms)
         land.setup[,initial.mb := sapply(1:num.regions, function(i) {
           sum(profile.abundance*marg.p.terms[,i]*(fullprods/unprot.terms[i,]))
         })]
         
         
         # Find end period. This is the first period in which, given rates of deforestation and protection,
         # if both sides did not collide, the available land at the end of that period would be non-positive.
         
         # Cumulative land bought per current allocation at end of each period
         # - calculate cumulative spend as of period t, then compute land purchased using that vector
         cum.spend = apply(cur.budget, 1, cumsum)
         initial.forestable.mat = matrix(area.forestable.init, ncol=num.regions, nrow=num.periods, byrow=T)
         initial.protected.mat = matrix(area.reserved.init, ncol=num.regions, nrow=num.periods, byrow=T)
         
         # calculate cumulative *new* area reserved by end of each period, 
         # temporarily ignoring running into bulldozers
         if(fixed.mc) {
           res.totals =  apply(cur.budget/t(ce.mat), 1, cumsum) 
         } else {
           if(cost.el!=-1) {
             ep.oneplusep = cost.el/(1+cost.el)
             # numerical approximations can lead to land conserved appearing to be slightly larger than is available, leading to NAs
             inner.term = pmax((initial.forestable.mat-initial.protected.mat)^(1/ep.oneplusep) - cum.spend/(ce.mat*ep.oneplusep),0)
             res.totals = initial.forestable.mat - initial.protected.mat - (inner.term)^ep.oneplusep 
           } else {
             stop("Cost elasticity of -1 not implemented")
           }
         }
         
         # deal with R automatically retyping as numeric for single period case
         if(num.periods==1) {
           res.totals = matrix(res.totals, ncol=num.regions, nrow=num.periods)
         }
         rem.land = forested.area.available.init - t(res.totals + cum.devel.mat) # How much available land is remaining at end of each period given current plan
         
         # Find out when the 'race' is over in each region. Earlier of:
         # a) when conservation meets development: pristine forest runs out by the end of that period
         # b) the end of the planning horizon
         meeting.periods = apply(rem.land, 1, FUN=function(rem.vec) { min(which(rem.vec<=0)) })
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
             area.avail.consendper = rem.land[r, cons.end.per-1]
           } else {
             area.avail.consendper = forested.area.available.init[r]
           }
           return(area.avail.consendper)
         })
         res.totals.at.consend = area.reserved.init +  # land initially conserved
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
         
         # If restoration is allowed, we need to know conditions for when the restoration phase ends.
         # It will be the earlier of 
         # - when all forested land is either conserved or restored
         # - the end of the planning horizon
         # Restoration begins when conservation stops.
         # As specified, note restoration costs the same as protection, but is just less effective per 
         # the restoration.value multiplier. This means we can use protected land totals computed earlier, but just
         # need to know how much is conserved and how much is restored. 
         if(restoration.value > 0) {
           area.newly.conserved = res.totals.at.consend-area.reserved.init
           area.restored = res.totals-area.newly.conserved
           area.restored[area.restored<0] = 0
           area.restored[is.nan(area.restored)] = 0
           # for restoration, there is no development to consider, so we only need to constrain
           # restoration totals by the total amount of forestable land initially available
           # (the sum of conservation and restoration can't exceed initially forestable land)
           area.avail.for.rest = matrix(area.forestable.init - res.totals.at.consend,
                                        ncol=num.regions,
                                        nrow=num.periods, 
                                        byrow=T)
           
           # can't restore more than is available for restoration
           area.restored[area.restored>area.avail.for.rest] = area.avail.for.rest[area.restored>area.avail.for.rest]
           
           # Find out when restoration 'finishes' (if it does): when
           # area available for restoration runs out.
           rest.end.periods = sapply(1:num.regions, FUN=function(r) {
             rest.end = min(which(area.restored[,r]>=area.avail.for.rest[r]))
           })
           # If we don't actually run out of land to restore, 
           # the 'end' period for restoration is the end of the planning horizon
           rest.end.periods = pmin(rest.end.periods, num.periods)
           
           # Find out ending area restored per region
           rest.totals.at.restend = apply(area.restored, 2, max)
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
         
         
         p.protect.init = (area.reserved.init/area.forestable.init)^ze      

         # for the sake of summary, we'll focus on potential spend at the end of the planning horizon:
         roi.period = num.periods
         
         if(restoration.value > 0) {
           # marginal benefits depend upon whether the marginal unit of area affected by the marginal 
           # unit of budget in the specified period is conserved or restored
           marg.p.protect = ze*((res.totals.at.consend + restoration.value * rest.totals.at.restend)^(ze-1))/(area.forestable.init^ze)
           # if the roi.period is later than the cons.end.period for a particular region OR
           # if the roi.period is the same as the cons.end.period but the budget is more than enough to conserve
           # available undeveloped, forested land in that period, 
           # then the marginal dollar would go toward restoration and we need to multiply the marginal.p.protect by the restoration.value
           roi.period.rest = (roi.period>cons.end.periods) | (roi.period==cons.end.periods & (res.totals[roi.period,]>res.totals.at.consend)) 
           marg.p.protect[roi.period.rest] = marg.p.protect[roi.period.rest]*restoration.value
         } else {
           marg.p.protect = ze*(res.totals.at.consend^(ze-1))/(area.forestable.init^ze)
         }

         # And build up terms indicating probability of each species type being unprotected in each ecoregion.
         # this will be a matrix with a column per species type, a row per ecoregion,
         # and the entry representing the probability the species of that type is unprotected
         # in that ecoregion.
         marg.p.terms = profiles * matrix(marg.p.protect, ncol=nrow(land.setup), nrow=length(profile.abundance), byrow=T)
         unprot.terms = apply(profiles, 1, function(sp.type.pa.vec) { (1-p.protect*sp.type.pa.vec)})
         
         unprot.terms.init = apply(profiles, 1, function(sp.type.pa.vec) { (1-p.protect.init*sp.type.pa.vec)})
         
         # now for each ecoregion x species type, we want the probability that
         # additional protected area would newly protect that species. We then
         # sum that across species in an ecoregion to get MB, and divide by MC to get ROI.
         
         # we can take full products of unprotected terms across all ecoregions once since that's expensive
         # then in the region-specific loop, we divide by the unprotected term from that focal region to remove its effect.
         fullprods = colprods(unprot.terms)
         fullprods.init = colprods(unprot.terms.init)
         mbs = sapply(1:num.regions, function(i) {
           sum(profile.abundance*marg.p.terms[,i]*(fullprods/unprot.terms[i,]))
         })
         
         # For each species, the probability of protection is one minus the probability it is
         # not protected anywhere. The probability it is not protected anywhere is the product
         # of probabilities of not being protected in each ecoregion. The probability of not being
         # protected in an ecoregion is one minus the probability it is protected there. To be
         # protected it must both be present and protected, so we multiply the P/A value by the
         # probability of protection in that region.
         prob.prot.by.region = profiles * matrix(p.protect, ncol=nrow(land.setup), nrow=length(profile.abundance), byrow=T)
         # this gives a n.profiles x n.regions matrix of probability of protecting that profile in that region
         # now we construct the probability that each profile is protected
         prob.profile.prot = apply(prob.prot.by.region, 1, FUN=function(profile.prot.probs) { 
           1 - (prod(1-profile.prot.probs))
         })
         
         prob.prot.by.region.init = profiles * matrix(p.protect.init, ncol=nrow(land.setup), nrow=length(profile.abundance), byrow=T)
         # this gives a n.profiles x n.regions matrix of probability of protecting that profile in that region
         # now we construct the probability that each profile is protected
         prob.profile.prot.init = apply(prob.prot.by.region.init, 1, FUN=function(profile.prot.probs) { 
           1 - (prod(1-profile.prot.probs))
         })
         

         exp.prot.spp = sum(prob.profile.prot*profile.abundance)
         exp.prot.spp.init = sum(prob.profile.prot.init*profile.abundance)
         
         
         # Might also want to know 
         # 1. Expected # species protected in an ecoregion
         exp.prot.spp.by.ecoregion = colSums(prob.prot.by.region * matrix(profile.abundance, nrow=length(profile.abundance), ncol=nrow(land.setup), byrow=F))
         exp.prot.spp.by.ecoregion.init = colSums(prob.prot.by.region.init * matrix(profile.abundance, nrow=length(profile.abundance), ncol=nrow(land.setup), byrow=F))
         
         # 2. Expected additional species protected by ecoregion
         #    This will be # species of a type in ecoregion times probability of protection of that type (profile) in that ecoregion and nowhere else
         #    Analogous to MB calcs but with probability of protection in focal region
         #    rather than marginal probability.
         
         exp.add.prot.spp.by.ecoregion = sapply(1:num.regions, function(i) {
           sum(profile.abundance*prob.prot.by.region[,i]*(fullprods/unprot.terms[i,]))
         })

         exp.add.prot.spp.by.ecoregion.init = sapply(1:num.regions, function(i) {
           sum(profile.abundance*prob.prot.by.region.init[,i]*(fullprods.init/unprot.terms.init[i,]))
         })        
         
         # Note we are really interested in how many *newly* protected species there 
         # are in expectation. So we should do these calculations 
         # for initial protection levels too, then subtract.
 
         # We'll now compute marginal costs given allocations up to and including roi.period
         if(fixed.mc) {
           period.mcs = ce
         } else {
           roi.period.res.totals = res.totals[roi.period,]
           roi.period.res.totals[cons.end.periods==roi.period] =  (res.totals.at.consend-area.reserved.init)[cons.end.periods==roi.period]
           period.mcs = ce.mat[roi.period,]*(area.forestable.init-roi.period.res.totals)^(1/cost.el)
         }
         rois = mbs/period.mcs
         
         # In the edge case in which both ze and res.totals.at.race.over are zero, this will be NaN.
         # If no land is conserved AND conserving is useless (such that ze=0), then the ROI must 
         # be zero since spending can't increase forested area conserved
         rois[is.nan(rois)] = 0
         
         # give back info about ending conditions, ROI, etc.
         if(restoration.value > 0) {
           race.end.periods = rest.end.periods
         } else {
           race.end.periods = cons.end.periods
         }
         
         # Finally, calculate two things about the final period during which 
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
         # If restoration is allowed, this is basically any initially forested land that is
         # not already conserved.
         # If restoration is not allowed, this is any undeveloped initially forested land that is not already conserved
         
         area.avail.raceendper = sapply(1:num.regions, function(r) {
           # find available land at start of period in which race ends.
           race.end.per = race.end.periods[r]
           if(race.end.per > 1) {
             if(restoration.value > 0) {
               area.avail.raceendper = (forested.area.available.init[r] - res.totals[race.end.per-1,r])
             } else {
               area.avail.raceendper = rem.land[r, race.end.per-1]
             }
           } else {
             area.avail.raceendper = forested.area.available.init[r]
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
                     rois = rois,
                     p.protect = p.protect,
                     exp.prot.spp = exp.prot.spp,
                     exp.prot.spp.by.ecoregion=exp.prot.spp.by.ecoregion,
                     exp.add.prot.spp.by.ecoregion=exp.add.prot.spp.by.ecoregion,
                     exp.prot.spp.init = exp.prot.spp.init,
                     exp.prot.spp.by.ecoregion.init=exp.prot.spp.by.ecoregion.init,
                     exp.add.prot.spp.by.ecoregion.init=exp.add.prot.spp.by.ecoregion.init,
                     exp.prot.spp.new = exp.prot.spp-exp.prot.spp.init,
                     exp.prot.spp.by.ecoregion.new=exp.prot.spp.by.ecoregion - exp.prot.spp.by.ecoregion.init,
                     exp.add.prot.spp.by.ecoregion.new=exp.add.prot.spp.by.ecoregion-exp.add.prot.spp.by.ecoregion.init,
                     initial.mb=land.setup$initial.mb,
                     initial.mc=land.setup$initial.mc,
                     initial.mb.mc.ratio=land.setup[,initial.mb/initial.mc]
                     )) 
         
       })
}
