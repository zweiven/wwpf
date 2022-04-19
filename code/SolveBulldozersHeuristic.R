# Function to solve a forest conservation problem
# using one of several myopic heuristics.
# Arguments: 
# problem.setup: data.table with info per ecoregion. Rows correspond to ecoregions; columns have characteristics of each ecoregion.
# num.periods: number of periods in optimization (50 in base case in paper)
# budget.per.period: max budget that can be spent per period
# profiles: matrix (rows: profiles, columns: ecoregions) of presence/absence data for species profiles
# profile.abundance: vector of abundances of different species profiles
# cost.el: cost elasticity of land (as a proxy for more general protection costs)
# restoration.value: multiplier indicating relative value of restored forest vs. intact forest for species protection
# heur.type: text indicating which of several heuristics is to be used. 
#            Options are:
#            hotspot: most unprotected species remaining in ecoregion
#            mb: highest marginal benefit
#            mc: lowest marginal cost
#            bcratio: highest marginal benefit to marginal cost ratio
#            threat: highest fraction of remaining forested land deforested per year
#            percprot: lowest percent of forestable land protected

# budget.incr: how large are the budget increments that are to be allocated sequentially as part of the heuristic solution?
# out.file: name of filename to save results to. If null, results simply returned to calling code.

SolveBulldozersProblemHeuristic =  function(problem.setup,
                                            num.periods,
                                            budget.per.period,
                                            profiles, 
                                            profile.abundance,
                                            cost.el=-999999999,
                                            restoration.value = 0,
                                            heur.type="hotspot",
                                            budget.incr = NA,
                                            out.file=NULL) {
  # Copy problem info for local manipulation
  land.setup = copy(problem.setup)
  num.regions = nrow(land.setup)
  
  # Problem initialization
  if(is.na(budget.incr)) {
    budget.incr = budget.per.period/1000
  }
  fixed.mc = is.na(cost.el)
  land.setup[,forestable.area.available.init:= area.forestable.init-area.reserved.init] # Initially Available Area  
  land.setup[,forestable.area.available.cur:=forestable.area.available.init]
  land.setup[,forested.area.available.init:= area.forested.init-area.reserved.init] # Initially Available Area
  land.setup[,forested.area.available.cur:=forested.area.available.init]
  land.setup[,area.reserved.cur:=area.reserved.init]
  land.setup[,area.restored.cur:=0]
  land.setup[,ze:=log(species/ae)/log(area.forestable.init)] # species area curve exponent
  
  if(restoration.value>0) {
    land.setup[,area.available.cur:=forestable.area.available.init]
  } else {
    land.setup[,area.available.cur:=forested.area.available.init]
  }
  land.setup[,orig.id:=1:nrow(land.setup)]
  land.setup[,threat:=def.rate/(forested.area.available.init-area.reserved.init)] #1/time to bulldoze initially unprotected forest
  land.setup[,percprot:=area.reserved.cur/area.forestable.init]
  

  # calculate initial metrics involving marginal cost and marginal benefit
  if(fixed.mc) {
    land.setup[,mc := ce]
  } else {
    # MC is a function of forestable land that is not yet protected via either conservation or reforestation
    land.setup[, mc := ce*(area.forestable.init-area.reserved.cur)^(1/cost.el)]
    # Note this is not well defined if all available land has been conserved via protection or reforestation, in which case 
    # MC should be infinite and no additional land can be conserved (we correct res.totals earlier).
    # Shouldn't happen, but want to handle edge case in case.
    land.setup[area.forestable.init<area.reserved.cur, mc:=Inf] 
  }

  p.protect = land.setup[,((area.reserved.cur + restoration.value*area.restored.cur)/area.forestable.init)^ze]
  p.protect[p.protect>1]=1
  
  land.setup[,exp.spp.unprotected:=(1-p.protect)*species]
  
  
  # Marginal probability of protection will depend upon whether current purchase
  # levels have put us into a restoration case. For initialization, that is simply
  # any area with <= 0 available forested land. For iteration below, we'll have to calculate
  # whether we're in a restoration case based on budget allocation.
  marg.p.protect = land.setup[,ze*((area.reserved.cur + restoration.value * area.restored.cur)^(ze-1))/(area.forestable.init^ze)]
  in.restoration = land.setup[,forested.area.available.cur<=0]
  marg.p.protect[in.restoration] = marg.p.protect[in.restoration]*restoration.value
  
  marg.p.terms = profiles * matrix(marg.p.protect, ncol=nrow(land.setup), nrow=length(profile.abundance), byrow=T)
  marg.p.terms[is.nan(marg.p.terms)] = 0
  unprot.terms = apply(profiles, 1, function(sp.type.pa.vec) { (1-p.protect*sp.type.pa.vec)})
  fullprods = colprods(unprot.terms)
  land.setup[,mb := sapply(1:num.regions, function(i) {
    sum(profile.abundance*marg.p.terms[,i]*(fullprods/unprot.terms[i,]))
  })]
  
  land.setup[,bcratio:=mb/mc]

  
  land.setup[,exp.spp.unprotected:=(1-p.protect)*species]
  
  # For tracking solution (keep identical components to proposed solution)
  budget = matrix(0, nrow=nrow(land.setup), ncol=num.periods)
  cons.end.periods = rep(num.periods,nrow(land.setup))
  rest.end.periods = rep(num.periods,nrow(land.setup))
  
  # for each period:
  # - Set remaining budget = budget for each period
  # - repeat the following until remaining budget = 0 
  # 1. Take a budget increment b=min($X, remaining budget), where $X is a fixed increment
  # 2. Rank (decreasing) regions by the specified heuristic and choose the top ranked area
  # 3. Find the maximum spendable budget in the area identified in the previous step (2)
  # 4. Update the budget increment to be b=min(b, maximum spendable budget in identified area)
  # 5. Spend b in the area identified in step 2 
  # 6. Update the remaining budget, protected area, etc.
  
  for(i in 1:num.periods) {
    print(paste0("Period: ", i))
    cur.budget.rem = budget.per.period
    while(cur.budget.rem > 0) {
      cur.budget.incr = min(budget.incr, cur.budget.rem)
      
      # Based on the supplied heuristic, find the top-ranked ecoregion
      # in which there is still positive area available for purchase
      if(heur.type=="hotspot") {
        # most unprotected species first - ignore complementarity
        target.region.id = land.setup[area.available.cur>0,][which.max(exp.spp.unprotected),orig.id]
      } else if(heur.type=="bcratio") {
        # highest marginal benefit to marginal cost ratio first
        target.region.id = land.setup[area.available.cur>0,][which.max(bcratio),orig.id]
      } else if(heur.type=="mb") {
        # highest marginal benefit first
        target.region.id = land.setup[area.available.cur>0,][which.max(mb),orig.id]
      } else if(heur.type=="threat") {
        # most threatened first
        target.region.id = land.setup[area.available.cur>0,][which.max(threat),orig.id]
      } else if(heur.type=="percprot") {
        # lowest % protected first
        target.region.id = land.setup[area.available.cur>0,][which.min(percprot),orig.id]
      } else if(heur.type=="cost") {
        # lowest marginal cost first
        target.region.id = land.setup[area.available.cur>0,][which.min(mc),orig.id]
      }
      
      
      # Find the highest priority area where land can still be purchased
      target.area = land.setup[orig.id==target.region.id,]
      target.rownum = which(land.setup$orig.id==target.region.id)
      target.forestable.available = target.area[["forestable.area.available.cur"]]
      target.forested.available = target.area[["forested.area.available.cur"]]
      target.reserved = target.area[["area.reserved.cur"]]
      target.restored = target.area[["area.restored.cur"]]
      target.purchased = target.reserved + target.restored
      target.available = target.area[["area.available.cur"]]
      target.ce = target.area[["ce"]]

      # Calculate cost of purchasing all available land
      if(fixed.mc) {
        cost.to.purchase.all.land = target.ce*target.available
      } else {
        if(cost.el!=-1) {
          ep.oneplusep = cost.el/(1+cost.el)
          if(restoration.value>0) {
            cost.to.purchase.all.land = ep.oneplusep*target.ce*(target.forestable.available^(1/ep.oneplusep)) 
          } else {
            cost.to.purchase.all.land = ep.oneplusep*target.ce*(target.forestable.available^(1/ep.oneplusep) - (target.forestable.available-target.forested.available)^(1/ep.oneplusep)) 
          }
        } else {
          stop("Cost elasticity of -1 not supported")
        }
      }
      purchasing.all.land = cost.to.purchase.all.land<=cur.budget.incr
      cur.budget.incr = min(cur.budget.incr, cost.to.purchase.all.land)
      
      # Calculate land purchased using that budget increment

      if(!is.na(cost.el)) {
        if(purchasing.all.land) { # avoid rounding issues
          if(restoration.value>0) {
            land.newly.purchased = target.forestable.available
          } else {
            land.newly.purchased = target.forested.available
          }
        } else {
          land.newly.purchased = target.forestable.available - (target.forestable.available^(1/ep.oneplusep) - cur.budget.incr/target.ce/ep.oneplusep)^ep.oneplusep 
        }
      } else {
        land.newly.purchased = cur.budget.incr/target.ce
      }

      # Update land bookkeeping 

      land.newly.reserved = max(0, min(land.newly.purchased, target.forested.available-target.purchased))
      land.newly.restored = land.newly.purchased-land.newly.reserved

      land.setup[target.rownum,
                 `:=`(area.reserved.cur = area.reserved.cur + land.newly.reserved,
                      area.restored.cur = area.restored.cur + land.newly.restored,
                      forested.area.available.cur = max(forested.area.available.cur - land.newly.reserved, 0),
                      forestable.area.available.cur = max(forestable.area.available.cur - land.newly.purchased, 0),
                      area.available.cur=ifelse(restoration.value>0,forestable.area.available.cur,forested.area.available.cur)
                      )]
      
      # Update threat index to reflect reforestation
      if(heur.type=="threat") {
        land.setup[target.rownum, threat:=def.rate/forested.area.available.cur]
      }
      
      # Update MC for target area
      if(!fixed.mc) {
        # MC is a function of forestable land that is not yet protected via either conservation or reforestation
        land.setup[target.rownum, 
                   mc := ce*(area.forestable.init-(area.reserved.cur+area.restored.cur))^(1/cost.el)]
        
        # Note this is not well defined if all available land has been conserved via protection or reforestation, in which case 
        # MC should be infinite and no additional land can be conserved (we correct res.totals earlier).
        # Shouldn't happen, but want to handle edge case in case.
        if(land.setup[target.rownum, forestable.area.available.cur <= 0]) {
          land.setup[target.rownum, mc:=Inf]
        } 
      }
      
      # Update MBs and MB/MC ratio for ALL areas
      p.protect = land.setup[,((area.reserved.cur + restoration.value*area.restored.cur)/area.forestable.init)^ze]
      p.protect[p.protect>1]=1
      # Marginal probability of protection will depend upon whether current purchase
      # levels have put us into a restoration case. For initialization, that is simply
      # any area with <= 0 available forested land. For iteration below, we'll have to calculate
      # whether we're in a restoration case based on budget allocation.
      marg.p.protect = land.setup[,ze*((area.reserved.cur + restoration.value * area.restored.cur)^(ze-1))/(area.forestable.init^ze)]
      in.restoration = land.setup[,forested.area.available.cur<=0]
      marg.p.protect[in.restoration] = marg.p.protect[in.restoration]*restoration.value
      
      marg.p.terms = profiles * matrix(marg.p.protect, ncol=nrow(land.setup), nrow=length(profile.abundance), byrow=T)
      marg.p.terms[is.nan(marg.p.terms)] = 0
      unprot.terms = apply(profiles, 1, function(sp.type.pa.vec) { (1-p.protect*sp.type.pa.vec)})
      fullprods = colprods(unprot.terms)
      land.setup[, mb := sapply(1:num.regions, function(i) {
        sum(profile.abundance*marg.p.terms[,i]*(fullprods/unprot.terms[i,]))
      })]
      land.setup[, bcratio := mb/mc]
      
      land.setup[is.infinite(mc),bcratio:=0]
      
      # update expected species unprotected, ignoring complementarity 
      if(heur.type=="hotspot"){
        land.setup[target.rownum,
                   exp.spp.unprotected := (1-p.protect[target.rownum])*species]
      }
      
      
      # update percent protected
      if(heur.type=="percprot") {
        land.setup[target.rownum,
                   percprot:=(area.reserved.cur+area.restored.cur)/area.forestable.init]
      }
      
      # - Budget
      cur.budget.rem = cur.budget.rem - cur.budget.incr
      budget[target.region.id,i] = budget[target.region.id,i] + cur.budget.incr
    }
    
    # deforestation happens in all areas
    land.setup[,forested.area.available.cur:=pmax(0,forested.area.available.cur-def.rate)]

    # update threat index
    if(heur.type=="threat") {
      land.setup[target.rownum, threat:=def.rate/forested.area.available.cur]
    }
    
    # if conservation opportunities end because forest is either reserved or 
    # deforested, note that. Find all regions with no available land,
    # and overwrite with current period if it is smaller than the end period that's currently there.
    # The latter condition avoids overwriting end periods for races that ended earlier.
    # 
    # We do something analogous for tracking when (if) we run out of land to reforest too
    cons.ended.regions = land.setup[forested.area.available.cur<=0,orig.id]
    cons.end.periods[cons.ended.regions]=pmin(cons.end.periods[cons.ended.regions],i)

    rest.ended.regions = land.setup[forestable.area.available.cur<=0,orig.id]
    rest.end.periods[rest.ended.regions]=pmin(rest.end.periods[rest.ended.regions],i)
    
  }
  race.end = list(cons.end.periods = cons.end.periods,
                  rest.end.periods = rest.end.periods,
                  end.areas.conserved = land.setup[,area.reserved.cur],
                  end.areas.restored = land.setup[,area.restored.cur],
                  end.rois = rep(NA, nrow(land.setup)))
  heur.soln = list(budget=budget,
                   race.end=race.end)
  if(!is.null(out.file)){
    save(heur.soln, file=out.file)
  } else {
    return(heur.soln)      
  }

}
