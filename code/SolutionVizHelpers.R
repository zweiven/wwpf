source(here("code","SummarizeResults.R"))
CombineSolutionData = function(sol.filename, merged.data, probabilistic=F, calc.traj=F) {
  load(sol.filename)
  loc.merged.data = data.table(copy(merged.data))

  budget.per.period = as.numeric(gsub(x=sol.filename, pattern=".*p_(.*)b_.*", replacement="\\1")) 
  
  if(is.null(global.soln$race.end$end.rois)) {
    if(is.null(global.soln$race.end$rois)) {
      global.soln$race.end$end.rois=global.soln$race.end$rois.left
    } else {
      global.soln$race.end$end.rois=global.soln$race.end$rois
    }
  }
  
  if(is.null(global.soln$race.end$end.periods)) {
    global.soln$race.end$end.periods=global.soln$race.end$race.end.periods
  }
  #tack on solution data to problem setup.
  
  loc.merged.data[!is.na(forest_gdp),`:=`(end.area.conserved=global.soln$race.end$end.areas.conserved,
                                          end.area.restored=global.soln$race.end$end.areas.restored,
                                          end.area.cons.or.rest=global.soln$race.end$end.areas.conserved+global.soln$race.end$end.areas.restored,
                                          end.roi=global.soln$race.end$end.rois,
                                          cons.end.period=global.soln$race.end$cons.end.periods,
                                          end.period=global.soln$race.end$end.periods,
                                          end.spend=rowSums(global.soln$budget),
                                          first5.spend=rowSums(global.soln$budget[,1:5]),
                                          first10.spend=rowSums(global.soln$budget[,1:10]),
                                          first50.spend=rowSums(global.soln$budget[,1:50])
                                          )]
  
  if(probabilistic) {
    prot.stats = GetFinalProtectionStats(global.soln)
    loc.merged.data[!is.na(forest_gdp), `:=`(exp.prot.spp.by.ecoregion.new=prot.stats$exp.prot.spp.by.ecoregion.new,
                                             exp.add.prot.spp.by.ecoregion.new=prot.stats$exp.add.prot.spp.by.ecoregion.new,
                                             initial.mb=prot.stats$initial.mb,
                                             initial.mc=prot.stats$initial.mc,
                                             initial.mb.mc.ratio=prot.stats$initial.mb.mc.ratio)]
    if(calc.traj) {
      loc.ecoregion.data = fread(here("data","alldata.csv"))
      loc.ecoregion.data[,forestable.area.init:=forest_area_2018+18*dozer_rate_2000_2018]
      loc.ecoregion.data=merge(loc.ecoregion.data, global.soln$land.setup[,.(eco_code, z=ze)], by="eco_code")
      loc.ecoregion.data[,ae:=exp(log(plant_spcs) - z*log(forestable.area.init))]
      soln.traj = ComputeTrajectories(budget = global.soln$budget, 
                                      ecoregion.data = loc.ecoregion.data, 
                                      profiles = global.soln$profiles, 
                                      profile.abundance = global.soln$profile.abundance,
                                      cost.el= global.soln$cost.el, 
                                      restoration.value = global.soln$restoration.value)
    }
  }
  
  # subtract off starting protected area to get newly protected area
  loc.merged.data[,new.area.conserved := end.area.conserved-protected_forest_area_2018]
  
  loc.merged.data[,new.area.cons.or.rest := end.area.cons.or.rest-protected_forest_area_2018]
  
  
  # Compute # of new species conserved in each ecoregion
  loc.merged.data[,new.species.conserved := ae*(end.area.conserved^z - protected_forest_area_2018^z)]
  loc.merged.data[,forest.info:= ifelse(!is.na(forest_area_2018), "Forested Ecoregion", "Not Forested Ecoregion")]
  
  
  loc.merged.data[,unprotected_forest_2018 := forest_area_2018 - protected_forest_area_2018]
  loc.merged.data[,threat:=dozer_rate_2000_2018/unprotected_forest_2018]
  loc.merged.data[,initial.species.conserved := ae*protected_forest_area_2018^z]
  
  loc.merged.data[,area.conserved.5:=protected_forest_area_2018 + first5.spend/GDPperm2_partial]
  loc.merged.data[,new.species.conserved.5 :=ae*(area.conserved.5^z - protected_forest_area_2018^z)]
  
  loc.merged.data[,area.conserved.10:=protected_forest_area_2018 + first10.spend/GDPperm2_partial]
  loc.merged.data[,new.species.conserved.10 :=ae*(area.conserved.10^z - protected_forest_area_2018^z)]
  
  loc.merged.data[,area.conserved.50:=protected_forest_area_2018 + first50.spend/GDPperm2_partial]
  loc.merged.data[,new.species.conserved.50 :=ae*(area.conserved.50^z - protected_forest_area_2018^z)]
  
  if(calc.traj) {
    return(list(soln.summary=loc.merged.data,
                soln.traj=soln.traj))
  } else {
    return(list(soln.summary=loc.merged.data))
  }

}

# Visualize
MakeMap = function(solution.dt, 
                   stat.to.show="new.species.conserved.10", 
                   stat.label = "Species protected \nin first 10 years",
                   png.path=NULL, pdf.path=NULL) {
  map.theme =   theme_minimal() + 
    theme(
      legend.justification=c(0,0),
      legend.position=c(0.01, 0.1),
      legend.key.height = unit(0.065, units="npc"),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_line(colour = "transparent"))
  
  overlaid.map.raceend = ggplot(solution.dt) + 
    geom_sf(aes(fill=forest.info, geometry=geometry), 
            colour="#dddddd") + 
    coord_sf(expand=F) + 
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0), 
                       limits = c(-56,78)) + # crop out Antarctica, extra arctic ocean, & part of Greenland (irrelevant)
    geom_path(data=coastlines.fortified, 
              aes(x=long, y=lat, group=group), 
              colour="#cccccc") +
    geom_point(data=solution.dt[solution.dt[,get(stat.to.show)]>0 & solution.dt$end.period<=10,],
               aes(x=lon, y=lat, size=get(stat.to.show), fill="End year <= 10"),
               pch=21) +
    geom_point(data=solution.dt[solution.dt[,get(stat.to.show)]>0 & solution.dt$end.period>10,],
               aes(x=lon, y=lat, size=get(stat.to.show), fill="End year > 10"),
               pch=21) +
    map.theme +
    scale_size_continuous(name=stat.label, range=c(1,10)) + 
    scale_fill_manual(name="",
                      values=c("Not Forested Ecoregion"=NA, 
                               "End year <= 10"="#228B22", 
                               "End year > 10"="#024B72", 
                               "Forested Ecoregion"="#efefef"),
                      labels=c("End year <= 10"=expression("End year" <= 10), 
                               "End year > 10"="End year > 10", 
                               "Forested Ecoregion"="Forested Ecoregion"),
                      limits=c("End year <= 10","End year > 10", "Forested Ecoregion")) +
    theme(legend.position = c(0.01,0.01),
          legend.spacing.y=unit(0.01, "cm"),
          legend.text.align = 0,
          legend.key = element_rect(color="white")) + 
    guides(fill = guide_legend(order = 2, override.aes=list(shape=c(NA,NA,NA))), 
           size = guide_legend(order = 1))
  
  if(!is.null(png.path)) {
    png(png.path, width=1100, height=500)
    overlaid.map.raceend
    dev.off()
  }
  if(!is.null(pdf.path)) {
    pdf(pdf.path, width=11, height=5)
    overlaid.map.raceend
    dev.off()
  }
  
  return(overlaid.map.raceend)
}


# Visualize
MakeMapShaded = function(solution.dt, 
                         stat.size="new.species.conserved", 
                         stat.shade="frac.conserved10",
                         stat.linecol=NULL,
                         size.label = "Additional species protected",
                         shade.label = "% protected in first 10 years",
                         stat.display.thresh=0.5, # due to tolerance there can be small protection levels that are really zeros. Filter for display.
                         png.path=NULL, pdf.path=NULL) {
  map.theme =   theme_minimal() + 
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_line(colour = "transparent"))
  
  overlaid.map.raceend = ggplot(solution.dt) + 
    geom_sf(data=solution.dt[forest.info=="Forested Ecoregion",],
            mapping=aes(geometry=geometry), 
            colour="#dddddd", fill="#cccccc") + 
    geom_sf(data=solution.dt[forest.info=="Not Forested Ecoregion",],
            mapping=aes(geometry=geometry), 
            colour="#dddddd", fill="#ffffff") + 
    coord_sf(expand=F) + 
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0), 
                       limits = c(-56,78)) + # crop out Antarctica, extra arctic ocean, & part of Greenland (irrelevant)
    geom_path(data=coastlines.fortified, 
              aes(x=long, y=lat, group=group), 
              colour="#cccccc")
  if(is.null(stat.linecol)) {
    overlaid.map.raceend = overlaid.map.raceend + geom_point(data=solution.dt[solution.dt[,get(stat.size)]>0 & solution.dt$end.spend>0,],
                                                             aes(x=lon, y=lat, 
                                                                 size=get(stat.size), 
                                                                 fill=get(stat.shade)),
                                                             pch=21)
  } else {
    overlaid.map.raceend = overlaid.map.raceend + geom_point(data=solution.dt[solution.dt[,get(stat.size)]>0 & solution.dt$end.spend>0,],
                                                             aes(x=lon, y=lat, 
                                                                 size=get(stat.size), 
                                                                 fill=get(stat.shade),
                                                                 color=get(stat.linecol)),
                                                             pch=21, stroke=1.25) + 
      scale_color_manual(values=c("#000000","#FFDA00"), guide="none")
  }
  overlaid.map.raceend = overlaid.map.raceend +
    map.theme +
    scale_size_continuous(name=size.label, range=c(1,10), breaks=c(1,100,1000)) + 
    theme(legend.position = "bottom",
          legend.key = element_rect(color="white")
          )
  
  if(class(solution.dt[,get(stat.shade)]) %in% c("character","factor")) {
    overlaid.map.raceend = overlaid.map.raceend+
    scale_fill_manual(name=shade.label,
                      breaks=c(F,T),
                      values=c("#ffffff","#014421"),
                      labels=c("0%","90-100%")) + 
      guides(fill = guide_legend(order = 2, 
                                 override.aes=list(size=c(4,4))), 
             size = guide_legend(order = 1))
  } else {
    overlaid.map.raceend = overlaid.map.raceend + 
      scale_fill_gradient(name=shade.label, low="#014421", high="#ffffff") +
      guides(size = guide_legend(order = 1),
             fill = guide_colorbar(title.vjust=0.8))
  }
  if(!is.null(png.path)) {
    png(png.path, width=1100, height=500)
    overlaid.map.raceend
    dev.off()
  }
  if(!is.null(pdf.path)) {
    pdf(pdf.path, width=11, height=5)
    overlaid.map.raceend
    dev.off()
  }
  
  return(overlaid.map.raceend)
}


MakeStatComparisonMap = function(solution.ref.dt, 
                                 solution.new.dt, 
                                 stat.to.show, 
                                 stat.label, terms="abs",
                                 stat.display.thresh=0, # due to tolerance there can be small differences that are immaterial (e.g., near zero). Can filter for display
                                 size.breaks = c(10000,100000,250000),
                                 size.labels = c("10,000", "100,000", "250,000"),
                                 size.limits = c(1,500000),
                                 diff.scale = 1E6,
                                 png.path=NULL, pdf.path=NULL) {
  map.dt = merge(solution.ref.dt, 
                 solution.new.dt[,.(ECO_CODE, new.stat=get(stat.to.show))],
                 by="ECO_CODE")
  map.dt[,old.stat:=get(stat.to.show)]
  
  if(terms=="abs") {
    map.dt[,stat.chg:=(new.stat-old.stat)/diff.scale]
  } else if(terms=="perc") {
    map.dt[,stat.chg:=(new.stat-old.stat)/old.stat]
  } else if(terms=="ratio") {
    map.dt[,stat.chg:=new.stat/old.stat]
  } else if(terms=="absbinned") {
    map.dt[,rawdiff:=new.stat-old.stat]
    map.dt[,stat.chg:=floor(log10(abs(rawdiff)))]
    map.dt[stat.chg<0,stat.chg:=0]
    map.dt[rawdiff<0,stat.chg:=-stat.chg]
  }
  map.dt[,incr:=stat.chg>=0]
  
  map.theme =   theme_minimal() + 
    theme(
      legend.key.height = unit(0.065, units="npc"),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_line(colour = "transparent"))
  
  
  overlaid.map.raceend = ggplot(map.dt) + 
    geom_sf(aes(
      geometry=geometry), 
      colour="#dddddd", fill="#ffffff") + 
    coord_sf(expand=F) + 
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0), 
                       limits = c(-56,78)) + # crop out Antarctica, extra arctic ocean, & part of Greenland (irrelevant)
    geom_path(data=coastlines.fortified, 
              aes(x=long, y=lat, group=group), 
              colour="#cccccc") +
    geom_point(data=map.dt[!is.na(incr) & abs(stat.chg)>stat.display.thresh,],   
               aes(x=lon, y=lat, 
                   size=abs(stat.chg),
                   color=incr,
                   #fill=incr,
                   shape=incr)) +
    map.theme +
    scale_size_continuous(name=stat.label,
                          breaks=size.breaks,
                          labels=size.labels,
                          limits=size.limits,
                          #trans="log",
                          range=c(1,6)) + 
    scale_shape_manual(name="Change",
                       values=c("TRUE"=24,
                                "FALSE"=25),
                       labels=c("TRUE"="Increase",
                                "FALSE"="Decrease")
    ) + 
    scale_color_manual(name="Change",
                       values=c("TRUE"="black",
                                "FALSE"="red"),
                       labels=c("TRUE"="Increase",
                                "FALSE"="Decrease"))+
    scale_fill_manual(name="Change",
                      values=c("TRUE"="black",
                               "FALSE"="red"),
                      labels=c("TRUE"="Increase",
                               "FALSE"="Decrease"))+
    theme(legend.position = "bottom",
          legend.key = element_rect(color="white")) +
    guides(fill = guide_legend(order = 2),#, override.aes=list(shape=c(24,25))), 
           size = guide_legend(order = 1, override.aes = list(shape=24)))
  
  if(!is.null(png.path)) {
    png(png.path, width=1100, height=500)
    overlaid.map.raceend
    dev.off()
  }
  if(!is.null(pdf.path)) {
    pdf(pdf.path, width=11, height=5)
    overlaid.map.raceend
    dev.off()
  }
  
  return(overlaid.map.raceend)
  
}

ComputeTrajectories = function(budget, 
                               ecoregion.data, 
                               profiles=NULL,
                               profile.abundance=NULL,
                               cost.el=-6, 
                               restoration.value=0.8) {
  ep.oneplusep = cost.el/(1+cost.el)
  budget.dt.wide = data.table(budget)
  setnames(budget.dt.wide, as.character(1:ncol(budget.dt.wide)))
  annot.budget.dt.wide = cbind(ecoregion.data[,.(wwf_mhtnam, wwf_realm2, eco_code, eco_num, eco_name, 
                                                 def.rate=dozer_rate_2000_2018,
                                                 init.mc=GDPperm2_partial, 
                                                 area.reserved.init=protected_forest_area_2018, 
                                                 area.forested.init=forest_area_2018,
                                                 area.forestable.init=forest_area_2018+18*dozer_rate_2000_2018,
                                                 ae, z)], budget.dt.wide)
  
  annot.budget.dt.long = melt(data=annot.budget.dt.wide, 
                              id.vars = c("wwf_mhtnam", "wwf_realm2", "eco_code", "eco_num", "eco_name", "def.rate",
                                          "init.mc","area.reserved.init", "area.forested.init", "area.forestable.init","ae","z"),
                              variable.name = "year",
                              value.name = "budget.allocated")
  setorder(annot.budget.dt.long, eco_code,year)
  annot.budget.dt.long[,ce:=init.mc/((area.forestable.init-area.reserved.init)^(1/cost.el))]
  # Now go through year by year and calculate
  # area reserved by end of year
  # area restored by end of year
  # prob protection by end of year
  # total spp by end of year
  annot.budget.dt.long[year==1,area.reserved.cur:=area.reserved.init]
  annot.budget.dt.long[year==1,area.restored.cur:=0]
  annot.budget.dt.long[year==1,area.purchased.cur:=area.reserved.cur+area.restored.cur]
  
  annot.budget.dt.long[year==1,area.forestable.available.cur:=area.forestable.init-area.reserved.init]
  annot.budget.dt.long[year==1,area.forested.available.cur:=area.forested.init-area.reserved.init]
  
  # calculate initial protection stats
  # Calculate spp protection stats at end of year
  annot.budget.dt.long[year==1,p.protect := ((area.reserved.cur + restoration.value*area.restored.cur)/area.forestable.init)^z]
  annot.budget.dt.long[year==1 & p.protect>1,p.protect:=1]
  init.prob.prot.by.region = profiles * matrix(annot.budget.dt.long[year==1,p.protect], ncol=nrow(budget), nrow=length(profile.abundance), byrow=T)
  # this gives a n.profiles x n.regions matrix of probability of protecting that profile in that region
  # now we construct the probability that each profile is protected
  init.prob.profile.prot = apply(init.prob.prot.by.region, 1, FUN=function(profile.prot.probs) { 
    1 - (prod(1-profile.prot.probs))
  })
  # calculate total protection 
  init.exp.prot.spp = sum(init.prob.profile.prot*profile.abundance)
  
  init.region.exp.prot=colSums(init.prob.prot.by.region * matrix(profile.abundance, nrow=length(profile.abundance), ncol=nrow(budget), byrow=F))
  
  for(y in 1:ncol(budget)) {
    annot.budget.dt.long[year==y,init.exp.prot.spp.by.ecoregion := colSums(init.prob.prot.by.region * matrix(profile.abundance, nrow=length(profile.abundance), ncol=nrow(budget), byrow=F))]
  }
  tot.exp.newly.prot.spp = rep(NA, ncol(budget))
  tot.exp.prot.spp = rep(NA, ncol(budget))
  for(y in 1:ncol(budget)) {
    # calculate increments: how much is purchased, broken down into reserved & restored
    annot.budget.dt.long[year==y, land.newly.purchased := area.forestable.available.cur - (pmax(0,area.forestable.available.cur^(1/ep.oneplusep) - budget.allocated/ce/ep.oneplusep))^ep.oneplusep] 
    annot.budget.dt.long[year==y & budget.allocated<=0, land.newly.purchased :=0] # fix rounding issues 
    annot.budget.dt.long[year==y, land.newly.reserved := pmin(land.newly.purchased, area.forested.available.cur)]
    annot.budget.dt.long[year==y, land.newly.restored := land.newly.purchased-land.newly.reserved]
    
    # Update tallies of forestable area available, forested area available, 
    # area purchased, broken down into reserved & restored
    annot.budget.dt.long[year==y,area.reserved.cur:=area.reserved.cur + land.newly.reserved]
    annot.budget.dt.long[year==y,area.restored.cur:=area.restored.cur + land.newly.restored]
    annot.budget.dt.long[year==y,area.purchased.cur:=area.purchased.cur + land.newly.purchased]
    annot.budget.dt.long[year==y,area.forestable.available.cur:=area.forestable.available.cur-land.newly.purchased]
    annot.budget.dt.long[year==y,area.forested.available.cur:=pmax(0,area.forested.available.cur-land.newly.reserved-def.rate)]
    
    # Propagate those to following year as starting values
    if(y<ncol(budget)) {
      annot.budget.dt.long[year==(y+1),area.reserved.cur:=annot.budget.dt.long[year==y,area.reserved.cur]]
      annot.budget.dt.long[year==(y+1),area.restored.cur:=annot.budget.dt.long[year==y,area.restored.cur]]
      annot.budget.dt.long[year==(y+1),area.purchased.cur:=annot.budget.dt.long[year==y,area.purchased.cur]]
      annot.budget.dt.long[year==(y+1),area.forestable.available.cur:=annot.budget.dt.long[year==y,area.forestable.available.cur]]
      annot.budget.dt.long[year==(y+1),area.forested.available.cur:=annot.budget.dt.long[year==y,area.forested.available.cur]]
    }
    
    # Calculate spp protection stats at end of year
    annot.budget.dt.long[year==y,p.protect := ((area.reserved.cur + restoration.value*area.restored.cur)/area.forestable.init)^z]
    annot.budget.dt.long[year==y & p.protect>1,p.protect:=1]
    prob.prot.by.region = profiles * matrix(annot.budget.dt.long[year==y,p.protect], ncol=nrow(budget), nrow=length(profile.abundance), byrow=T)
    # this gives a n.profiles x n.regions matrix of probability of protecting that profile in that region
    # now we construct the probability that each profile is protected
    prob.profile.prot = apply(prob.prot.by.region, 1, FUN=function(profile.prot.probs) { 
      1 - (prod(1-profile.prot.probs))
    })
    # calculate total protection 
    tot.exp.prot.spp[y] = sum(prob.profile.prot*profile.abundance)
    tot.exp.newly.prot.spp[y] = tot.exp.prot.spp[y]-init.exp.prot.spp
    annot.budget.dt.long[year==y,exp.prot.spp.by.ecoregion:=colSums(prob.prot.by.region * matrix(profile.abundance, nrow=length(profile.abundance), ncol=nrow(budget), byrow=F))]
    annot.budget.dt.long[year==y,exp.newly.prot.spp.by.ecoregion:=exp.prot.spp.by.ecoregion-init.exp.prot.spp.by.ecoregion]
  }
  
  return(list(annot.budget.dt.long=annot.budget.dt.long,
              tot.exp.prot.spp=tot.exp.prot.spp,
              tot.exp.newly.prot.spp=tot.exp.newly.prot.spp))
}





kMapNACol = "#e5e5e5"
kBorderCol = "#e5e5e5"
coastlines.fortified = fortify(coastline110)
# Load spatial data
load(here("data","ecoregions_simplified.Rdata"))

# Load ecoregion data
ecoregion.data = fread(here("data","alldata.csv"))
ecoregion.data[,GDPperm2:=forest_gdp/forest_area_2018]
ecoregion.data[,GDPperm2_partial:=forest_partial_gdp/forest_area_2018]
ecoregion.data[,z:=0.2]
ecoregion.data[,ae:=exp(log(plant_spcs) - z*log(forest_area_2018))]

# Apply same filtering used in solution code
ecoregion.data = ecoregion.data[!is.na(GDPperm2) & GDPperm2>0,]

# and unprotected forest info
ecoregion.data[,unprotected_forest_2018 := forest_area_2018 - protected_forest_area_2018]

sim.periods = 50
ecoregion.data[,initial.sl.dozers:=ae*unprotected_forest_2018^z - ae*(pmax(unprotected_forest_2018-dozer_rate_2000_2018,0))^z]
ecoregion.data[,overall.sl.dozers:=ae*unprotected_forest_2018^z - ae*(pmax(unprotected_forest_2018-sim.periods*dozer_rate_2000_2018,0))^z]



# merge with spatial data
merged.data = merge(ecoregions.simple, ecoregion.data, by.x="ECO_CODE", by.y="eco_code", all.x=T)

# Add a measure of threat:
# Deforestation rate doesn't fully capture the race, since the severity of the race
# depends also upon how much land is unprotected. An intuitive idea is # of years
# it would take for bulldozers to consume all unprotected forest. To turn that
# into a "threat index" for which higher #s indicate higher threat, we can simply
# take the inverse, i.e. (deforestation rate)/(unprotected forest) rather than
# (unprotected forest)/(deforestation rate).
merged.data$unprotected_forest_2018 = merged.data$forest_area_2018 - merged.data$protected_forest_area_2018
merged.data$threat = merged.data$dozer_rate_2000_2018/merged.data$unprotected_forest_2018

threat.quants = quantile(merged.data$threat, probs=seq(from=0,to=1,by=0.2),na.rm=T)
merged.data$threat.quantile = cut(merged.data$threat, breaks=threat.quants)

# SF doesn't like the TNC shapefile; revert to earlier processing
# just for the mapping of centroids. Any issues 
# with projections affect visualization in a minor way only.
sf::sf_use_s2(FALSE)
centroids = st_coordinates(st_centroid(merged.data$geometry))
sf::sf_use_s2(TRUE)
merged.data$lon = centroids[,1]
merged.data$lat = centroids[,2]



