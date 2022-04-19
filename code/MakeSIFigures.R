# Code to make figures in supplementary information
library(sf)
library(data.table)
library(ggplot2)
library(cowplot)
library(rmapshaper)
library(rnaturalearthdata)
library(Rfast)
source(here("code","SolutionVizHelpers.R"))

profile.data = fread(here("data","mammal_based_profiles.csv"))
profile.abundance = profile.data[,filt.abundances]               
profiles = as.matrix(profile.data[,-ncol(profile.data), with=F])

#########
# Alternate budget plot
#########
ComputeSppConserved = function(budget) {
  load(here("results",paste0("altbudget_",budget,".Rdata")))
  budget.traj = ComputeTrajectories(global.soln$budget, ecoregion.data, profiles=profiles,profile.abundance = profile.abundance)
  budget.traj$tot.exp.newly.prot.spp[50]
}
budgets = c(1E7,5E7,1E8,2.5E8,5E8,7.5E8,1E9,1E10,2.5E10)
budget.text = c("1e+07","5e+07","1e+08","2.5e+08","5e+08","7.5e+08","1e+09","1e+10", "2.5e+10")
budget.res.dt = data.table(budget=budgets,
                           new.exp.spp.conserved = sapply(budget.text, FUN=ComputeSppConserved))

budget.exp.plot = ggplot(budget.res.dt, 
                         aes(x=budget/1E6, 
                             y=new.exp.spp.conserved))  + 
  geom_point() + 
  geom_line() +
  xlab("Annual budget (Million $)") +
  ylab("Expected species newly conserved") +
  ylim(0,max(budget.res.dt$new.exp.spp.conserved)) + 
  theme_minimal()

pdf(here("figures","SI_BudgetExpansion.pdf"), width=4, height=4)
budget.exp.plot
dev.off()

# Bound cost savings from using optimal solution vs heuristics

# Load main solution
load(here("results","main.Rdata"))

# Load heuristic solutions
load(here("results","heur_bcratio.Rdata"))
heur.soln.bcratio = heur.soln

load(here("results","heur_hotspot.Rdata"))
heur.soln.hotspot = heur.soln

load(here("results","heur_mb.Rdata"))
heur.soln.mb = heur.soln

load(here("results","heur_threat.Rdata"))
heur.soln.threat = heur.soln

load(here("results","heur_cost.Rdata"))
heur.soln.cost = heur.soln

load(here("results","heur_percprot.Rdata"))
heur.soln.percprot = heur.soln
rm(heur.soln)

# Compute trajectories of (expected) species conservation
soln.trajectories = ComputeTrajectories(global.soln$budget, ecoregion.data, profiles=profiles,profile.abundance = profile.abundance)
hotspot.trajectories = ComputeTrajectories(heur.soln.hotspot$budget, ecoregion.data, profiles=profiles,profile.abundance = profile.abundance)
bcratio.trajectories = ComputeTrajectories(heur.soln.bcratio$budget, ecoregion.data, profiles=profiles,profile.abundance = profile.abundance)
mb.trajectories = ComputeTrajectories(heur.soln.mb$budget, ecoregion.data, profiles=profiles,profile.abundance = profile.abundance)
threat.trajectories = ComputeTrajectories(heur.soln.threat$budget, ecoregion.data, profiles=profiles,profile.abundance = profile.abundance)
cost.trajectories = ComputeTrajectories(heur.soln.cost$budget, ecoregion.data, profiles=profiles,profile.abundance = profile.abundance)
percprot.trajectories = ComputeTrajectories(heur.soln.percprot$budget, ecoregion.data, profiles=profiles,profile.abundance = profile.abundance)

new.spp.data = rbind(data.table(soln="Optimal Solution",year=1:50,tot.protected.spp.in.year=soln.trajectories$tot.exp.newly.prot.spp),
                     data.table(soln="Species Richness Heuristic",year=1:50,tot.protected.spp.in.year=hotspot.trajectories$tot.exp.newly.prot.spp),
                     data.table(soln="Benefit:Cost Heuristic",year=1:50,tot.protected.spp.in.year=bcratio.trajectories$tot.exp.newly.prot.spp),
                     data.table(soln="Marginal Benefit Heuristic",year=1:50,tot.protected.spp.in.year=mb.trajectories$tot.exp.newly.prot.spp),
                     data.table(soln="Threat Heuristic",year=1:50,tot.protected.spp.in.year=threat.trajectories$tot.exp.newly.prot.spp),
                     data.table(soln="Cost Heuristic",year=1:50,tot.protected.spp.in.year=cost.trajectories$tot.exp.newly.prot.spp), 
                     data.table(soln="% Protection Heuristic",year=1:50,tot.protected.spp.in.year=percprot.trajectories$tot.exp.newly.prot.spp))
new.spp.data[,soln:=factor(soln, levels=c("Optimal Solution",
                                          "Benefit:Cost Heuristic",
                                          "Marginal Benefit Heuristic",
                                          "Cost Heuristic",
                                          "% Protection Heuristic",
                                          "Species Richness Heuristic",
                                          "Threat Heuristic"))]


# find minimal budget from optimal solution that dominates heuristic
budget.indices = sapply(new.spp.data[year==50,tot.protected.spp.in.year], 
                        FUN=function(spp.prot) {
                          cut(x=spp.prot, 
                              breaks=budget.res.dt$new.exp.spp.conserved, 
                              labels=F)
                        })
# we'll express savings in terms of the smallest budget we consider for the 
# optimal solution that would outperform the heuristic (i.e., a lower bound on savings)
min.dominant.budget = budget.res.dt[budget.indices+1,budget]
min.dominant.budget[is.na(min.dominant.budget)] = budget.res.dt$budget[1]
savings.dt = new.spp.data[year==50,.(soln, tot.protected.spp.in.year)][,min.dominant.budget:=min.dominant.budget]
savings.dt[,min.percent.savings:=(1E9-min.dominant.budget)/1E9]

savings.dt

#######
# Series of two-panel figures comparing the solution under baseline assumptions 
# to that under alternate assumptions.
#######
# Load baseline solution information for comparison
soln.base = CombineSolutionData(here("results","main.Rdata"), 
                                merged.data,
                                probabilistic = T)

# Function to create two-panel map
MakeSIAltAssumptionFig = function(alt.soln.fname, 
                                  soln.base, 
                                  main.stat.col="exp.prot.spp.by.ecoregion.new",
                                  calc.traj=T) {
  base.plot.data = soln.base$soln.summary
  soln.alt = CombineSolutionData(alt.soln.fname, 
                                 merged.data,
                                 probabilistic = T,
                                 calc.traj=calc.traj)

  # Identify areas in which conservation of intact forest
  # is done in period 1
  alt.plot.data = soln.alt$soln.summary
  alt.plot.data[,cons.end.1p:=cons.end.period==1]
  
  if(calc.traj==T) {
    yr.50.prot = soln.alt$soln.traj$annot.budget.dt.long[year==50,
                                                         .(eco_code, 
                                                           "exp.prot.spp.by.ecoregion.new.50"=exp.newly.prot.spp.by.ecoregion)]
    alt.plot.data = merge(alt.plot.data, yr.50.prot, by.x="ECO_CODE", by.y="eco_code") 
    base.plot.data[,exp.prot.spp.by.ecoregion.new.50:=exp.prot.spp.by.ecoregion.new]
  }
  setorder(alt.plot.data, -cons.end.period)
  new.spp.map = MakeMapShaded(solution.dt = alt.plot.data, 
                              stat.size = main.stat.col,
                              stat.shade = "cons.end.period",
                              stat.linecol= "cons.end.1p",
                              shade.label = "Year last intact\n forest protected",
                              size.label = "Expected species\n newly protected")
  
  change.map = MakeStatComparisonMap(base.plot.data,
                                     alt.plot.data, 
                                     stat.to.show = main.stat.col, 
                                     stat.label = "Change in expected \nspecies conserved",
                                     terms="abs",
                                     stat.display.thresh = 1E-5,
                                     size.breaks = c(1,10,100),
                                     size.labels = c("1", "10", "100"),
                                     size.limits = c(1E-5,1000),
                                     diff.scale = 1
  ) 
  comb.plot = plot_grid(new.spp.map, 
                        change.map, 
                        labels=c("A","B"), 
                        ncol=1, nrow=2)
  return(list(comb.plot = comb.plot, soln.alt=soln.alt, alt.plot.data=alt.plot.data))
}

##################
# Alternate land cost
##################
altcost.plot = MakeSIAltAssumptionFig(here("results","altcost.Rdata"),
                                      soln.base)

pdf(here("figures","SI_altcost.pdf"), width=11, height=11)
altcost.plot$comb.plot
dev.off()




##################
# Different SAR exponents (z values)
##################
altz_0.1plot = MakeSIAltAssumptionFig(here("results","altz_0.1.Rdata"),
                                   soln.base)

pdf(here("figures","SI_altz_0.1.pdf"), width=11, height=11)
altz_0.1plot$comb.plot
dev.off()


altz_0.3plot = MakeSIAltAssumptionFig(here("results","altz_0.3.Rdata"),
                                      soln.base)

pdf(here("figures","SI_altz_0.3.pdf"), width=11, height=11)
altz_0.3plot$comb.plot
dev.off()

altz_storchplot = MakeSIAltAssumptionFig(here("results","altz_storch.Rdata"),
                                      soln.base)

pdf(here("figures","SI_altz_storch.pdf"), width=11, height=11)
altz_storchplot$comb.plot
dev.off()

altz_kierplot = MakeSIAltAssumptionFig(here("results","altz_kier.Rdata"),
                                         soln.base)

pdf(here("figures","SI_altz_kier.pdf"), width=11, height=11)
altz_kierplot$comb.plot
dev.off()


##################
# Different forest loss data
##################
# Hansen *net* forest loss
hansen.net.plot = MakeSIAltAssumptionFig(here("results","altdef_hansennet.Rdata"),
                                         soln.base)

pdf(here("figures","SI_hansennet.pdf"), width=11, height=11)
hansen.net.plot$comb.plot
dev.off()

# Copernicus forest loss
copernicus.plot = MakeSIAltAssumptionFig(here("results","altdef_copernicus.Rdata"),
                                         soln.base)

pdf(here("figures","SI_copernicus.pdf"), width=11, height=11)
copernicus.plot$comb.plot
dev.off()

# SSP-3-based forest loss
ssp.3.plot = MakeSIAltAssumptionFig(here("results","altdef_ssp3.Rdata"),
                                    soln.base)

pdf(here("figures","SI_ssp3.pdf"), width=11, height=11)
ssp.3.plot$comb.plot
dev.off()

# Simply halving deforestation rates, assuming improvement in the future
halfdef.plot = MakeSIAltAssumptionFig(here("results","half_def_rate.Rdata"),
                                    soln.base)

pdf(here("figures","SI_halfdef.pdf"), width=11, height=11)
halfdef.plot$comb.plot
dev.off()

# HalvingSSP-3-based forest loss
ssp.3.halfdef.plot = MakeSIAltAssumptionFig(here("results","altdef_ssp3_halved.Rdata"),
                                    soln.base)

pdf(here("figures","SI_ssp3_halfdef.pdf"), width=11, height=11)
ssp.3.halfdef.plot$comb.plot
dev.off()

# Alternate planning horizons
T100.plot = MakeSIAltAssumptionFig(here("results","althorizon_100.Rdata"),
                                   soln.base,
                                   main.stat.col = "exp.prot.spp.by.ecoregion.new.50",
                                   calc.traj=T)

pdf(here("figures","SI_T100.pdf"), width=11, height=11)
T100.plot$comb.plot
dev.off()

T200.plot = MakeSIAltAssumptionFig(here("results","althorizon_200.Rdata"),
                                   soln.base,
                                   main.stat.col = "exp.prot.spp.by.ecoregion.new.50",
                                   calc.traj=T)

pdf(here("figures","SI_T200.pdf"), width=11, height=11)
T200.plot$comb.plot
dev.off()

##################################################
# Look at total spend as an alternate outcome
##################################################
soln.base$soln.summary[,cons.end.1p:=cons.end.period==1]
spend.map = MakeMapShaded(soln.base$soln.summary, 
                          stat.size="end.spend", 
                          stat.shade = "cons.end.period",
                          stat.linecol= "cons.end.1p",
                          shade.label = "Year last intact\n forest protected",
                          png.path=NULL, pdf.path=NULL)

pdf(here("figures","SI_spend.pdf"), width=10, height=4)
spend.map
dev.off()



# Visualize effects of some of our key assumptions on conservation:
# - Complementarity vs. endemism
# - Restoration vs. not
# - Endogenous vs. exogenous land costs

# No species complementarity
nocomp.plot = MakeSIAltAssumptionFig(here("results","nocomplementarity.Rdata"),
                                         soln.base)

pdf(here("figures","SI_nocomplementarity.pdf"), width=11, height=11)
nocomp.plot$comb.plot
dev.off()

# No restoration
norestoration.plot = MakeSIAltAssumptionFig(here("results","norestoration.Rdata"),
                                            soln.base)

pdf(here("figures","SI_norestoration.pdf"), width=11, height=11)
norestoration.plot$comb.plot
dev.off()

# No rising land costs (Fixed marginal cost of land)
norisingcosts.plot = MakeSIAltAssumptionFig(here("results","norisingcosts.Rdata"),
                                     soln.base)

pdf(here("figures","SI_norisingcosts.pdf"), width=11, height=11)
norisingcosts.plot$comb.plot
dev.off()



# Compare to fully simplified model (none of the three key features)
noextensions.plot = MakeSIAltAssumptionFig(here("results","noextensions.Rdata"),
                                            soln.base)

pdf(here("figures","SI_noextensions.pdf"), width=11, height=11)
noextensions.plot$comb.plot
dev.off()


##########################
# Tables 
##########################

# S2: correlations in optimal solution across different 
# assumptions.

# Comparison of species protection vectors in base model 
# vs. under alternate models
CompareSolutionsToBase = function(soln.base.cordata, 
                                  soln.list, 
                                  period=50) {
  all.data = soln.base.cordata$soln.traj$annot.budget.dt.long[year==period,.(eco_code, base=exp.newly.prot.spp.by.ecoregion)]
  for(i in 1:length(soln.list)) {
    to.merge = soln.list[[i]]$soln.alt$soln.traj$annot.budget.dt.long[year==period,.(eco_code, exp.newly.prot.spp.by.ecoregion)]
    setnames(to.merge, old=2, new=names(soln.list)[i])
    all.data=merge(all.data, to.merge, by="eco_code")
  }
  all.data = all.data[complete.cases(all.data),]
  return(list("comp.dt"=all.data, cors=cor(all.data[,-1])))
}

soln.base.cordata =  CombineSolutionData(here("results","main.Rdata"), 
                                         merged.data,
                                         probabilistic = T,
                                         calc.traj = T)

first10.comp = CompareSolutionsToBase(soln.base.cordata, 
                                      list("Hansen Net Loss" = hansen.net.plot,
                                           "Copernicus Loss" = copernicus.plot,
                                           "SSP Projections" = ssp.3.plot, 
                                           "Alt. Cost" = altcost.plot,
                                           "z=0.1" = altz_0.1plot,
                                           "z=0.3" = altz_0.3plot,
                                           "z from Storch et al." = altz_storchplot,
                                           "z from Kier et al." = altz_kierplot,
                                           "100 year horizon" = T100.plot,
                                           "200 year horizon" = T200.plot
                                      ),
                                      period=10)

first10.cormat = round(first10.comp$cors, 3)
first10.cormat[upper.tri(first10.cormat, diag=F)]=""
fwrite(first10.cormat, here("results","altassump_first10_cormat.csv"))

all50.comp = CompareSolutionsToBase(soln.base.cordata, 
                                    list("Hansen Net Loss" = hansen.net.plot,
                                         "Copernicus Loss" = copernicus.plot,
                                         "SSP Projections" = ssp.3.plot, 
                                         "Alt. Cost" = altcost.plot,
                                         "z=0.1" = altz_0.1plot,
                                         "z=0.3" = altz_0.3plot,
                                         "z from Storch et al." = altz_storchplot,
                                         "z from Kier et al." = altz_kierplot,
                                         "100 year horizon" = T100.plot,
                                         "200 year horizon" = T200.plot
                                    ),
                                   period=50)

all50.cormat = round(all50.comp$cors, 3)
all50.cormat[upper.tri(all50.cormat, diag=F)]=""
fwrite(all50.cormat, here("results","altassump_all50_cormat.csv"))




# S3 : detailed results
soln.base = CombineSolutionData(here("results","main.Rdata"), 
                                merged.data,
                                probabilistic = T)
all.spp.data = soln.base$soln.summary[order(-exp.prot.spp.by.ecoregion.new),
                                               .("Ecoregion" = eco_name, 
                                                 "# species" = plant_spcs,
                                                 "Initial cost ($/ha)" = GDPperm2_partial*10000, #original costs are per m^2
                                                 "SAR alpha" = ae,
                                                 "Initial forest (kha)"=(forest_area_2018)/1e3/1e4,
                                                 "Initial prot. forest (kha)"=(protected_forest_area_2018)/1e3/1e4,
                                                 "Deforest. rate (%)" = dozer_rate_2000_2018/unprotected_forest_2018*100,
                                                 "Conservation end period" = cons.end.period, # when conservation ends, even if restoration continues
                                                 "Overall end period" = end.period,
                                                 "Total spend (M$)" = first50.spend/1e6,
                                                 "Expected new spp. conserved" = exp.prot.spp.by.ecoregion.new
                                                 #"Expected new spp. uniquely conserved" = exp.add.prot.spp.by.ecoregion.new
                                               )
]


fwrite(all.spp.data, file=here("tables","allspp_allext.csv"))

# round for display version
all.spp.data.rounded = soln.base$soln.summary[order(-new.species.conserved),
                                                       .("Ecoregion" = eco_name, 
                                                         "# species" = plant_spcs,
                                                         "Initial cost ($/ha)" = round(GDPperm2_partial*10000, digits = 2), #original costs are per m^2
                                                         "SAR alpha" = round(ae, digits=2),
                                                         "Initial forest (kha)"=round((forest_area_2018)/1e3/1e4, digits=1),
                                                         "Initial prot. forest (kha)"=round((protected_forest_area_2018)/1e3/1e4, digits=1),
                                                         "Deforest. rate (%)" = round(dozer_rate_2000_2018/unprotected_forest_2018*100, digits=2),
                                                         "Conservation end period" = cons.end.period, # when conservation ends, even if restoration continues
                                                         "Overall end period" = end.period,
                                                         "Total spend (M$)" = round(first50.spend/1e6, digits=1),
                                                         "Expected new spp. conserved" = round(exp.prot.spp.by.ecoregion.new, digits=0)#,
                                                         #"Expected new spp. uniquely conserved" = round(exp.add.prot.spp.by.ecoregion.new, digits=0)
                                                         
                                                       )
]
fwrite(all.spp.data.rounded, file=here("tables","allspp_allext_rounded.csv"))

# Summarize same columns across ecoregions (mean + sd)
# for a separate panel in the same table:
sapply(colnames(all.spp.data.rounded)[2:ncol(all.spp.data.rounded)], FUN=function(cname) {
  cdata = all.spp.data.rounded[,get(cname)]
  c(mean(cdata, na.rm=T), sd(cdata, na.rm=T))
})
