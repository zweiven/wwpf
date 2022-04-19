# File to make main text figures 1 and 2, assuming
# all solution files with results of simulations are available
# in specified locations.

library(here)
library(sf)
library(data.table)
library(ggplot2)
library(cowplot)
library(rmapshaper)
library(rnaturalearthdata)
library(Rfast)

# load helper functions
source(here("code","SolutionVizHelpers.R"))

profile.data = fread(here("data","mammal_based_profiles.csv"))
profile.abundance = profile.data[,filt.abundances]               
profiles = as.matrix(profile.data[,-ncol(profile.data), with=F])


############
# Figure 1
############
soln.allext = CombineSolutionData(here("results","main.Rdata"), 
                                  merged.data,
                                  probabilistic = T)
# Identify areas in which conservation of intact forest
# is done in period 1
allext.plot.data = soln.allext$soln.summary
allext.plot.data[,cons.end.1p:=cons.end.period==1]
setorder(allext.plot.data, -cons.end.period)

# Make plot
fig1 = MakeMapShaded(solution.dt = allext.plot.data,
                     stat.size = "exp.prot.spp.by.ecoregion.new",
                     stat.shade = "cons.end.period",
                     stat.linecol= "cons.end.1p",
                     shade.label = "Year last intact\n forest protected",
                     size.label = "Expected species\n newly protected")
pdf(here("figures","fig1.pdf"), width=10, height=5)
fig1
dev.off()


############
# Figure 2
############
# Panel 2A: benefit-cost ratio driving solution
fig2a = ggplot(soln.allext$soln.summary[is.finite(initial.mb.mc.ratio),],
                               aes(x=initial.mc*10000,
                                   y=initial.mb*10000,
                                   fill=threat,
                                   size=round(exp.prot.spp.by.ecoregion.new))) + 
  geom_abline(slope=1, intercept=-15.9, colour="#888888", linetype="dashed") +
  geom_point(pch=21) + 
  theme_minimal() + 
  scale_x_continuous(
    trans="log",
    breaks=c(10, 100, 1000, 10000, 100000)) + 
  scale_y_continuous(
    trans="log", 
    breaks=c(0.00001,0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000)) + 
  scale_fill_gradient(name=expression(Threat~(yr^-1)),
                      low="#ffffff", high="#ff0000", trans="log",
                      breaks=c(0.00001, 0.0001, 0.001, 0.01, 0.1),
                      guide = guide_colorbar(barwidth=8)) + 
  scale_size_continuous("Expected # species protected",
                        breaks=c(0,100,500,1000),
                        range = c(1,8),
                        limits=c(0,1000)) +
  theme(axis.line = element_line(), 
        panel.grid = element_blank(),
        legend.position = c(0.17, 0.675)) + 
  xlab("Initial protection cost ($/hectare)") + 
  ylab("Initial benefit (expected additional species/hectare)") + 
  guides(size=guide_legend(order=1))

# Panel 2B: comparison with heuristics

# Load main solution
load(here("results","main.Rdata"))

# Load heuristic solutions
load(here("results","heur_bcratio.Rdata"))
heur.soln.bcratio = heur.soln

load(here("results", "heur_hotspot.Rdata"))
heur.soln.hotspot = heur.soln

load(here("results","heur_mb.Rdata"))
heur.soln.mb = heur.soln

load(here("results","heur_threat.Rdata"))
heur.soln.threat = heur.soln

load(here("results","heur_cost.Rdata"))
heur.soln.cost = heur.soln

load(here("results", "heur_percprot.Rdata"))
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

fig2b = ggplot(new.spp.data, 
               aes(x=as.integer(as.character(year)), 
                  y=tot.protected.spp.in.year, group=soln, 
                   colour=soln, 
                   linetype=soln)) + 
  geom_line(size=1) + 
  theme_minimal() + 
  xlab("Year") +
  ylab("Total expected protected species") +
  scale_x_continuous(expand=c(0,0), breaks=seq(from=0,to=50,by=5)) + 
  scale_y_continuous(expand=c(0,0)) + 
  scale_colour_viridis_d(name="Allocation rule", begin=0.1, end=0.75) + 
  guides(color=guide_legend(nrow=4, byrow=TRUE, title.position="top")) +
  scale_linetype(name="Allocation rule") +
  theme(legend.position = c(0.2, 0.85),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size=0.5),
        legend.key.width = unit(4,"line"))

# Make Fig 2 as two-panel figure
bottom.leg = theme(legend.position="bottom", legend.box="vertical")
pdf(here("figures","fig2.pdf"), width=11, height=6.25)
plot_grid(fig2a+bottom.leg + theme(plot.margin=unit(c(0.2,0,1.07,0),"cm")), 
          fig2b+bottom.leg, nrow=1, ncol=2, labels = c("A","B"),
          align = "v")
dev.off()



################
# Calculation of other numbers included in text
################
# Load main solution
load(here("results","main.Rdata"))

# Increase in expected species protected during
# the planning horizon
print(paste0("Increase in expected species protected: ",
             round(soln.trajectories$tot.exp.newly.prot.spp[50])))

# Number of regions in which $ is spent
print(paste0("# ecoregions in which $ is spent: ",
             sum(rowsums(global.soln$budget)>0)))

# ... during the first year
print(paste0("# ecoregions in which $ is spent in the first year: ",
             sum(global.soln$budget[,1]>0)))

# ... during the first 10
print(paste0("# ecoregions in which $ is spent in the first 10 years: ",
             sum(rowsums(global.soln$budget[,1:10])>0)))

# In how many of those regions where we spend in the first year
# do we protect all remaining intact forest?
# Note we check both that all available forest is conserved and the
# last period of conservation is period 1. There are two regions in which
# the last period of conservation is period 1 but we run out of budget, and so 
# deforestation removes remaining forest before conservation can resume in period 2
print(paste0("# ecoregions in which all intact forest is protected in year 1: ",
             sum(global.soln$race.end$end.areas.conserved==global.soln$land.setup$area.forested.init & global.soln$race.end$cons.end.periods==1)))

# If we instead want to know all of the areas where there's no intact 
# forest after the first period, whether or not that's due to conservation
# or protection:
print(paste0("# ecoregions in which no unprotected intact forest remains after year 1: ",
             sum(global.soln$race.end$cons.end.periods==1)))

# How much conservation occurs in these ecoregions?
yr1.end.areas = global.soln$land.setup$eco_code[which(global.soln$race.end$cons.end.periods==1)]
yr1.end.areas.stats.yr50 = soln.trajectories$annot.budget.dt.long[year==50,.(tot.new.exp.spp=sum(exp.newly.prot.spp.by.ecoregion)), by=eco_code %in% yr1.end.areas]
                                       
# what fraction of gains come from these areas (ignoring complementarity)
print(paste0("Fraction of expected increase in spp protection in priority ecoregions after 50 years:",
             yr1.end.areas.stats.yr50[eco_code==T,tot.new.exp.spp]/sum(yr1.end.areas.stats.yr50$tot.new.exp.spp)))


# Increase in expected species protected that occurs
# during the first 10 years
print(paste0("Increase in expected species protected in first 10 years: ",
             round(soln.trajectories$tot.exp.newly.prot.spp[10])))
# ...as a fraction of the protection that occurs during the 
# entire horizon
print(paste0("Increase in expected species protected in first 10 years as fraction of gains over 50 year horizon: ",
             round(soln.trajectories$tot.exp.newly.prot.spp[10]/soln.trajectories$tot.exp.newly.prot.spp[50],3)))




# Comparing expected increases in species protection
# under the optimal solution with heuristics:
new.spp.data[year==50,.(soln,frac.of.optimal=tot.protected.spp.in.year/tot.protected.spp.in.year[soln=="Optimal Solution"])]

# Compare heruistics with optimal solution in terms of cost savings.
# we do this using the simulations using varying budgets; we can find
# the smallest budget which lets us achieve (under the optimal solution) 
# at least the same increase in expected species protection as the heuristics.
ComputeSppConserved = function(budget) {
  load(here("results",paste0("altbudget_",budget,".Rdata")))
  budget.traj = ComputeTrajectories(global.soln$budget, ecoregion.data, profiles=profiles,profile.abundance = profile.abundance)
  budget.traj$tot.exp.newly.prot.spp[50]
}
budgets = c(1E7,5E7,1E8,2.5E8,5E8,7.5E8,1E9,1E10,2.5E10)
budget.text = c("1e+07","5e+07","1e+08","2.5e+08","5e+08","7.5e+08","1e+09","1e+10", "2.5e+10")
budget.res = data.table(budget=budgets,
                        new.exp.spp.conserved = sapply(budget.text, FUN=ComputeSppConserved))

# How much better could we do if we had all budget up front?
load(here("results","upfrontbudget.Rdata"))
upfront.traj = ComputeTrajectories(global.soln$budget, ecoregion.data, profiles=profiles,profile.abundance = profile.abundance)
print(paste0("Species protection improvement with up front budget: ",
             "total: ", upfront.traj$tot.exp.newly.prot.spp,
             ", delta: ", upfront.traj$tot.exp.newly.prot.spp-soln.trajectories$tot.exp.newly.prot.spp[50],
             ", % imp: ", (upfront.traj$tot.exp.newly.prot.spp-soln.trajectories$tot.exp.newly.prot.spp[50])/soln.trajectories$tot.exp.newly.prot.spp[50]))
