# This script generates and analyse site- & season-specific food webs based on RiverDNA sp. occurrence & (functional group-based) metaweb data.

# (1) Analyse food webs as they are, put NA to those metrics cannot be derived with isolated nodes.
# (2) Remove isolated nodes, if they present.
# (3) Add the minimum number of sp to link them to the web.
# for (2) and (3) if a consumer (isolated or not) has no resource, add 1sp of its resource to make its trophic level make sense.


# Packages & functions
#####
library(bipartite) # for networklevel
library(igraph) # for modularity.m
library(vegan)  # for nestednodf

# Calculate unweighted TL (prey-averaged, see Williams & Martinez 2004) and Omnivory index & coherence index based on this TL.
Unw.TL.Omn.q = function(DietMatrix, Spe.Richness. = Spe.Richness, basal. = basal){
  TL1 = -t(DietMatrix)/apply(DietMatrix, 2, sum)
  TL1[basal.,] = 0
  diag(TL1) = 1
  TL2 = c(rep(1, Spe.Richness.))
  if(det(TL1) < 1e-14) {TL = NA; Omn = NA; q = NA # Avoid error from computational singular matrix. 
  } else {
    TL = solve(TL1, TL2)
    
    Omn = rep(NA, Spe.Richness.)
    for(s in 1:Spe.Richness.){
      if(s %in% basal.) next
      if(apply(DietMatrix, 2, sum)[s] == 1){ Omn[s] = 0
      }else{
        Omn[s] = sd(TL[which(DietMatrix[, s]==1)])
      }
    }
    
    Distance1 = matrix(TL, Spe.Richness., Spe.Richness., byrow = T)
    Distance2 = t(Distance1)
    Distance = DietMatrix * (Distance2 - Distance1)
    q = sd(Distance[which(Distance != 0)])
  } # end of else.
  
  return(list(TL = TL, Omn = Omn, q = q))
}

# Calculate (unipartite networks') robustness based on Dunne & Williams 2009 paper --- here, use random primary extinction sequence, and quantify R_50.
### As these food webs are mostly small and oligo-resource-supported, letting the basal nodes also vulnerable to primary removal may be a biased move.
### Thus, by default, basal species are exempted from the primary removal (i.e., Basal.rm = F, unless there are no consumers in the system, where this is no longer a food web anyway...)
### But also calculate Basal.rm = T cases.
### Also note that by definition, basal nodes should not suffer secondary extinction (as no apparent competition considered)
Uni.Robust = function(DietMatrix, No.rep = 50, Basal.rm = FALSE,
                      Basal.List = Func$Genus[Func$ppt_slide_grouping %in% c("resource", "Autotroph_bac")]){
  reading = c()
  for(i in 1:No.rep){
  Temp.D = DietMatrix
  pri = 0  # counts for primary removal of a node

    repeat{  # repeat loop for replicates of node-removing sequence.
      Basal.Node = which(rownames(Temp.D) %in% Basal.List)
      if(length(Basal.Node) == 0){break} # this is possible when Basal.rm = TRUE, and if there's no basal nodes the system is doomed, thus stop.
      pri = pri + 1
      No.Node = nrow(Temp.D)  # No.Node at the beginning of each removal step.
      if(No.Node <= 0.5*nrow(DietMatrix)){break}  # if 50% of nodes have gone extincted, stop and record the pri needed.
      if(Basal.rm == FALSE) {rm.node = sample((1:No.Node)[-Basal.Node], 1)}
      if(Basal.rm == TRUE | No.Node == length(Basal.Node)) {rm.node = sample(1:No.Node, 1)}
      print(paste("rm", rm.node, sep =" "))
      Temp.D = Temp.D[-rm.node, -rm.node] # Primary removal
      Basal.Node = which(rownames(Temp.D) %in% Basal.List) # update basal nodes location after primary removal
    
      while(length(setdiff(which(apply(Temp.D, 2, sum) == 0), Basal.Node)) != 0){ # while loop for secondary extinction following each primary node removal (if there's secondary extinction gonna happen).
        Sec.Ext = setdiff(which(apply(Temp.D, 2, sum) == 0), Basal.Node)
        print(paste("sec", Sec.Ext, sep = " "))
        Temp.D = Temp.D[-Sec.Ext, -Sec.Ext]  # update Temp.D and basal nodes location in the Temp.D following secondary extinction(s).
        Basal.Node = which(rownames(Temp.D) %in% Basal.List)
        No.Node = nrow(Temp.D)
        if(is.null(No.Node)){break}  # Null is when the above Sec.Ext happens to kill all nodes.
        if(No.Node <= 0.5*nrow(DietMatrix)){break} # When 50% ext is already achieved.
        print(paste("basal", Basal.Node, sep = " "))
      } # if no more secondary extinction is happening, exit while loop, prepare to start the next primary removal
    } # end of primary removal repeat loop.
  reading = append(reading, pri/nrow(DietMatrix))
  } # end of i for-loop. 
  return(mean(reading))
}
#####

setwd("~/Desktop/") 

# Load data
River_DNA = read.csv("Rosie Aquatic Invertebrates/Data/RiverDNA_PA_20210113.csv", header = T, check.names = F)
DMat = read.csv("Rosie Aquatic Invertebrates/Data/diet_matrix_data_20210113.csv", header = T, check.names = F)
Func = read.csv("Rosie Aquatic Invertebrates/Data/Function_genus_20210107.csv", header = T, check.names = F)
Link = read.csv("Rosie Aquatic Invertebrates/Data/trophic_links_20210113.csv", header = T, check.names = F)
# sum(DMat[,-1]) == nrow(Link) # this confirms the Diet matrix & Links datasets are consistent.

All.Func = unique(Func$ppt_slide_grouping)
DMat2 = as.matrix(DMat[, -1]) # convert DMat to matrix with taxa name as rownames.
rownames(DMat2) = DMat[, 1]
# all(rownames(DMat2) == colnames(DMat2)) # this confirms that the metaweb and its subset will have taxa arranged identically in order.
# any(apply(DMat2, 1, sum) == 0 & apply(DMat2, 2, sum) == 0) # this confirms that there's no isolated node in the metaweb.
# all(colnames(River_DNA)[-c(1:4)] %in% rownames(DMat2)) # this confirms that the genus names are consistent between the metaweb and the occurrence dataset.
# for(f in unique(Func$ppt_slide_grouping)){print(f); print(length(which(Func$ppt_slide_grouping == f)))}  # this counts the no. of genus in each functional group.

# Empty dataframes & lists for storing community/network metrics of each scheme.
Data.1 = data.frame()
Data.2 = data.frame()
Data.3 = data.frame()
D.List.1 = list()
D.List.2 = list()
D.List.3 = list()


for(scheme in 1:3){  # 1, 2, or 3

for(i in 1:nrow(River_DNA)){ #nrow(River_DNA) 
  Site_no = River_DNA[i, 2]
  Season = River_DNA[i, 3]
  River = River_DNA[i, 4]
  Sp = colnames(River_DNA)[-(1:4)][which(River_DNA[i, 5:length(River_DNA)] == 1)]
  
  # First thing, derive local Diet matrix, and see if isolated nodes present.
  ### If yes, take required modification then derive the rest based on the modified community. 
  Local.D = DMat2[which(rownames(DMat2) %in% Sp), which(colnames(DMat2) %in% Sp)]
  Iso.Node =  names(which(apply(Local.D, 1, sum) == 0 & apply(Local.D, 2, sum) == 0))
  NoRes.Con = names(which(apply(Local.D, 2, sum) == 0 & !colnames(Local.D) %in% Func$Genus[Func$ppt_slide_grouping %in% c("resource", "Autotroph_bac")]))

  # For scheme 2, remove isolated nodes then reevalute the NoRes condition; add sp if neccessary.
  # For scheme 3, add sp to resolve iso & NoRes anyway.
  ### update Sp and Local.D for running through the below.
  if(scheme == 1){ Added.Sp = c() }
  if(scheme == 2 & length(Iso.Node) != 0){
    Sp = setdiff(Sp, Iso.Node)
    Added.Sp = c()
    repeat{
      Local.D = DMat2[which(rownames(DMat2) %in% Sp), which(colnames(DMat2) %in% Sp)]
      New.NoRes.Con = names(which(apply(Local.D, 2, sum) == 0 & !colnames(Local.D) %in% Func$Genus[Func$ppt_slide_grouping %in% c("resource", "Autotroph_bac")]))
      if(length(New.NoRes.Con) != 0){ # if still NoRes after Iso-removal
        Add.Sp = names(sample(which(apply(DMat2[, New.NoRes.Con, drop=FALSE], 1, sum) == max(apply(DMat2[, New.NoRes.Con, drop=FALSE], 1, sum))), 1))
        Sp = c(Sp, Add.Sp)
        Added.Sp = append(Added.Sp, Add.Sp)
      } else {break}
    } # End of repeat loop.
  } # End of scheme 2 if statement.
  if(scheme == 3 & length(Iso.Node) != 0){
    Added.Sp = c()
    repeat{
      Local.D = DMat2[which(rownames(DMat2) %in% Sp), which(colnames(DMat2) %in% Sp)]
      New.NoRes.Con = names(which(apply(Local.D, 2, sum) == 0 & !colnames(Local.D) %in% Func$Genus[Func$ppt_slide_grouping %in% c("resource", "Autotroph_bac")]))
      New.NoCon.Res = names(which(apply(Local.D, 1, sum) == 0 & colnames(Local.D) %in% Func$Genus[Func$ppt_slide_grouping %in% c("resource", "Autotroph_bac")]))
      if(length(New.NoRes.Con) != 0){ 
        Add.Sp = names(sample(which(apply(DMat2[, New.NoRes.Con, drop = FALSE], 1, sum) == max(apply(DMat2[, New.NoRes.Con, drop = FALSE], 1, sum))), 1))
        Sp = c(Sp, Add.Sp)
        Added.Sp = append(Added.Sp, Add.Sp)
      } else if(length(New.NoCon.Res) != 0){
        Add.Sp = names(sample(which(apply(DMat2[New.NoCon.Res, , drop = FALSE], 2, sum) == max(apply(DMat2[New.NoCon.Res, , drop = FALSE], 2, sum))), 1))
        Sp = c(Sp, Add.Sp)
        Added.Sp = append(Added.Sp, Add.Sp)
      } else {break}
    } # End of repeat loop.
  } # End of scheme 3 if statement.
  
  # fundamental & function-relevant properties
  Sp.Richness = length(Sp) # Should detritus be included here?
  No.resource = length(which(Sp %in% Func$Genus[Func$ppt_slide_grouping == "resource"]))
  No.Heterotroph_bac = length(which(Sp %in% Func$Genus[Func$ppt_slide_grouping == "Heterotroph_bac"]))
  No.Autotroph_bac = length(which(Sp %in% Func$Genus[Func$ppt_slide_grouping == "Autotroph_bac"]))
  No.Zooplantkon = length(which(Sp %in% Func$Genus[Func$ppt_slide_grouping == "Zooplantkon"]))
  No.sessile_filterers = length(which(Sp %in% Func$Genus[Func$ppt_slide_grouping == "sessile_filterers"]))
  No.Collector_Filterer = length(which(Sp %in% Func$Genus[Func$ppt_slide_grouping == "Collector_Filterer"]))
  No.small_preds = length(which(Sp %in% Func$Genus[Func$ppt_slide_grouping == "small_preds"]))
  No.Shredder = length(which(Sp %in% Func$Genus[Func$ppt_slide_grouping == "Shredder"]))
  No.Grazer_scraper = length(which(Sp %in% Func$Genus[Func$ppt_slide_grouping == "Grazer_scraper"]))
  No.Invert_predator = length(which(Sp %in% Func$Genus[Func$ppt_slide_grouping == "Invert_predator"]))
  No.Omnivorous_fish = length(which(Sp %in% Func$Genus[Func$ppt_slide_grouping == "Omnivorous_fish"]))
  No.Invert_eating_fish = length(which(Sp %in% Func$Genus[Func$ppt_slide_grouping == "Invert_eating_fish"]))
  No.Piscivore_fish = length(which(Sp %in% Func$Genus[Func$ppt_slide_grouping == "Piscivore_fish"]))
  No.Parasite = length(which(Sp %in% Func$Genus[Func$ppt_slide_grouping == "Parasite"]))
  ### Summary: basal nodes vs consumer nodes.
  No.Basal = No.resource + No.Autotroph_bac # Heterotroph bacs are consumers now.
  No.Consumer = Sp.Richness - No.Basal
  ### Below according to cheddar categories
  No.Detritus = No.resource
  No.Decomposer = No.Heterotroph_bac
  No.Producer = No.Autotroph_bac
  No.Detritivore = No.Zooplantkon + No.Shredder + No.sessile_filterers
  No.Herbivore = No.Collector_Filterer + No.Grazer_scraper
  No.Predator = No.small_preds + No.Invert_predator + No.Parasite + No.Omnivorous_fish + No.Invert_eating_fish + No.Piscivore_fish
  No.Vert = No.Omnivorous_fish + No.Invert_eating_fish + No.Piscivore_fish # vertebrates, in this dataset, fish.
  ### G: based on our functional grouping. C: based on cheddar functional categories.
  Func.G.Diversity = length(which(c(No.resource, No.Heterotroph_bac, No.Autotroph_bac, No.Zooplantkon, No.sessile_filterers,
                              No.Collector_Filterer, No.small_preds, No.Shredder, No.Grazer_scraper, No.Invert_predator,
                              No.Omnivorous_fish, No.Invert_eating_fish, No.Piscivore_fish, No.Parasite) > 0))
  Func.G.Redundancy = Sp.Richness/Func.G.Diversity # Detritus is now included. But should it?
  Func.C.Diversity = length(which(c(No.Detritus, No.Decomposer, No.Producer, No.Detritivore, No.Herbivore, No.Predator) > 0))
  Func.C.Redundancy = Sp.Richness/Func.C.Diversity # Detritus is now included. But should it?
  
  # Food-web metrics        
  No.Link = sum(Local.D)
  Link.Density = No.Link/Sp.Richness
  Connectance = No.Link/(Sp.Richness^2)
  Nic.over.Con = networklevel(Local.D, index = "niche overlap")[1] # viewed as consumers vs. resources but actually the same set of taxa.
  Nic.over.Res = networklevel(Local.D, index = "niche overlap")[2]
  Mean.Gen = mean(apply(Local.D, 2, sum))
  SD.Gen = sd(apply(Local.D, 2, sum))
  Mean.Vul = mean(apply(Local.D, 1, sum))
  SD.Vul = sd(apply(Local.D, 1, sum))  
  ### Trophic level relevant
  TL.info = Unw.TL.Omn.q(Local.D, Sp.Richness, which(apply(Local.D, 2, sum)==0))
  Max.TL = max(TL.info$TL)
  Mean.TL = mean(TL.info$TL)
  Omnivory = mean(TL.info$Omn, na.rm = T) # In terms of TL, not primary producer
  Coherence = TL.info$q # Lower q, higher trophic coherence.
  ### Topological
  Nestedness = unname(nestednodf(Local.D)$statistic)[3]
  # Modularity = modularity(cluster_walktrap(graph.adjacency(Local.D)))
  Modularity = modularity(multilevel.community(graph.adjacency(Local.D, mode = "undirected")))
  R50.Basal.excl = Uni.Robust(Local.D)
  R50.Basal.incl = Uni.Robust(Local.D, Basal.rm = T)
  
  # Fill in dataframe
  temp.data = data.frame(
    Site_no = Site_no,
    Season = Season,
    River = River,
    No.Iso.Node = length(Iso.Node),
    No.NoRes.Con = length(NoRes.Con),
    Sp.Richness = Sp.Richness,
    No.resource = No.resource,
    No.Heterotroph_bac = No.Heterotroph_bac,
    No.Autotroph_bac = No.Autotroph_bac,
    No.Zooplantkon = No.Zooplantkon,
    No.sessile_filterers = No.sessile_filterers,
    No.Collector_Filterer = No.Collector_Filterer,
    No.small_preds = No.small_preds,
    No.Shredder = No.Shredder,
    No.Grazer_scraper = No.Grazer_scraper,
    No.Invert_predator = No.Invert_predator,
    No.Omnivorous_fish = No.Omnivorous_fish,
    No.Invert_eating_fish = No.Invert_eating_fish,
    No.Piscivore_fish = No.Piscivore_fish,
    No.Parasite = No.Parasite,
    No.Basal = No.Basal,
    No.Consumer = No.Consumer,
    No.Detritus = No.Detritus,
    No.Decomposer = No.Decomposer,
    No.Producer = No.Producer,
    No.Detritivore = No.Detritivore,
    No.Herbivore = No.Herbivore,
    No.Predator = No.Predator,
    No.Vert = No.Vert,
    Func.G.Diversity = Func.G.Diversity,
    Func.G.Redundancy = Func.G.Redundancy,
    Func.C.Diversity = Func.C.Diversity,
    Func.C.Redundancy = Func.C.Redundancy,
    No.Link = No.Link,
    Link.Density = Link.Density,
    Connectance = Connectance,
    Nic.over.Con = Nic.over.Con,
    Nic.over.Res = Nic.over.Res,
    Mean.Gen = Mean.Gen,
    SD.Gen = SD.Gen,
    Mean.Vul = Mean.Vul,
    SD.Vul = SD.Vul, 
    Max.TL = Max.TL,
    Mean.TL = Mean.TL,
    Omnivory = Omnivory,
    Coherence = Coherence,
    Nestedness = Nestedness,
    Modularity = Modularity,
    R50.Basal.excl,
    R50.Basal.incl)
  if(scheme == 1){
    Data.1 = rbind(Data.1, temp.data)
    D.List.1[[i]] = list(
      Site_no = Site_no,
      Season = Season,
      River = River,
      Iso.Node = Iso.Node,  # before modification
      NoRes.Con = NoRes.Con, # before modification
      Added.Sp = Added.Sp, # modification
      Local.D = Local.D) # after modification
  }
  if(scheme == 2){
    Data.2 = rbind(Data.2, temp.data)
    D.List.2[[i]] = list(
      Site_no = Site_no,
      Season = Season,
      River = River,
      Iso.Node = Iso.Node,  # before modification
      NoRes.Con = NoRes.Con, # before modification
      Added.Sp = Added.Sp, # modification
      Local.D = Local.D) # after modification
  }
  if(scheme == 3){
    Data.3 = rbind(Data.3, temp.data)
    D.List.3[[i]] = list(
      Site_no = Site_no,
      Season = Season,
      River = River,
      Iso.Node = Iso.Node,  # before modification
      NoRes.Con = NoRes.Con, # before modification
      Added.Sp = Added.Sp, # modification
      Local.D = Local.D) # after modification
  }
} # End of i for-loop.

} # End of scheme for-loop.


# Export data
# rm(list = setdiff(ls(), c("Data.1", "Data.2", "Data.3", "D.List.1", "D.List.2", "D.List.3")))
# dir.create("Rosie Aquatic Invertebrates/Output/", showWarnings = F)
# setwd("Rosie Aquatic Invertebrates/Output/")
# write.table(Data.1, col.names=T, row.names=F, sep =",", "Data_1_20210521.csv")
# write.table(Data.2, col.names=T, row.names=F, sep =",", "Data_2_20210521.csv")
# write.table(Data.3, col.names=T, row.names=F, sep =",", "Data_3_20210521.csv")
# save.image("Data_AllSchemes_20210521.RData")
