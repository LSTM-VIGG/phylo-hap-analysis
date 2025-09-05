### Define input ####

# input files
library(ape)
library(phytools)
library(stringr)
library(data.table)
library(dplyr)

# edit these paths accordingly 
prefix    = "results/phylo/"
phy_list  = c("results/phylo/coeae1f_focal.fasta.treefile", "results/phylo/coeae1f_upstream.fasta.treefile", "results/phylo/coeae1f_downstream.fasta.treefile")
phy_name  = c("Focal haplotype","Upstream","Downstream")

phi = "results/phylo/coeae1f_focal.fasta.treefile"
meta_path = "results/phylo/coeae1f_focal.metadata.tsv"
cou = 0

# tree to distance matrix
cou    = cou+1
phy    = read.tree(phi)
dis    = cophenetic.phylo(phy)
meta = fread(meta_path) %>% as.data.frame()
phy_p = midpoint.root(phy)
meta = meta %>% arrange(factor(hap, levels = phy_p$tip.label))

phy_p$tip_label_sep  = gsub("_"," ",phy_p$tip.label)
phy_p$tip_species    = meta$aim_species
phy_p$tip_country = meta$country
phy_p$tip_cluster   = ""
phy_p$tip_karyotype  = ""
#phy_p$tip_color     = brewe
phy_p$tip.label      = rep("·", length(phy$tip.label))
phy_p$edge.length[phy_p$edge.length >  0.005] = 0.005
phy_p$edge.length[phy_p$edge.length == 0]    = 5e-5


plot_phylo = function(meta, var){
  if (var == 'karyotype'){
    meta = meta %>% mutate("tip_color" = case_when(karyotype == '2l+a' ~ 'bisque2',
                                                   karyotype == '2la' ~ 'dodgerblue', TRUE ~ 'grey'))
  } else if (var == 'Sweep IDs'){
    meta$tip_color = rainbow(length(unique(meta[,var])))[as.factor(meta[, var])]
    meta = meta %>% mutate(tip_color = case_when(`Sweep IDs` == 'wt' ~ 'grey',
                                                 TRUE ~ tip_color))
  } else if (var == 'aim_species'){
    meta = meta %>% mutate("tip_color" = case_when(aim_species == 'gambiae' ~ 'indianred',
                                                   aim_species == 'coluzzii' ~ 'dodgerblue',
                                                   aim_species == 'arabiensis' ~ 'aquamarine',
                                                   aim_species == 'mela' ~ 'cornsilk',
                                                   aim_species == 'mera' ~ 'cornsilk4',
                                                   aim_species == 'quad' ~ 'darkolivegreen',
                                                   TRUE ~ 'grey'))
  }
  print(plot.phylo(phy_p, type="unr",
             use.edge.length=T, show.tip.label=T, show.node.label=F,
             tip.color = meta$tip_color, lab4ut = "axial",
             edge.color = "slategray3",
             font = 1, edge.width = 2, node.depth=1, cex=5,
             main=phy_name[cou], no.margin = TRUE))
}

meta = meta %>% mutate("tip_color" = case_when(aim_species == 'gambiae' ~ 'indianred',
                                            aim_species == 'coluzzii' ~ 'dodgerblue',
                                            aim_species == 'arabiensis' ~ 'aquamarine',
                                            aim_species == 'mela' ~ 'cornsilk',
                                            aim_species == 'mera' ~ 'cornsilk4',
                                            aim_species == 'quad' ~ 'darkolivegreen',
                                        TRUE ~ 'grey'))

var = 'karyotype' # ' Sweep IDs'

meta = meta %>% mutate("tip_color" = case_when(karyotype == '2l+a' ~ 'bisque2',
                                               karyotype == '2la' ~ 'dodgerblue', TRUE ~ 'grey'))

meta$tip_color = rainbow(length(unique(meta[,var])))[as.factor(meta[, var])]
meta = meta %>% mutate(tip_color = case_when(`Sweep IDs` == 'wt' ~ 'grey',
                                              TRUE ~ tip_color))

plot.phylo(phy_p, type="unr",
           use.edge.length=T, show.tip.label=T, show.node.label=F,
           tip.color = meta$tip_color, lab4ut = "axial",
           edge.color = "slategray3",
           font = 1, edge.width = 2, node.depth=1, cex=5,
           main=phy_name[cou], no.margin = TRUE)




phy_p$tip_colors

#  plot phylogenies
pdf(file=paste(prefix,"dists_phylo.pdf",sep="."),height=6,width=6)
cou      = 0
for (phi in phy_list) {
  
  # tree to distance matrix
  cou    = cou+1
  phy    = read.tree(phi)
  dis    = cophenetic.phylo(phy)
  
  phy_p = midpoint.root(phy)
  phy_p$tip.label_sep  = gsub("_"," ",phy_p$tip.label)
  phy_p$tip.species    = gsub("^[A-Z][A-Z]","",word(phy_p$tip.label_sep,2))
  phy_p$tip.species[phy_p$tip.species == ""] = "gamcol"
  phy_p$tip.population = word(phy_p$tip.label_sep,2)
  phy_p$tip.genotype   = word(phy_p$tip.label_sep,3)
  phy_p$tip.karyotype  = as.factor(word(phy_p$tip.label_sep,4))
  phy_p$tip.ktgt       = as.factor(paste(phy_p$tip.karyotype,paste(phy_p$tip.genotype,phy_p$tip.species,sep="_"),sep="_"))
  phy_p$tip.colors     = c("orangered3","orangered3","orangered3","deeppink3","deeppink4",         # kt0_gt0_col, kt0_gt0_gam, kt0_gt0_gamcol, kt0_gt0_mel, kt0_gt0_qua
                           "orange","orange","deeppink",                                           # kt0_gt1_col, kt0_gt1_gam, kt0_gt2_qua
                           "springgreen4","dodgerblue3","dodgerblue3","dodgerblue3","dodgerblue4", # kt2_gt0_ara, kt2_gt0_col, kt2_gt0_gam, kt2_gt0_gamcol, kt2_gt0_mer
                           "purple","purple",                                                      # kt2_gt1_col, kt2_gt1_gam, 
                           "springgreen3","turquoise2"                                             # kt2_gt2_ara, kt2_gt2_col
                           )[phy_p$tip.ktgt]
  phy_p$tip.label      = rep("·", length(phy$tip.label))
  phy_p$edge.length[phy_p$edge.length >  0.005] = 0.005
  phy_p$edge.length[phy_p$edge.length == 0]    = 5e-5
  
  plot.phylo(phy_p, type="unr",
             use.edge.length=T, show.tip.label=T, show.node.label=F,
             tip.color = phy_p$tip.colors, lab4ut = "axial",
             edge.color = "slategray3",
             font = 0.5, edge.width = 0.5, node.depth=1,
             main=phy_name[cou])
  
  
  
}

dev.off()


pdf(file=paste(prefix,"dists_phylo_llarg.pdf",sep="."),height=150,width=30)
cou      = 0
for (phi in phy_list) {
  
  # tree to distance matrix
  cou    = cou+1
  phy    = read.tree(phi)
  dis    = cophenetic.phylo(phy)
  
  phy_p = midpoint.root(phy)
  phy_p$tip.label_sep  = gsub("_"," ",phy_p$tip.label)
  phy_p$tip.species    = gsub("^[A-Z][A-Z]","",word(phy_p$tip.label_sep,2))
  phy_p$tip.species[phy_p$tip.species == ""] = "gamcol"
  phy_p$tip.population = word(phy_p$tip.label_sep,2)
  phy_p$tip.genotype   = word(phy_p$tip.label_sep,3)
  phy_p$tip.karyotype  = as.factor(word(phy_p$tip.label_sep,4))
  phy_p$tip.ktgt       = as.factor(paste(phy_p$tip.karyotype,paste(phy_p$tip.genotype,phy_p$tip.species,sep="_"),sep="_"))
  phy_p$tip.colors     = c("orangered3","orangered3","orangered3","deeppink3","deeppink4",         # kt0_gt0_col, kt0_gt0_gam, kt0_gt0_gamcol, kt0_gt0_mel, kt0_gt0_qua
                           "orange","orange","deeppink",                                           # kt0_gt1_col, kt0_gt1_gam, kt0_gt2_qua
                           "springgreen4","dodgerblue3","dodgerblue3","dodgerblue3","dodgerblue4", # kt2_gt0_ara, kt2_gt0_col, kt2_gt0_gam, kt2_gt0_gamcol, kt2_gt0_mer
                           "purple","purple",                                                      # kt2_gt1_col, kt2_gt1_gam, 
                           "springgreen3","turquoise2"                                             # kt2_gt2_ara, kt2_gt2_col
  )[phy_p$tip.ktgt]
  
  plot.phylo(phy_p, type="phy",
             use.edge.length=T, show.tip.label=T, show.node.label=T,
             tip.color = phy_p$tip.colors, lab4ut = "axial",
             edge.color = "slategray3",
             font = 0.5, node.depth=1, underscore=F,
             main=phy_name[cou])
  
  
  
}

dev.off()


stop("ara")

