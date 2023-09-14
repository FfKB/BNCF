# B.Network

require(bnlearn)
require(Rgraphviz)

# Set Works Spaces
#setwd("")

# load the data as dataframe
# This Dataset is a integration of metadata with all the clinical variable and the discretized gene expression.
# For the discretization A = <25% quantile, B = between 25% and 50%, C = between 50% and 75% and E = >75% quantile.
# The metadata Legend:
# For CF, A = CF, B = Non-CF
# Genotype, A  = Presence of the Deletion and L = absence.
# Virus and Bacteria, Virus, A = Presence, C = absence, Bacteria, B = Presence, C = absence.
# Modulator, D = No Treatment, C = Using modulators

Total <- read.delim("Datasets/Discretized_Fulldata.txt", header = T, sep=" ")
Total <- as.data.frame(unclass(Total),stringsAsFactors=TRUE) # The data must be a factor
Total <- Total[,complete.cases(Total)] # Remove Missing data, Sanity Check.

# Peforms the  bayesian network bootstrap

arcs = boot.strength(Total[,c(1:400)],  # Use only the 40 Columns in the dataset... Full data is extreme slow
                     algorithm = "hc", # gs, iamb, fast.iamb, inter.iamb, hc, tabu, mmhc, rsmax, mmpc, si.hiton.pc, chow.liu, aracne, naive.bayes and tree.bayes
                     R=100) # Number of replacement 100. 

# make the consensus network, choose your threshold wisdomly
abn = averaged.network(arcs, threshold = 0) # You can set the strenght threshold

# merge the edges strengts from bootstrap in arcs variables
mergeboot = as.data.frame(merge(abn$arcs, arcs, sort= F, by=c("from", "to")))

# creates a graphnell directed aciclic graph network
gnellnet =  graphviz.plot(abn, shape = "ellipse")

# creates a variable with edge strengh (bootstrap values)
ew <- as.character(mergeboot$strength*100)

# removes non bootstraped edges from graphnell network "caga o pipe se tirar"
ew <- ew[setdiff(seq(along = ew), removedEdges(gnellnet))]

# Associates the edge names with bootstrap values
names(ew) <- edgeNames(gnellnet)

# creates a edge attribute variable
eAttrs <- list()

# fulfill the attribute variable with bootstrap values
eAttrs$label <- ew

# generates an attribute variable for graph, as ellipse
attrs <- list(node = list(shape = "ellipse", fixedsize = F), edge = list(arrowsize = "0"))

# generates the plot with bootstrap values
plot(gnellnet, edgeAttrs = eAttrs, attrs = attrs)

# improved plot with bootstrap values. the nodelabels are best suited 
plot(gnellnet, edgeAttrs = eAttrs, attrs = list(node=list(label="foo", shape = "ellipse", fixedsize = F),
                                                edge=list(color="black", arrowsize = 0),
                                                graph=list(rankdir="UD")))

# Remove the important conections 
CF_arcs <- rbind(arcs[arcs$from == "CF",],
                 arcs[arcs$from == "Treatment",],
                 arcs[arcs$from == "Bacteria",],
                 arcs[arcs$from == "Virus",])
                 
CF_arcs <- CF_arcs[!duplicated(CF_arcs),]

CF_arcs <- CF_arcs[CF_arcs$strength >0,]

library(qgraph)
library(igraph)

Netdata = graph_from_data_frame(CF_arcs, 
                                directed = F)
#Plotar Network de uma tabela From To
plot(Netdata, 
     vertex.label.color="black",
     vertex.size=6,
     vertex.label.cex =0.6,
     edge.arrow.size=1, 
     vertex.label.dist=0.8,
     layout=layout_with_dh,
     rescale=T)
