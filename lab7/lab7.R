library(ggplot2)
library(RColorBrewer)
library(rgl)

source("/home/pietromalagoli/PCN/common.R")
mypal <- viridis::viridis_pal()(20)
set.seed(12345)

g <- watts.strogatz.game(dim=2, size = 20, p = 0.05, nei=1)
lay2D <- layout_with_fr(g, dim=2)
lay3D <- layout_with_fr(g, dim=3)
deg <- degree(g)
centr_deg <- deg
cols_deg <- vec2pal(centr_deg, mypal)
clear3d()
rglplot(g, layout=lay3D, vertex.color=cols_deg, vertex.size=deg, vertex.label=NA, edge.color="gray80", vertex.frame.color=NA, edge.width=0.5)


centr_betw <- betweenness(g, normalized = TRUE)
cols_betw <- vec2pal( log(1 + centr_betw), mypal)
clear3d()
rglplot(g, layout=lay3D, vertex.color=cols_deg, vertex.size=deg, vertex.label=NA, edge.color="gray80", vertex.frame.color=NA, edge.width=0.5)

centr_close <- closeness(g, normalized = TRUE)

centr_pr <- page_rank(g)$vector