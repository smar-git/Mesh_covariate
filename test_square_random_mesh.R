#' Here we simulate three PP using  a LGCP where the log intensity
#' depends on only an intercept and a covariate
#' 
#' $$
#' \beta_0 + beta_1\ x(s)
#' $$
#' 
#' We then estimate the model using different meshes, all meshes have 
#' similar parameters but different locations for the mesh nodes.


#' build the mesh with points instead of shape, use a lattice input to build the mesh
#' 
#' check correlation range of the covariate and the size of the mesh, ratio between correlation of covariate
#' and mesh length
#' 
#' in slide 24 have both size and color depend on weights
#' scale_size_area
#' 
#' 
#' slide 31 flip int scheme and spde scheme
#' 
#' dual mesh (voronii cells) and integration cells
#' 
#' paper guy from oslo taht looks at the mesh


#+ echo = FALSE
knitr::opts_chunk$set(collapse = TRUE, echo = FALSE,
                      message = FALSE, warning = FALSE,
                      fig.width=7, fig.height=7)



library(tidyverse)
library(inlabru)
library(patchwork)
library(INLA)
library(sf)
library(scico)
library(tidyterra)
library(scico)
library(terra)
rm(list = ls())
source("book_mesh_dual.R")
set.seed(1234)

plot_save_dir  = "presentation/plots/random_mesh/"
theme_maps = list(theme(axis.title.x=element_blank(),
                        axis.text.x=element_blank(),
                        axis.ticks.x=element_blank(),
                        axis.title.y=element_blank(),
                        axis.text.y=element_blank(),
                        axis.ticks.y=element_blank()),
                  scale_fill_scico(na.value = "transparent"),
                  scale_color_scico(na.value = "transparent"))



# define domain -----------------------------------------------------------



max_dom  = 500
res = 1

max_dom  = 500
res = 1

box = c(xmin = 0 , 
        ymin = 0,
        xmax = max_dom,
        ymax = max_dom)
poly <- data.frame(x = box[c(1,3)], y = box[c(2,4)]) %>% 
  st_as_sf(coords = c("x", "y")) %>% 
  st_bbox() %>% 
  st_as_sfc()

# polys = list()
# polys[[1]] = poly
# polys[[2]]  = data.frame(x = box[c(1,3)]-5, y = box[c(2,4)]+10) %>% 
#   st_as_sf(coords = c("x", "y")) %>% 
#   st_bbox() %>% 
#   st_as_sfc()
# polys[[3]]  = data.frame(x = box[c(1,3)]+10, y = box[c(2,4)]+10) %>% 
#   st_as_sf(coords = c("x", "y")) %>% 
#   st_bbox() %>% 
#   st_as_sfc()
# polys[[4]]  = data.frame(x = box[c(1,3)]+20, y = box[c(2,4)]-20) %>% 
#   st_as_sf(coords = c("x", "y")) %>% 
#   st_bbox() %>% 
#   st_as_sfc()
# polys[[5]]  = data.frame(x = box[c(1,3)]+7, y = box[c(2,4)]+10) %>% 
#   st_as_sf(coords = c("x", "y")) %>% 
#   st_bbox() %>% 
#   st_as_sfc()
# polys[[6]]  = data.frame(x = box[c(1,3)]-20, y = box[c(2,4)]+0) %>% 
#   st_as_sf(coords = c("x", "y")) %>% 
#   st_bbox() %>% 
#   st_as_sfc()

cov1 <- rast(nrows = max_dom/res, ncols=max_dom/res, xmin=0, xmax=max_dom, ymin = 0, ymax = max_dom)
cov2 <- rast(nrows = max_dom/res, ncols=max_dom/res, xmin=0, xmax=max_dom, ymin = 0, ymax = max_dom)
cov3 <- rast(nrows = max_dom/res, ncols=max_dom/res, xmin=0, xmax=max_dom, ymin = 0, ymax = max_dom)



n_meshes = 6

xy0 = expand.grid(seq(-100,600,50),
                 seq(-100,600,50))


## create meshes and int points -----------------------------------------------------------
#meshes = list()

#for(i in 1:n_meshes)
#    meshes[[i]] =  fm_mesh_2d_inla(boundary = polys[[i]],
#                               max.edge = c(75,80))


meshes = list()
cc = cbind(c(0,20, 10, 0,-20,-20),
           c(0,20, -10, 30,0,-10))
for(i in 1:n_meshes)
{
  xy = xy0
  xy[,1] = xy[,1]+cc[i,1]
  xy[,2] = xy[,2]+cc[i,2]
  
    meshes[[i]] =  fm_mesh_2d_inla(loc = xy)
}

ggplot() + gg(meshes[[1]]) +  geom_sf(data = poly, alpha = 0, color = "red") +
  coord_sf(xlim = c(-150,650), 
           ylim = c(-150,650)
    ) + 
  ggtitle("Mesh 1") +
  xlab("") + ylab("") + theme_maps +
ggplot() + gg(meshes[[2]]) +  geom_sf(data = poly, alpha = 0, 
                                      color = "red", size = 2)+
  coord_sf(xlim = c(-150,650), ylim = c(-150,650))+ ggtitle("Mesh 2") +
  xlab("") + ylab("") +theme_maps +
ggplot() + gg(meshes[[3]]) +  geom_sf(data = poly, alpha = 0, color = "red")+
  coord_sf(xlim = c(-150,650), ylim = c(-150,650))+ ggtitle("Mesh 3") +
  xlab("") + ylab("") +theme_maps +
ggplot() + gg(meshes[[4]]) +  geom_sf(data = poly, alpha = 0, color = "red")+
  coord_sf(xlim = c(-150,650), ylim = c(-150,650))+ ggtitle("Mesh 4") +
  xlab("") + ylab("") +theme_maps +
ggplot() + gg(meshes[[5]]) +  geom_sf(data = poly, alpha = 0, color = "red")+
  coord_sf(xlim = c(-150,650), ylim = c(-150,650))+ ggtitle("Mesh 5") +
  xlab("") + ylab("") +theme_maps +
ggplot() + gg(meshes[[6]]) +  geom_sf(data = poly, alpha = 0, color = "red")+
  coord_sf(xlim = c(-150,650), ylim = c(-150,650))+ ggtitle("Mesh 6") +
  xlab("") + ylab("")+  theme_maps 

 ggsave(paste(plot_save_dir,"meshes.png", sep = ""))


sapply(meshes, function(x)x$n)


int_points = list()
for(i in 1:n_meshes)
{
  print(paste("compute int. point  for mesh", i))
  int_points[[i]] = fm_int(meshes[[i]], samplers = poly)
}





# simulate a covariate ----------------------------------------------------

sim_mesh <- fm_mesh_2d(loc.domain = (poly), 
                         max.edge = c(10, 40))
sim_matern <- inla.spde2.pcmatern(sim_mesh,
                                  prior.sigma = c(1, 0.5),
                                  prior.range = c(1, 0.5))
A = inla.spde.make.A(mesh = sim_mesh, loc  = crds(cov1))


range1 = 10
sigma1 = 1
Q <- inla.spde2.precision(spde = sim_matern,
                          theta = c(log(range1),log(sigma1)))
samp1 <- inla.qsample(n = 2,
                     Q = Q,
                     mu = rep(0,dim(Q)[1]))
val =  as.vector((A%*%samp1[,1]))
#val =   ifelse(val<1.5,0,val)
values(cov1) = val #- mean(val[val>0])
names(cov1) =  "val"

#ggplot() + geom_spatraster(data = cov1) + scale_fill_scico( direction = -1) +
#  gg(meshes[[1]])


range2 = 50
sigma2 = 1
Q <- inla.spde2.precision(spde = sim_matern,
                          theta = c(log(range2),log(sigma2)))
samp2 <- inla.qsample(n = 2,
                      Q = Q,
                      mu = rep(0,dim(Q)[1]))
val =  as.vector((A%*%samp2[,1]))
#val = ifelse(val<1,0,val)
values(cov2) = val# - mean(val[val>0])
names(cov2) =  "val"
#ggplot() + geom_spatraster(data = cov2) + scale_fill_scico( direction = -1) +
#  gg(meshes[[1]]) + coord_equal()


range3 = 600
sigma3 = 1
Q <- inla.spde2.precision(spde = sim_matern,
                          theta = c(log(range3),log(sigma3)))
samp3 <- inla.qsample(n = 2,
                     Q = Q,
                     mu = rep(0,dim(Q)[1]))
val =  as.vector((A%*%samp3[,1]))
values(cov3) = val - mean(val)
names(cov3) =  "val"



ggplot() + geom_spatraster(data = cov1) + geom_sf(data = poly,
                                                  alpha = 0) +
  theme_maps + theme(legend.position  = "none") +
  ggtitle("Covariate 1") + 
  ggplot() + geom_spatraster(data = cov2) + geom_sf(data = poly,
                                                    alpha = 0) +
  theme_maps + theme(legend.position   = "none") +
  ggtitle("Covariate 2") + 
  ggplot() + geom_spatraster(data = cov3) + geom_sf(data = poly,
                                                    alpha = 0) +
  theme_maps + theme(legend.position   = "none") +
  ggtitle("Covariate 3") + plot_layout(ncol = 2)
ggsave(paste(plot_save_dir,"covariates.png",sep =""))
  


#ggplot() + geom_spatraster(data = cov3) + scale_fill_scico( direction = -1) +
#  gg(meshes[[1]])

# covariates at mesh nodes ------------------------------------------------



mesh_nodes = c()
for(i in 1:n_meshes)
{
  mesh_nodes  = rbind(mesh_nodes, 
                      data.frame(x = meshes[[i]]$loc[,1],
                                 y = meshes[[i]]$loc[,2],
                                 mesh = i))
 
}

mesh_nodes = st_as_sf(mesh_nodes, coords = c("x","y"))
mesh_nodes$cov1 = extract(cov1, mesh_nodes)$val
mesh_nodes$cov2 = extract(cov2, mesh_nodes)$val
mesh_nodes$cov3 = extract(cov3, mesh_nodes)$val


mesh_nodes$cov1 =  bru_fill_missing(cov1, mesh_nodes, mesh_nodes$cov1)
mesh_nodes$cov2 =  bru_fill_missing(cov2, mesh_nodes, mesh_nodes$cov2)
mesh_nodes$cov3 =  bru_fill_missing(cov3, mesh_nodes, mesh_nodes$cov3)

head(mesh_nodes)
dual = list()
int = list()
weights = c()
for(i in 1:n_meshes)
{
  
  dual[[i]] = st_as_sf(book.mesh.dual(meshes[[i]]))
  dual[[i]]$cov1 = mesh_nodes %>% filter(mesh==i) %>% pull(cov1)
  dual[[i]]$cov2 = mesh_nodes %>% filter(mesh==i) %>% pull(cov2)
  dual[[i]]$cov3 = mesh_nodes %>% filter(mesh==i) %>% pull(cov3)
  
  int[[i]] = st_intersection(dual[[i]], poly)
  int[[i]]$weights = st_area(int[[i]])/st_area(poly)
  
  int[[i]] = int[[i]] %>%
    mutate(int1 = cov1*weights,
           int2 = cov2*weights,
           int3 = cov3*weights)
  int[[i]]$mesh = i
  
}

data.frame(mesh = paste("mesh",1:6,sep = ""),
  cov1 = sapply(int, function(x) sum(x$int1, na.rm = T)),
  cov2 = sapply(int, function(x) sum(x$int2, na.rm = T)),
  cov3 = sapply(int, function(x) sum(x$int3, na.rm = T))) %>%
  pivot_longer(-mesh) %>%
  ggplot() + geom_line(aes(mesh, value, group = name, color = name))

ggsave(paste(plot_save_dir, "integrated_cov.png",sep = ""))


p2 = rbind(cbind(dual[[1]], mesh = "Mesh 1"),
      cbind(dual[[2]], mesh = "Mesh 2"),
      cbind(dual[[3]], mesh = "Mesh 3"),
      cbind(dual[[4]], mesh = "Mesh 4"),
      cbind(dual[[5]], mesh = "Mesh 5"),
      cbind(dual[[6]], mesh = "Mesh 6")) %>%
  pivot_longer(-c(geometry, mesh)) %>%
  ggplot() + geom_sf(aes(fill = value)) +
  facet_grid(name~mesh) + theme_maps

p1 = ggplot() + geom_spatraster(data = cov1)+ theme_maps +
  coord_equal() + ggtitle("Cov1") +
ggplot() + geom_spatraster(data = cov2)+ theme_maps +
  coord_equal() + ggtitle("Cov2") +
ggplot() + geom_spatraster(data = cov3)+ theme_maps +
  coord_equal() + ggtitle("Cov3") + plot_layout(ncol = 3)

p1
p2


ints = rbind(int[[1]], int[[2]], int[[3]], int[[4]], int[[5]], int[[6]])
#+ fig.width=7, fig.height=7
ggplot() + geom_sf(data = mesh_nodes, aes(color = cov1)) + facet_wrap(.~mesh)  +
  xlim(c(0,500)) + ylim(c(0,500)) +
  theme_maps + theme(legend.position = "none")+
  geom_sf(data = poly, alpha = 0) + 
  geom_sf(data = ints, alpha = 0) + 
  ggtitle("Covariate 1")
ggsave(paste(plot_save_dir,"cov1_points.png", sep = ""))


#+ fig.width=7, fig.height=7
ggplot() + geom_sf(data = mesh_nodes, aes(color = cov2)) + facet_wrap(.~mesh)  +
  xlim(c(0,500)) + ylim(c(0,500)) + 
  theme_maps + theme(legend.position = "none")+
  geom_sf(data = poly, alpha = 0) + 
  geom_sf(data = ints, alpha = 0) + 
  
  theme(legend.position = "none") + ggtitle("Covariate 2")
ggsave(paste(plot_save_dir,"cov2_points.png", sep = ""))

#+ fig.width=7, fig.height=7
ggplot() + geom_sf(data = mesh_nodes, aes(color = cov3)) + facet_wrap(.~mesh)  +
  xlim(c(0,500)) + ylim(c(0,500)) + 
  theme_maps + theme(legend.position = "none") + 
  geom_sf(data = poly, alpha = 0) + 
  geom_sf(data = ints, alpha = 0) + 
  
  theme(legend.position = "none") + ggtitle("Covariate 3")
ggsave(paste(plot_save_dir,"cov3_points.png", sep = ""))


# simulate a point process ------------------------------------------------


simulate_PP = function(loglambda)
{
  wmax <- max(values(loglambda))
  Npoints <- rpois(1, lambda = max_dom^2 * exp(wmax))
  pointsX <- runif(n = Npoints, min = 0, max = max_dom)
  pointsY <- runif(n = Npoints, min = 0, max = max_dom)
  print(length(pointsX))
  
  s = extract(loglambda, st_as_sf(data.frame(x = pointsX, y = pointsY), coords = c("x","y")))
  lambda_ratio <- exp(s$val- wmax)
  keep = runif(Npoints) < lambda_ratio
  print(sum(keep))
  points = st_as_sf(data.frame(x = pointsX, y =  pointsY)[keep,],
                    coords = c("x","y"))
  return(points)
}


int1 = -6.5
int2 = -6.5
int3 = -6.5
beta = -0.8

loglambda1 = int1 + beta * cov1
loglambda2 = int2 + beta * cov2
loglambda3 = int3 + beta * cov3

points1 = simulate_PP(loglambda = loglambda1)
points2 = simulate_PP(loglambda = loglambda2)
points3 = simulate_PP(loglambda = loglambda3)





ggplot() + geom_spatraster(data = cov1)+
  geom_sf(data = points1, size = 0.3) +
  theme_maps + 
  scale_fill_scico(direction = -1) + 
  theme(legend.position = "none") +
ggplot() + geom_spatraster(data = cov2)+
  geom_sf(data = points2, size = 0.3) +
  theme_maps + 
  scale_fill_scico(direction = -1) + 
  theme(legend.position = "none") +
  
ggplot() + geom_spatraster(data = cov3)+
  geom_sf(data = points3, size = 0.3) +
  theme_maps + 
  scale_fill_scico(direction = -1) + 
  theme(legend.position = "none") +
  
  plot_layout(ncol = 2, ) +
  plot_annotation(
    title = 'Simulated PP'
  )

ggsave(paste(plot_save_dir, "simulatedPP.png"))

if(1)

  {
  data.frame( terra::extract(cov1, rbind(st_coordinates(points1), crds(cov1)))) %>% 
  mutate(p = c(rep(1, nrow(points1)), rep(0,dim(crds(cov1))[1]))) %>%
  group_by(p) %>%
  summarise(m = mean(val))

  data.frame( terra::extract(cov2, rbind(st_coordinates(points2), crds(cov2)))) %>% 
    mutate(p = c(rep(1, nrow(points2)), rep(0,dim(crds(cov2))[1]))) %>%
    group_by(p) %>%
    summarise(m = mean(val))


  data.frame( terra::extract(cov3, rbind(st_coordinates(points3), crds(cov3)))) %>% 
    mutate(p = c(rep(1, nrow(points3)), rep(0,dim(crds(cov3))[1]))) %>%
    group_by(p) %>%
    summarise(m = mean(val))
}



# FIT MODELS ----------------------------------------------------------------

bru_options_set(bru_verbose = 1)



## MESH 1  --------------------------------------------------

models1 = list()
models2 = list()
models3 = list()


cmp1 = ~ Intercept(1, model = "linear", prec.linear = 0.01)

for(i in 1:n_meshes)
{
  print(paste("MESH 1", " Int points ",i))
  lik = like(geometry ~ .,
             data = st_as_sf(points1),
             family = "cp",
             samplers = poly,
             ips = int_points[[i]],
             #domain = list(geometry = meshes[[i]])
             )
  
  models1[[i]] = bru(update.formula(cmp1, .~ . + 
                                      cov(cov1, model = "linear")), 
                     lik)
  
  
  lik = like(geometry ~ .,
             data = st_as_sf(points2),
             family = "cp",
             samplers = poly,
             ips = int_points[[i]],
             #domain = list(geometry = meshes[[i]])
             )
  
  models2[[i]] = bru(update.formula(cmp1, .~ . + 
                                           cov(cov2, model = "linear")), 
                          lik)
  
  lik = like(geometry ~ .,
             data = st_as_sf(points3),
             family = "cp",
             samplers = poly,
             ips = int_points[[i]],
             #domain = list(geometry = meshes[[i]])
             )
  
  models3[[i]] = bru(update.formula(cmp1, .~ . + 
                                      cov(cov3, model = "linear")), 
                     lik)
  
  
}





aa1 = purrr::map_dfr(models1,function(x) x$summary.fixed )[,c(1,3,5)]  %>%
  mutate(names  = rep(c("Int", "cov"), n_meshes),
         int_points = rep(paste("mesh", 1:n_meshes), each = 2),
         cov = "Cov1",
         true = rep(c(int1, beta),n_meshes ))


aa2 = purrr::map_dfr(models2,function(x) x$summary.fixed )[,c(1,3,5)]  %>%
  mutate(names  = rep(c("Int", "cov"), n_meshes),
         int_points = rep(paste("mesh", 1:n_meshes), each = 2),
         cov = "Cov2",
         true = rep(c(int2, beta),n_meshes ))


aa3 = purrr::map_dfr(models3,function(x) x$summary.fixed )[,c(1,3,5)]  %>%
  mutate(names  = rep(c("Int", "cov"), n_meshes),
         int_points = rep(paste("mesh", 1:n_meshes), each = 2),
         cov = "Cov3",
         true = rep(c(int3, beta),n_meshes ))


rbind(aa1, aa2,aa3 ) %>% 
  ggplot() + geom_point(aes(x = mean, y = int_points)) +
  geom_segment(aes(y = int_points, yend  = int_points, x = `0.025quant`, xend = `0.975quant`)) +
  geom_vline(aes(xintercept = true), linetype = "dashed") + 
  facet_grid(cov~ names, scales = "free") + xlab("") + ylab("")
ggsave(paste(plot_save_dir,"results.png",sep =""))

  
min_dist = numeric(n_meshes)
for(i in 1:n_meshes)
  min_dist[i] = min(as.vector(dist(meshes[[1]]$loc[,-3])))

range1/min_dist
range2/min_dist
range3/min_dist

