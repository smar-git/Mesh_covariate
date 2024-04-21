#' Here we simulate three PP using  a LGCP where the log intensity
#' depends on only an intercept and a covariate
#' 
#' $$
#' \beta_0 + beta_1\ x(s)
#' $$
#' 
#' We then estimate the model using different meshes, all meshes have 
#' similar parameters but different locations for the mesh nodes.

#+ echo = FALSE
knitr::opts_chunk$set(collapse = TRUE, echo = FALSE,
                      message = FALSE, warning = FALSE,
                      fig.width=5, fig.height=5)



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

set.seed(123)
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

polys = list()
polys[[1]] = poly
polys[[2]]  = data.frame(x = box[c(1,3)]-5, y = box[c(2,4)]+10) %>% 
  st_as_sf(coords = c("x", "y")) %>% 
  st_bbox() %>% 
  st_as_sfc()
polys[[3]]  = data.frame(x = box[c(1,3)]+10, y = box[c(2,4)]+10) %>% 
  st_as_sf(coords = c("x", "y")) %>% 
  st_bbox() %>% 
  st_as_sfc()
polys[[4]]  = data.frame(x = box[c(1,3)]+20, y = box[c(2,4)]-20) %>% 
  st_as_sf(coords = c("x", "y")) %>% 
  st_bbox() %>% 
  st_as_sfc()
polys[[5]]  = data.frame(x = box[c(1,3)]+7, y = box[c(2,4)]+10) %>% 
  st_as_sf(coords = c("x", "y")) %>% 
  st_bbox() %>% 
  st_as_sfc()
polys[[6]]  = data.frame(x = box[c(1,3)]-20, y = box[c(2,4)]+0) %>% 
  st_as_sf(coords = c("x", "y")) %>% 
  st_bbox() %>% 
  st_as_sfc()

cov1 <- rast(nrows = max_dom/res, ncols=max_dom/res, xmin=0, xmax=max_dom, ymin = 0, ymax = max_dom)
cov2 <- rast(nrows = max_dom/res, ncols=max_dom/res, xmin=0, xmax=max_dom, ymin = 0, ymax = max_dom)
cov3 <- rast(nrows = max_dom/res, ncols=max_dom/res, xmin=0, xmax=max_dom, ymin = 0, ymax = max_dom)



n_meshes = 6



## create meshes and int points -----------------------------------------------------------
meshes = list()

for(i in 1:n_meshes)
    meshes[[i]] =  fm_mesh_2d_inla(boundary = polys[[i]],
                               max.edge = c(75,80))




ggplot() + gg(meshes[[1]]) +  geom_sf(data = poly, alpha = 0, color = "red") +
  coord_sf(xlim = c(-100,600), ylim = c(-100,600)) + ggtitle("Mesh 1") +
  xlab("") + ylab("") +
ggplot() + gg(meshes[[2]]) +  geom_sf(data = poly, alpha = 0, color = "red")+
  coord_sf(xlim = c(-100,600), ylim = c(-100,600))+ ggtitle("Mesh 2") +
  xlab("") + ylab("") +
ggplot() + gg(meshes[[3]]) +  geom_sf(data = poly, alpha = 0, color = "red")+
  coord_sf(xlim = c(-100,600), ylim = c(-100,600))+ ggtitle("Mesh 3") +
  xlab("") + ylab("") +
ggplot() + gg(meshes[[4]]) +  geom_sf(data = poly, alpha = 0, color = "red")+
  coord_sf(xlim = c(-100,600), ylim = c(-100,600))+ ggtitle("Mesh 4") +
  xlab("") + ylab("") +
ggplot() + gg(meshes[[5]]) +  geom_sf(data = poly, alpha = 0, color = "red")+
  coord_sf(xlim = c(-100,600), ylim = c(-100,600))+ ggtitle("Mesh 5") +
  xlab("") + ylab("") +
ggplot() + gg(meshes[[6]]) +  geom_sf(data = poly, alpha = 0, color = "red")+
  coord_sf(xlim = c(-100,600), ylim = c(-100,600))+ ggtitle("Mesh 6") +
  xlab("") + ylab("") 

 


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


range = 10
sigma = 1
Q <- inla.spde2.precision(spde = sim_matern,
                          theta = c(log(range),log(sigma)))
samp1 <- inla.qsample(n = 2,
                     Q = Q,
                     mu = rep(0,dim(Q)[1]))
val =  as.vector((A%*%samp1[,1]))
#val =   ifelse(val<1.5,0,val)
values(cov1) = val #- mean(val[val>0])
names(cov1) =  "val"

#ggplot() + geom_spatraster(data = cov1) + scale_fill_scico( direction = -1) +
#  gg(meshes[[1]])


range = 50
sigma = 1
Q <- inla.spde2.precision(spde = sim_matern,
                          theta = c(log(range),log(sigma)))
samp2 <- inla.qsample(n = 2,
                      Q = Q,
                      mu = rep(0,dim(Q)[1]))
val =  as.vector((A%*%samp2[,1]))
#val = ifelse(val<1,0,val)
values(cov2) = val# - mean(val[val>0])
names(cov2) =  "val"
#ggplot() + geom_spatraster(data = cov2) + scale_fill_scico( direction = -1) +
#  gg(meshes[[1]]) + coord_equal()


range = 600
sigma = 1
Q <- inla.spde2.precision(spde = sim_matern,
                          theta = c(log(range),log(sigma)))
samp3 <- inla.qsample(n = 2,
                     Q = Q,
                     mu = rep(0,dim(Q)[1]))
val =  as.vector((A%*%samp3[,1]))
values(cov3) = val - mean(val)
names(cov3) =  "val"

#ggplot() + geom_spatraster(data = cov3) + scale_fill_scico( direction = -1) +
#  gg(meshes[[1]])

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
  
  plot_layout(ncol = 2)


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
  facet_grid(cov~ names, scales = "free") 


  
