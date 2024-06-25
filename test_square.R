
#' Here we simulate three PP using  a LGCP where the log intensity
#' depends on only an intercept and a covariate
#' 
#' $$
#' \beta_0 + beta_1\ x(s)
#' $$
#' 

#+ echo = FALSE
knitr::opts_chunk$set(collapse = TRUE, echo = FALSE,
                      message = FALSE, warning = FALSE)

library(tidyverse)
library(inlabru)
library(patchwork)
library(INLA)
library(sf)
library(scico)
library(tidyterra)
library(terra)
library(fmesher)

rm(list = ls())


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

box = c(xmin = 0 , 
        ymin = 0,
        xmax = max_dom,
        ymax = max_dom)
poly <- data.frame(x = box[c(1,3)], y = box[c(2,4)]) %>% 
  st_as_sf(coords = c("x", "y")) %>% 
  st_bbox() %>% 
  st_as_sfc()

cov1 <- rast(nrows = max_dom/res, ncols=max_dom/res, xmin=0, xmax=max_dom, ymin = 0, ymax = max_dom)
cov2 <- rast(nrows = max_dom/res, ncols=max_dom/res, xmin=0, xmax=max_dom, ymin = 0, ymax = max_dom)
cov3 <- rast(nrows = max_dom/res, ncols=max_dom/res, xmin=0, xmax=max_dom, ymin = 0, ymax = max_dom)


# create meshes and int points -----------------------------------------------------------
meshes = list()
meshes[[1]] =  fm_mesh_2d_inla(boundary = poly,
                               max.edge = c(95,100))


ggplot() + gg(meshes[[1]]) + coord_equal()

n_meshes = 4

for(i in 2:n_meshes) {
  meshes[[i]] = fmesher:::fm_subdivide(meshes[[i-1]])
  
  
  sapply(meshes, function(x)x$n) %>% plot()

} 

p1 <- ggplot() + gg(meshes[[1]]) + coord_equal()
p2 <- ggplot() + gg(meshes[[2]]) + coord_equal()
p3 <- ggplot() + gg(meshes[[3]]) + coord_equal()
p4 <- ggplot() + gg(meshes[[4]]) + coord_equal()

gridExtra::grid.arrange(p1, p2, p3, p4, nrow = 2)

int_points = list()

for(i in 1:n_meshes){
  print(paste("compute int. point  for mesh", i))
  int_points[[i]] = fm_int(meshes[[i]], samplers = poly)
}


# simulate a covariate ----------------------------------------------------


sim_mesh <- fm_mesh_2d_inla(boundary = poly, 
                            max.edge = c(5, 10))

sim_matern <- inla.spde2.pcmatern(sim_mesh,
                                  prior.sigma = c(1, 0.5),
                                  prior.range = c(1, 0.5))
A = inla.spde.make.A(mesh = sim_mesh, loc  = crds(cov1))

## cov 1, short range
range = 8
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

ggplot() + 
  geom_spatraster(data = cov1) + 
  scale_fill_scico(palette = 'imola', direction = -1) +
  gg(meshes[[1]]) + 
  coord_equal()

## cov 2, medium range
range = 40
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

ggplot() + 
  geom_spatraster(data = cov2) + 
  scale_fill_scico(palette = 'imola', direction = -1) +
  gg(meshes[[1]]) + 
  coord_equal()

## cov 3, long range
range = 600
sigma = 1

Q <- inla.spde2.precision(spde = sim_matern,
                          theta = c(log(range),log(sigma)))

samp3 <- inla.qsample(n = 2,
                     Q = Q,
                     mu = rep(0,dim(Q)[1]))

val =  as.vector((A%*%samp3[,1]))

values(cov3) = val -  mean(val)

names(cov3) =  "val"

ggplot() + 
  geom_spatraster(data = cov3) + 
  scale_fill_scico(palette = 'imola', direction = -1) +
  gg(meshes[[1]]) + 
  coord_equal()

# simulate a point process ------------------------------------------------


simulate_PP = function(loglambda) {
  
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
beta = -1

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


if(0) {
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

bru_options_set(bru_verbose = 2)



## MESH 1  --------------------------------------------------

models1 = list()
models2 = list()
models3 = list()

cmp1 = ~ Intercept(1, model = "linear", prec.linear = 0.01)

for(i in 1:n_meshes) {
  
  print(paste("MESH 1", " Int points ",i))
  lik = like(geometry ~ .,
             data = st_as_sf(points1),
             family = "cp",
             samplers = poly,
             ips = int_points[[i]],
             #domain = list(geometry = meshes[[3]])
             )
  
  models1[[i]] = bru(update.formula(cmp1, .~ . + 
                                      cov(cov1, model = "linear")), 
                     lik)
  
  
  lik = like(geometry ~ .,
             data = st_as_sf(points2),
             family = "cp",
             samplers = poly,
             ips = int_points[[i]],
             #domain = list(geometry = meshes[[3]])
             )
  
  models2[[i]] = bru(update.formula(cmp1, .~ . + 
                                           cov(cov2, model = "linear")), 
                          lik)
  
  lik = like(geometry ~ .,
             data = st_as_sf(points3),
             family = "cp",
             samplers = poly,
             ips = int_points[[i]],
             #domain = list(geometry = meshes[[3]])
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


rbind(aa1, aa2,aa3 ) %>% 
  filter(names=="Int") %>%
  ggplot() + geom_point(aes(x = mean, y = int_points)) +
  geom_segment(aes(y = int_points, yend  = int_points, x = `0.025quant`, xend = `0.975quant`)) +
  geom_vline(aes(xintercept = true), linetype = "dashed") + 
  facet_grid(cov~ ., scales = "free") +
  coord_cartesian(xlim = c(-7,-3)) +
  
  ggtitle("Intercept")

  
rbind(aa1, aa2, aa3 ) %>% 
  filter(names=="cov") %>%
  ggplot() + geom_point(aes(x = mean, y = int_points)) +
  geom_segment(aes(y = int_points, yend  = int_points, x = `0.025quant`, xend = `0.975quant`)) +
  geom_vline(aes(xintercept = true), linetype = "dashed") + 
  facet_grid(cov~ ., scales = "free") +
  ggtitle("covariate")
