library(tidyverse)
library(inlabru)
library(patchwork)
library(INLA)
library(sf)
library(scico)
library(tidyterra)
library(scico)
library(terra)
library(dplyr)
library(ggplot2)

rm(list = ls())

set.seed(12455)

theme_maps = list(theme(axis.title.x=element_blank(),
                        axis.text.x=element_blank(),
                        axis.ticks.x=element_blank(),
                        axis.title.y=element_blank(),
                        axis.text.y=element_blank(),
                        axis.ticks.y=element_blank()),
                  scale_fill_scico(na.value = "transparent"),
                  scale_color_scico(na.value = "transparent"))

#### SS1 - Shifted meshes ####
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

## meshes
n_meshes = 6

xy0 = expand.grid(seq(-100,600,50),
                  seq(-100,600,50))

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


xy_list <- list()
for (i in 1:n_meshes) {
  xy <- xy0
  xy[,1] <- xy[,1] + cc[i,1]
  xy[,2] <- xy[,2] + cc[i,2]
  xy_list[[i]] <- data.frame(x = xy[,1], y = xy[,2], mesh = factor(i))
}
xy_all <- bind_rows(xy_list)
ggplot(xy_all, aes(x = x, y = y, color = mesh)) +
  geom_point(size = 1.5) +
  coord_equal() +
  labs(x = "X", y = "Y") +
  theme_minimal()

# meshes plot
ggplot() + gg(meshes[[1]]) +  geom_sf(data = poly, alpha = 0, color = "red") +
  coord_sf(xlim = c(-150,650), ylim = c(-150,650)) + ggtitle("Mesh 1") +
  xlab("") + ylab("") + theme_maps +
  ggplot() + gg(meshes[[2]]) +  geom_sf(data = poly, alpha = 0, color = "red") +
  coord_sf(xlim = c(-150,650), ylim = c(-150,650))+ ggtitle("Mesh 2") +
  xlab("") + ylab("") + theme_maps +
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

sapply(meshes, function(x)x$n)

int_points = list()
for(i in 1:n_meshes)
{
  print(paste("compute int. point  for mesh", i))
  int_points[[i]] = fm_int(meshes[[i]], samplers = poly)
}

## covariates
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
values(cov1) = val 
names(cov1) =  "val"

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

ggplot() + geom_spatraster(data = cov1) + geom_sf(data = poly, alpha = 0) +
  theme_maps + theme(legend.position  = "none") + #ggtitle("Covariate 1") +
  scale_fill_viridis_c(direction = -1) +
  ggplot() + geom_spatraster(data = cov2) + geom_sf(data = poly, alpha = 0) +
  theme_maps + theme(legend.position   = "none") + #ggtitle("Covariate 2") + 
  scale_fill_viridis_c(direction = -1) +
  ggplot() + geom_spatraster(data = cov3) + geom_sf(data = poly, alpha = 0) +
  theme_maps + theme(legend.position   = "none") + #ggtitle("Covariate 3") + 
  scale_fill_viridis_c(direction = -1) + plot_layout(ncol = 3)

## point patterns
set.seed(827372)
simulate_PP = function(loglambda)
{
  wmax <- max(values(loglambda))
  Npoints <- rpois(1, lambda = max_dom^2 * exp(wmax))
  pointsX <- runif(n = Npoints, min = 0, max = max_dom)
  pointsY <- runif(n = Npoints, min = 0, max = max_dom)
  print(length(pointsX))
  
  s = raster::extract(loglambda, st_as_sf(data.frame(x = pointsX, y = pointsY), coords = c("x","y")))
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

# cov + points figure
ggplot() + geom_spatraster(data = cov1)+
  geom_sf(data = points1, size = 1) +
  theme_maps + 
  scale_fill_viridis_c(direction = -1) + 
  theme(legend.position = "none") +
  ggplot() + geom_spatraster(data = cov2)+
  geom_sf(data = points2, size = 1) +
  theme_maps + 
  scale_fill_viridis_c(direction = -1) + 
  theme(legend.position = "none") +
  ggplot() + geom_spatraster(data = cov3)+
  geom_sf(data = points3, size = 1) +
  theme_maps + 
  scale_fill_viridis_c(direction = -1) + 
  theme(legend.position = "none") +
  plot_layout(ncol = 3)

## model fitting
bru_options_set(bru_verbose = 1)

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

## results of ss1 figure
rbind(aa1, aa2,aa3 ) %>% 
  ggplot() + geom_point(aes(x = mean, y = int_points)) +
  geom_segment(aes(y = int_points, yend  = int_points, x = `0.025quant`, xend = `0.975quant`)) +
  geom_vline(aes(xintercept = true), linetype = "dashed") + 
  facet_grid(cov~ names, scales = "free") + xlab("") + ylab("")

#### SS2 - 3 meshes + spde ####
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

## create meshes and int points
meshes = list()
xy0 = expand.grid(seq(-100,600,50),
                  seq(-100,600,50))
meshes[[1]] =  fm_mesh_2d_inla(loc = xy0)#fm_mesh_2d_inla(boundary = poly,
# max.edge = c(60,80))

#ggplot() + gg(meshes[[1]]) + coord_equal()
n_meshes = 3
for(i in 2:n_meshes)
  meshes[[i]] = fmesher:::fm_subdivide(meshes[[i-1]], n=2) #n=2?

sapply(meshes, function(x)x$n)

int_points = list()
for(i in 1:n_meshes)
{
  print(paste("compute int. point  for mesh", i))
  int_points[[i]] = fm_int(meshes[[i]], samplers = poly)
}

ggplot() + gg(meshes[[1]]) + ggtitle("Mesh 1") +
  coord_sf(xlim = c(-150,650), ylim = c(-150,650)) + xlab("") + ylab("") + theme_maps +
  ggplot() + gg(meshes[[2]]) + ggtitle("Mesh 2") +
  coord_sf(xlim = c(-150,650), ylim = c(-150,650)) + xlab("") + ylab("") + theme_maps +
  ggplot() + gg(meshes[[3]]) + ggtitle("Mesh 3") +
  coord_sf(xlim = c(-150,650), ylim = c(-150,650)) + xlab("") + ylab("") + theme_maps +
  plot_layout(ncol = 3)


# simulate a covariate 
# here we use a fine mesh to simulate the RF

sim_mesh <- inla.mesh.2d(loc.domain = poly, 
                         max.edge = c(5, 10))
sim_matern <- inla.spde2.pcmatern(sim_mesh,
                                  prior.sigma = c(1, 0.5),
                                  prior.range = c(1, 0.5))
A = inla.spde.make.A(mesh = sim_mesh, loc  = crds(cov1))

# covariate that varies in space with a short range!
covariate_ranges = c(10, 50, 600)
range = covariate_ranges[1]
sigma = 1
Q <- inla.spde2.precision(spde = sim_matern,
                          theta = c(log(range),log(sigma)))
samp1 <- inla.qsample(n = 2,
                      Q = Q,
                      mu = rep(0,dim(Q)[1]))
val =  as.vector((A%*%samp1[,1]))
val = round(val)
values(cov1) = val
names(cov1) =  "val"

# covariate that varies in space with a larger range!
range = covariate_ranges[2]
sigma = 1
Q <- inla.spde2.precision(spde = sim_matern,
                          theta = c(log(range),log(sigma)))
samp2 <- inla.qsample(n = 2,
                      Q = Q,
                      mu = rep(0,dim(Q)[1]))
val =  as.vector((A%*%samp2[,1]))
val = val
values(cov2) = val 
names(cov2) =  "val"

# covariate that varies very smoothly in space
range = covariate_ranges[3]
sigma = 1
Q <- inla.spde2.precision(spde = sim_matern,
                          theta = c(log(range),log(sigma)))
samp3 <- inla.qsample(n = 2,
                      Q = Q,
                      mu = rep(0,dim(Q)[1]))
val =  as.vector((A%*%samp3[,1]))
values(cov3) = val- mean(val)
names(cov3) =  "val"

# Simulate data from PP

# Simulate also a realization of the RF, here we choose a rather large range
range_matern = 300 #100
sigma_matern = 1
Q <- inla.spde2.precision(spde = sim_matern,
                          theta = c(log(range),log(sigma)))
samp <- inla.qsample(n = 1,
                     Q = Q,
                     mu = rep(0,dim(Q)[1]))

matern <- rast(nrows = max_dom/res, ncols=max_dom/res, 
               xmin=0, xmax=max_dom,
               ymin = 0, ymax = max_dom)
matern$val = as.vector((A%*%samp[,1]))

ggplot() + geom_spatraster(data = cov1) + geom_sf(data = poly, alpha = 0) +
  theme_maps + theme(legend.position  = "none") + #ggtitle("Covariate 1") +
  scale_fill_viridis_c(direction = -1) +
  ggplot() + geom_spatraster(data = cov2) + geom_sf(data = poly, alpha = 0) +
  theme_maps + theme(legend.position   = "none") + #ggtitle("Covariate 2") + 
  scale_fill_viridis_c(direction = -1) +
  ggplot() + geom_spatraster(data = cov3) + geom_sf(data = poly, alpha = 0) +
  theme_maps + theme(legend.position   = "none") + #ggtitle("Covariate 3") + 
  scale_fill_viridis_c(direction = -1) + 
  ggplot() + geom_spatraster(data = matern) + geom_sf(data = poly, alpha = 0) + 
  theme_maps + theme(legend.position = "none") + 
  scale_fill_viridis_c( direction = -1) +
  plot_layout(ncol = 4)

## point pattern
simulate_PP = function(loglambda)
{
  wmax <- max(values(loglambda))
  Npoints <- rpois(1, lambda = max_dom^2 * exp(wmax))
  pointsX <- runif(n = Npoints, min = 0, max = max_dom)
  pointsY <- runif(n = Npoints, min = 0, max = max_dom)
  print(length(pointsX))
  
  s = raster::extract(loglambda, st_as_sf(data.frame(x = pointsX, y = pointsY), coords = c("x","y")))
  lambda_ratio <- exp(s$val- wmax)
  keep = runif(Npoints) < lambda_ratio
  print(sum(keep))
  points = st_as_sf(data.frame(x = pointsX, y =  pointsY)[keep,],
                    coords = c("x","y"))
  return(points)
}


int1 = -7.5
int2 = -7.5
int3 = -7.5
beta = -1.3

loglambda1 = int1 + beta * cov1 + 1 * matern
loglambda2 = int2 + beta * cov2 + 1 * matern
loglambda3 = int3 + beta * cov3 + 1 * matern

points1 = simulate_PP(loglambda = loglambda1)
points2 = simulate_PP(loglambda = loglambda2)
points3 = simulate_PP(loglambda = loglambda3)

# covs, spde and points figure
ggplot() + geom_spatraster(data = beta * cov1 + matern)+
  geom_sf(data = points1, size = 1) +
  scale_fill_viridis_c() + 
  theme(legend.position = "none") +
  ggplot() + geom_spatraster(data = beta * cov2 + matern)+
  geom_sf(data = points2, size = 1) +
  scale_fill_viridis_c() + 
  theme(legend.position = "none") +
  ggplot() + geom_spatraster(data = beta * cov3 + matern)+
  geom_sf(data = points3, size = 1) +
  scale_fill_viridis_c() + 
  theme(legend.position = "none") +
  ggplot() + geom_spatraster(data = matern) + 
  scale_fill_viridis_c() + 
  theme(legend.position = "none") +
  plot_layout(ncol = 4)

## model fitting
bru_options_set(bru_verbose = 1)

# 3 meshes
# 3 covs + spde
# 3 pps
## MODEL 1 = COV + PP 1: rep for all 3 meshes
spde <- inla.spde2.pcmatern(mesh = meshes[[3]],          # change [[1/2/3]]
                            prior.range = c(200, 0.5), 
                            prior.sigma = c(1, 0.5))

cmp1 = ~ Intercept(1, model = "linear", prec.linear = 0.01) +  # change cmp1/2/3 and cov1/2/3
  matern(geometry, model = spde) + 
  cov(cov1, model = "linear", prec.linear = 0.01)

model1.3 <- bru(cmp1,                                     # change 1.1,1.2.,,,
                like(geometry ~.,                         # change cmp1/2/3
                     data = points1,                      # change points1/2/3
                     family = "cp", 
                     samplers = poly,
                     #domain = list(geometry = meshes[[3]])))
                     ips = int_points[[3]]))                   # change [[1/2/3]]

model1 <- list(model1.1, model1.2, model1.3) 
model2 <- list(model2.1, model2.2, model2.3)
model3 <- list(model3.1, model3.2, model3.3)

library(dplyr)
library(ggplot2)
library(purrr)
library(tibble)

# Combine fixed effect summaries
fixed_all <- map2_dfr(model3, 1:3, ~ {
  as.data.frame(.x$summary.fixed) %>%
    rownames_to_column("parameter") %>%
    mutate(mesh = paste0("Mesh ", .y))
})

# Combine hyperparameter summaries
hyper_all <- map2_dfr(model3, 1:3, ~ {
  as.data.frame(.x$summary.hyperpar) %>%
    rownames_to_column("parameter") %>%
    mutate(mesh = paste0("Mesh ", .y))
})

# Fixed effects plot
plot_fixed <- ggplot(fixed_all, aes(x = parameter, y = mean, color = mesh)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = `0.025quant`, ymax = `0.975quant`),
                width = 0.5, position = position_dodge(width = 0.5)) +
  geom_segment(aes(x = 0.7, xend = 1.3, y = -1.3, yend = -1.3), linetype = "dashed", color = "darkgray") +  # true beta
  geom_segment(aes(x = 1.7, xend = 2.3, y = -7.5, yend = -7.5), linetype = "dashed", color = "darkgray") +  # true intercept
  theme_minimal() +
  labs(title = "Point pattern 2: Fixed Effects", y = "Mean (95% CI)", x = "") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "none") 

# Hyperparameters plot
plot_hyper <- ggplot(hyper_all, aes(x = parameter, y = mean, color = mesh)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = `0.025quant`, ymax = `0.975quant`),
                width = 0.5, position = position_dodge(width = 0.5)) +
  geom_segment(aes(x = 0.7, xend = 1.3, y = 300, yend = 300), 
               linetype = "dashed", color = "darkgray") +  # true range
  geom_segment(aes(x = 1.7, xend = 2.3, y = 1, yend = 1), 
               linetype = "dashed", color = "darkgray") +  # true sd
  theme_minimal() +
  labs(title = "Hyperparameters", y = "", x = "") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) #legend.position = "none"

# Print or combine
pdf()
wrap_plots(plot_fixed, plot_hyper)
dev.off()
