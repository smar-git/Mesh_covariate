#' Simulate PP data from  a LGCP
#' The log likelihood is defined as
#' 
#' $$
#' \beta0 + beta_1 * x(s) + u(s)
#' $$
#' 
#' where   $u(s)$ is a SPDE model
#' 
#' The data are simulated using 3 different covariates with different 
#' spatial characteristics
#' 
#' The model is then estimated using different meshes


#+ setup, include=FALSE
knitr::opts_chunk$set(collapse = TRUE, echo = FALSE,
                      message = FALSE, warning = FALSE)


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

#' ## Create the meshes and compute integration points
#+

## create meshes and int points -----------------------------------------------------------
meshes = list()
meshes[[1]] =  fm_mesh_2d_inla(boundary = poly,
                               max.edge = c(60,100))

#ggplot() + gg(meshes[[1]]) + coord_equal()
n_meshes = 4
for(i in 2:n_meshes)
  meshes[[i]] = fmesher:::fm_subdivide(meshes[[i-1]])


sapply(meshes, function(x)x$n)

ggplot() + gg(meshes[[n_meshes]]) + coord_equal()

int_points = list()
for(i in 1:n_meshes)
{
  print(paste("compute int. point  for mesh", i))
  int_points[[i]] = fm_int(meshes[[i]], samplers = poly)
}

#' ## Simulate covariates and matern field
#' 
#' 
#' 
#+


# simulate a covariate ----------------------------------------------------

# here we use a fine mesh to simulate the RF

sim_mesh <- inla.mesh.2d(loc.domain = poly, 
                         max.edge = c(5, 10))
sim_matern <- inla.spde2.pcmatern(sim_mesh,
                                  prior.sigma = c(1, 0.5),
                                  prior.range = c(1, 0.5))
A = inla.spde.make.A(mesh = sim_mesh, loc  = crds(cov1))


# covariate that varies in space with a short range!
range = 7
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

ggplot() + geom_spatraster(data = cov1) + scale_fill_scico( direction = -1) +
  gg(meshes[[1]]) + coord_equal() + ggtitle("Covariate 1") + theme_maps


# covariate that varies in space with a larger range!
range = 40
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
ggplot() + geom_spatraster(data = cov2) + scale_fill_scico( direction = -1) +
  gg(meshes[[1]]) + coord_equal()+ ggtitle("Covariate 2") + theme_maps

# covariate that varies very smoothly in space
range = 600
sigma = 1
Q <- inla.spde2.precision(spde = sim_matern,
                          theta = c(log(range),log(sigma)))
samp3 <- inla.qsample(n = 2,
                     Q = Q,
                     mu = rep(0,dim(Q)[1]))
val =  as.vector((A%*%samp3[,1]))
values(cov3) = val- mean(val)
names(cov3) =  "val"

ggplot() + geom_spatraster(data = cov3) + scale_fill_scico( direction = -1) +
  gg(meshes[[1]])+ coord_equal()+ ggtitle("Covariate 3") + theme_maps


#' \newpage
#' ## Simulate data from PP
#+


# Simulate also a realization of the RF, here we choose a rather 
# large range
range_matern = 300
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

ggplot() + geom_spatraster(data = matern) + scale_fill_scico( direction = -1) +
  gg(meshes[[1]])+ coord_equal()+ ggtitle("Matern field ") + theme_maps



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


int1 = -7
int2 = -7
int3 = -7
beta = -1.3

loglambda1 = int1 + beta * cov1 + 1 * matern
loglambda2 = int2 + beta * cov2 + 1 * matern
loglambda3 = int3 + beta * cov3 + 1 * matern

points1 = simulate_PP(loglambda = loglambda1)
points2 = simulate_PP(loglambda = loglambda2)
points3 = simulate_PP(loglambda = loglambda3)





ggplot() + geom_spatraster(data = beta * cov1 + matern)+
  geom_sf(data = points1) +
  theme_maps + 
  scale_fill_scico(direction = -1) + 
  theme(legend.position = "none") +
ggplot() + geom_spatraster(data = beta * cov2 + matern)+
  geom_sf(data = points2) +
  theme_maps + 
  scale_fill_scico(direction = -1) + 
  theme(legend.position = "none") +
  
ggplot() + geom_spatraster(data = beta * cov3 + matern)+
  geom_sf(data = points3) +
  theme_maps + 
  scale_fill_scico(direction = -1) + 
  theme(legend.position = "none") +

ggplot() + geom_spatraster(data = matern) +   
  theme_maps + 
  scale_fill_scico(direction = -1) + 
  theme(legend.position = "none") +
    plot_layout(ncol = 2)


if(0)
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

#' \newpage
#' ## Fit the models

# FIT MODELS ----------------------------------------------------------------

bru_options_set(bru_verbose = 1)


# COVARIATE 3 -------------------------------------------------------------
models3 = list()

for(idx_mesh in 1:n_meshes)
{
  print(idx_mesh)
  spde =  inla.spde2.pcmatern(mesh = meshes[[idx_mesh]],
                              # PC-prior on range: P(practic.range < 0.05) = 0.01
                              prior.range = c(200, 0.5),
                              # PC-prior on sigma: P(sigma > 1) = 0.01
                              prior.sigma = c(1, 0.5)) 
  
  cmp3 = ~ Intercept(1, model = "linear", prec.linear = 0.01) +
    matern(geometry, model= spde) + 
    cov(cov3, model = "linear", prec.linear = 0.01)
  
  models3[[idx_mesh]] = list()
  for(j in idx_mesh:n_meshes)
  {
    lik = like(geometry ~ .,
               data = points3,
               family = "cp",
               samplers = poly,
               ips = int_points[[j]])
    models3[[idx_mesh]][[j]] = bru(cmp3, lik)
  }
}


# COVARIATE 2 -------------------------------------------------------------
models2 = list()

for(idx_mesh in 1:n_meshes)
{
  print(idx_mesh)
  spde =  inla.spde2.pcmatern(mesh = meshes[[idx_mesh]],
                              # PC-prior on range: P(practic.range < 0.05) = 0.01
                              prior.range = c(200, 0.5),
                              # PC-prior on sigma: P(sigma > 1) = 0.01
                              prior.sigma = c(1, 0.5)) 
  
  cmp2 = ~ Intercept(1, model = "linear", prec.linear = 0.01) +
    matern(geometry, model= spde) + 
    cov(cov2, model = "linear", prec.linear = 0.01)
  
  models2[[idx_mesh]] = list()
  for(j in idx_mesh:n_meshes)
  {
    lik = like(geometry ~ .,
               data = points2,
               family = "cp",
               samplers = poly,
               ips = int_points[[j]])
    models2[[idx_mesh]][[j]] = bru(cmp2, lik)
  }
}

# COVARIATE 1 -------------------------------------------------------------
models1 = list()

for(idx_mesh in 1:n_meshes)
{
  print(idx_mesh)
  spde =  inla.spde2.pcmatern(mesh = meshes[[idx_mesh]],
                              # PC-prior on range: P(practic.range < 0.05) = 0.01
                              prior.range = c(200, 0.5),
                              # PC-prior on sigma: P(sigma > 1) = 0.01
                              prior.sigma = c(1, 0.5)) 
  
  cmp1 = ~ Intercept(1, model = "linear", prec.linear = 0.01) +
    matern(geometry, model= spde) + 
    cov(cov1, model = "linear", prec.linear = 0.01)
  models1[[idx_mesh]] = list()
  for(j in idx_mesh:n_meshes)
  {
    lik = like(geometry ~ .,
             data = points1,
             family = "cp",
             samplers = poly,
             ips = int_points[[j]])
    models1[[idx_mesh]][[j]] = bru(cmp1, lik)
  }
}


# RESULTS -----------------------------------------------------------------

#' ## Integration points
#' 
#+ 
list1 = vector("list",n_meshes*n_meshes )

kk =1
for(ii in 1:n_meshes)
  for(jj in 1:n_meshes)
  {
    if(jj<ii)
    {
      list1[[kk]] = plot_spacer()
      kk = kk + 1
    }
    else{
      list1[[kk]] = ggplot() + gg(meshes[[ii]]) + 
        #ggtitle(paste("int_point", idx_mesh)) +
        coord_equal() + theme_maps + 
        geom_sf(data = int_points[[jj]], aes(color = weight)) +
        theme(legend.position  = "none") 
      kk = kk + 1
    }
  }
wrap_plots(list1)  



#'\newpage
#' ## Estimated Parameters
#' 
#+

create_plots = function(mods)
{
  plots_int = list()
  plots_cov = list()
  plots_range = list()
  plots_sd = list()
  for(idx_mesh in 1:n_meshes)
  { 
    err = purrr::map(mods[[idx_mesh]],function(x) x$error )
    null = purrr::map_int(mods[[idx_mesh]],function(x) !is.null(x) )
    err = data.frame(idx = 1:n_meshes, valid = sapply(err,function(x) is.null(x)))
    err$null = null
    err = err$idx[null>0 & err$valid]
    
    
    df_fixed = data.frame(idx = rep(1:n_meshes,each = 2), mean = NA, q1 = NA, q2 = NA)
    df_hyper = data.frame(idx = rep(1:n_meshes,each = 2), mean = NA, q1 = NA, q2 = NA)
    
    # fixed effects
    
    vals = purrr::map_dfr(mods[[idx_mesh]][err],function(x) x$summary.fixed )[,c(1,3,5)] 
    df_fixed[df_fixed$idx %in% err,-1] = vals
    
    df_fixed = df_fixed %>%
      mutate(param = rep(c("int", "cov"), n_meshes),
             mesh = rep(1:n_meshes, each = 2),
             mesh_spde = rep(paste("int ",1:n_meshes),each = 2),
             true = rep(c(int2,beta),n_meshes)) 
    
    p1 = df_fixed %>% filter(param == "int") %>%
        ggplot() + geom_segment(aes(y = mesh_spde, yend = mesh_spde,
                                  x = q1, xend = q2)) +
      geom_point(aes(y = mesh_spde, x = mean)) +
      geom_vline(aes(xintercept = true)) +
        xlab("") + ylab("") 
      p1
      
      p2 = df_fixed %>% filter(param == "cov") %>%
        ggplot() + geom_segment(aes(y = mesh_spde, yend = mesh_spde,
                                    x = q1, xend = q2)) +
        geom_point(aes(y = mesh_spde, x = mean)) +
        geom_vline(aes(xintercept = true)) +
        xlab("") + ylab("") 
    p2   
   plots_cov[[idx_mesh]] =    p2
   plots_int[[idx_mesh]] =    p1
   
   
   # hyperpar
   
   
   vals = purrr::map_dfr(mods[[idx_mesh]][err],function(x) x$summary.hyperpar )[,c(1,3,5)] 
   df_hyper[df_hyper$idx %in% err,-1] = vals
   
   df_hyper = df_hyper %>%
     mutate(param = rep(c("range", "sd"), n_meshes),
            mesh = rep(1:n_meshes, each = 2),
            mesh_spde = rep(paste("int ",1:n_meshes),each = 2),
            true = rep(c(int2,beta),n_meshes)) 
   
   p3 = df_hyper %>% filter(param == "range") %>%
     ggplot() + geom_segment(aes(y = mesh_spde, yend = mesh_spde,
                                 x = q1, xend = q2)) +
     geom_point(aes(y = mesh_spde, x = mean)) +
     geom_vline(aes(xintercept = range_matern)) +
     xlab("") + ylab("") 
   p3
   
   p4 = df_hyper %>% filter(param == "sd") %>%
     ggplot() + geom_segment(aes(y = mesh_spde, yend = mesh_spde,
                                 x = q1, xend = q2)) +
     geom_point(aes(y = mesh_spde, x = mean)) +
     geom_vline(aes(xintercept = sigma_matern)) +
     xlab("") + ylab("") 
   p4 
   plots_range[[idx_mesh]] =    p3
   plots_sd[[idx_mesh]] =    p4
  }
return(list(covariate  = plots_cov,
        intercept = plots_int,
        range = plots_range,
        sd = plots_sd))
}

plots1 = create_plots(models1)
plots2 = create_plots(models2)
plots3 = create_plots(models3)

wrap_plots(c(plots1$covariate, plots2$covariate, plots3$covariate))+
  plot_layout(axes ="collect_y", tag_level = "new") +
  plot_annotation(title = "Estimated Covariate")



wrap_plots(c(plots1$intercept, plots2$intercept, plots3$intercept))+
  plot_layout(axes ="collect_y") +
  plot_annotation(title = "Estimated Intercept")

wrap_plots(c(plots1$range, plots2$range, plots3$range))+
  plot_layout(axes ="collect_y") +
  plot_annotation(title = "Estimated Range")

wrap_plots(c(plots1$sd, plots2$sd, plots3$sd))+
  plot_layout(axes ="collect_y") +
  plot_annotation(title = "Estimated Sd")













































































































































