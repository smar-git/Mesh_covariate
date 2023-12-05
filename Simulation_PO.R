#' ---
#' title: "Simulation study about covariates and PP modeling in INLA"
#' author: "Sara Martino"
#' ---
#' 
#' We simulate data from a poisson process with two different covariates, one that is smooth 
#' in space and one that is not smooth.
#' We then estimate the PP model using 5 different meshes going from fine to coarses ones.





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
source("dual_mesh.R")
#set.seed(10055)
set.seed(4567)


# read data ---------------------------------------------------------------
load(here::here("Data/data_serotine1.RData"))

# read shapefile ----------------------------------------------------------
outline = st_simplify(st_as_sf(outline), dTolerance = 3)

box = c(xmin = 134.8978 , 
        ymin = 13.5963,
        xmax = 655.9502,
        ymax = 975.3520 )
      #  ymax = 300)

outline = st_crop(outline, box)


poly <- data.frame(x = box[c(1,3)], y = box[c(2,4)]) %>% 
  st_as_sf(coords = c("x", "y"), 
           crs = st_crs(outline)) %>% 
  st_bbox() %>% 
  st_as_sfc()

theme_maps = list(theme(axis.title.x=element_blank(),
                        axis.text.x=element_blank(),
                        axis.ticks.x=element_blank(),
                        axis.title.y=element_blank(),
                        axis.text.y=element_blank(),
                        axis.ticks.y=element_blank()),
                  scale_fill_scico(na.value = "transparent"),
                  scale_color_scico(na.value = "transparent"),
                  coord_sf(xlim = box[c(1,3)],
                           ylim = box[c(2,4)]) )
#' ## Read covariates
#' there are two covariates. One is smooth the other is not
#+

# read covariates ---------------------------------------------------------
cov = terra::rast(here::here("covariates/landuse_spatRaster.tif"))
Temp = terra::rast(here::here("covariates/temperature.tif"))
Temp = terra::project(Temp,crs(outline))
names(Temp) = "temp"
values(Temp) = scale(values(Temp))
Temp = mask(Temp,outline)
Temp  = focal(Temp, w=15, fun="mean", 
              expand = TRUE, na.rm = T)

Broadlead = cov$Broadleaf
values(Broadlead) = scale(values(Broadlead))
Broad_smooth = focal(Broadlead, w=13, fun="mean", 
                     expand = TRUE,
                     na.rm =T)
Temp = (project(Temp, Broad_smooth))



p1 = ggplot() + geom_spatraster(data = Temp) +
  geom_sf(data= outline, alpha = 0)+ theme_maps+
  ggtitle("Temp") 
p2 =   ggplot() + geom_spatraster(data = Broadlead) +
  geom_sf(data= outline, alpha = 0)+ 
  theme_maps+
  ggtitle("Broadlead") 


p2a =   ggplot() + geom_spatraster(data = Broad_smooth) +
  geom_sf(data= outline, alpha = 0)+
  theme_maps+
  ggtitle("Broadlead smoothed")



p1 + p2 + p2a + plot_annotation(
  title = 'used covariates',
  caption = 'covariates used in the simulation'
)


#' ## Simulate intensiity surfaces
#' These are simulated as:
#' $$
#' \log(\lambda(s)) = \beta_0 + \beta_1 \text{Cov}(s) + \omega(s)
#' $$
#' where $\omega(s)$ is a SPDE field
#' The two simulated surfaces have the same gaussian field but each uses one covariate
#+

# create mesh for simulations-------------------------------------------------------------

sim_mesh <- fm_mesh_2d_inla(boundary = outline,
                        max.edge = c(5,150),
                        cutoff = 5,
                        crs= crs(outline))
matern <- inla.spde2.pcmatern(sim_mesh,
                              prior.sigma = c(1, 0.5),
                              prior.range = c(100, 0.5))


# simulate a random field -------------------------------------------------

true_range = 180
true_sd = 0.6
Q <- inla.spde2.precision(spde = matern,
                          theta = c(log(true_range),log(true_sd)))
samp <- inla.qsample(n = 1,
                     Q = Q,
                     mu = rep(0,dim(Q)[1]))


points = terra::crds(Broadlead, na.rm = F)
A = inla.spde.make.A(mesh = sim_mesh, points )

true_GF = Broadlead
names(true_GF) = "spde"
values(true_GF) = (A%*%samp)[,1]

beta = 0.8
int = -6.5

#' We set the value of the intercept and the covariate to
#' $\beta_0 =$`r beta` and $\beta_1 =$`r int`.
#' 
log_lambda0 = int + true_GF+
  beta * Broadlead


log_lambda1 = int + true_GF+
  beta * Temp

p1 = ggplot() + geom_spatraster(data = (log_lambda0))+
  theme_maps+
  ggtitle("lambda with broadlead")
p2 = ggplot() + geom_spatraster(data = (log_lambda1))+
  theme_maps+
  ggtitle("lambda with temp")
p3 = ggplot() + geom_spatraster(data = true_GF)+
  geom_sf(data=outline, alpha = 0) +
  theme_maps+
    ggtitle("latent field")
p1 + p2 + p3 + plot_annotation(title = "True log intensity fields")

#' ## Simulate data
#' We create two datasets one for each simulated intensity surfaces
#+


# simulate PP -------------------------------------------------------------


mat0 =  (matrix(values(log_lambda0)[,1],
               nrow = length(unique(points[,1]))))
for(i in 1:dim(mat0)[1])
  mat0[i,] = rev(mat0[i,])
mat0 = t(mat0)
pp0 = spatstat.random::rpoispp(spatstat.geom::im(exp(mat0)), nsim = 1)
#plot(spatstat.geom::im((mat0)))

mat1 =  (matrix(values(log_lambda1)[,1],
                nrow = length(unique(points[,1]))))
for(i in 1:dim(mat1)[1])
  mat1[i,] = rev(mat1[i,])
mat1 = t(mat1)
pp1 = spatstat.random::rpoispp(spatstat.geom::im(exp(mat1)), nsim = 1)
#plot(spatstat.geom::im((mat1)))


PO0 = data.frame(x = pp0$x, y = pp0$y) %>% st_as_sf(coords = c("x","y"))
st_crs(PO0) = st_crs(outline)
PO0 = st_intersection(PO0, outline)

PO1 = data.frame(x = pp1$x, y = pp1$y) %>% st_as_sf(coords = c("x","y"))
st_crs(PO1) = st_crs(outline)
PO1 = st_intersection(PO1, outline)

ggplot() + geom_sf(data = outline) + 
  geom_spatraster(data = log_lambda0) + 
  geom_sf(data = PO0, size = 0.5) +
  ggtitle("PP with broadlead") + theme_maps+
  ggplot() + geom_sf(data = outline) +
  geom_spatraster(data = log_lambda1) + 
  
  geom_sf(data = PO1, size = 0.5) +
  ggtitle("PP with temperature") +
theme_maps


# model fit----------------------------------------------------------------

bru_options_set(bru_verbose = 2)


models0 = list()
models0a = list()
models1 = list()


## create meshes -----------------------------------------------------------
edges = c(15,25,45,60)
cutoff = c(5,5,20,20)
meshes = list()
for(i in 1:4)
  meshes[[i]] = fm_mesh_2d_inla(boundary = outline,
                                max.edge = c(edges[i],80),
                                cutoff = cutoff[i],
                                crs= crs(outline))


n_meshes = 4
## models for temperature --------------------------------------------------
for(i in 1:n_meshes)
{
     matern <-  inla.spde2.pcmatern(meshes[[i]],
                                  prior.sigma = c(1, 0.5),
                                  prior.range = c(100, .5))
    cmp_PO1 = ~ InterceptPO(1, model = "linear", prec.linear = 0.01) +
      maternPO(geometry, model = matern ) +
      cov(Temp, model = "linear")

    likPO1 = like(geometry ~ .,
                 data = PO1,
                 family = "cp",
                 samplers = outline,
                 domain = list(geometry = meshes[[i]]))

  models1[[i]] = bru(cmp_PO1, likPO1)
}


## models with broadlead ---------------------------------------------------

for(i in 1:n_meshes)
{
  matern <- inla.spde2.pcmatern(meshes[[i]],
                                prior.sigma = c(1, 0.5),
                                prior.range = c(100, .5))
  cmp_PO0 = ~ InterceptPO(1, model = "linear", prec.linear = 0.01) +
    maternPO(geometry, model = matern ) +
  cov(Broadlead, model = "linear")
  likPO0 = like(geometry ~ .,
           data = PO0,
           family = "cp",
           samplers = outline,
           domain = list(geometry = meshes[[i]]))
  models0[[i]] = bru(cmp_PO0, likPO0)
}


## models with smoothed broadlead ------------------------------------------
for(i in 1:n_meshes)
{
  print(i)
  matern <- inla.spde2.pcmatern(meshes[[i]],
                                prior.sigma = c(1, 0.5),
                                prior.range = c(100, .5))
  cmp_PO0 = ~ InterceptPO(1, model = "linear", prec.linear = 0.01) +
    maternPO(geometry, model = matern ) +
    cov(Broad_smooth, model = "linear")
  likPO0 = like(geometry ~ .,
                data = PO0,
                family = "cp",
                samplers = outline,
                domain = list(geometry = meshes[[i]]))
  models0a[[i]] = bru(cmp_PO0, likPO0)
}

#' # Results
#' ## Fixed effects
#+

# fixed effects -----------------------------------------------------------


model_list = c("models0","models0a", "models1")
covariate = c("broad", "broad_smooth", "temp")
fixed = c()
for(k in 1:length(model_list))
{
  assign("mm", get(model_list[[k]]))
  aa = purrr::map(mm,function(x) x$summary.fixed )
  names = c("Int",
            "cov")
  res = data.frame(aa[[1]][,c(1,3,5)], mesh = 1, names = names)
  colnames(res) = c("mean","0.025quant", "0.975quant", "mesh", "names")

  for(i in 2:n_meshes)
  {
    if(!is.null(dim(aa[[i]])[1]))
    {
      bb = aa[[i]][,c(1,3,5)]
      res = rbind(res, as.data.frame(cbind(bb, mesh = i, names = names)))
    }
  }
  res$model = covariate[[k]]
  fixed = rbind(fixed, res)
}

p1 =fixed  %>% dplyr::filter(names=="cov") %>%
  ggplot() + geom_point(aes(x = mean, y = mesh)) +
  geom_errorbar(aes(y = mesh, xmin = `0.025quant`, xmax = `0.975quant`,
                    color = factor(mesh))) +
  facet_grid(model~.) +
  geom_vline(xintercept = beta) +
 # coord_cartesian(xlim = c(-3,3)) +
  ggtitle("estimanted covariate") +  xlab("") + ylab("")

p2 = fixed  %>% dplyr::filter(names=="Int") %>%
  ggplot() + geom_point(aes(x = mean, y = mesh)) +
  geom_errorbar(aes(y = mesh, xmin = `0.025quant`, xmax = `0.975quant`,
                    color = factor(mesh))) +
  facet_grid(model~.) +
  geom_vline(xintercept = int) +
  coord_cartesian(xlim = c(-10,-4)) +
  ggtitle("Estimated Intercept")+  xlab("") + ylab("")

p1 + p2 + plot_layout(guides = "collect") 
  


p1 =fixed  %>% dplyr::filter(names=="cov") %>%
  ggplot() + geom_point(aes(x = mean, y = mesh)) +
  geom_errorbar(aes(y = mesh, xmin = `0.025quant`, xmax = `0.975quant`,
                    color = factor(mesh))) +
  facet_grid(model~.) +
  geom_vline(xintercept = beta) +
   coord_cartesian(xlim = c(-3,3)) +
  ggtitle("estimanted covariate") +  xlab("") + ylab("")

p2 = fixed  %>% dplyr::filter(names=="Int") %>%
  ggplot() + geom_point(aes(x = mean, y = mesh)) +
  geom_errorbar(aes(y = mesh, xmin = `0.025quant`, xmax = `0.975quant`,
                    color = factor(mesh))) +
  facet_grid(model~.) +
  geom_vline(xintercept = int) +
  coord_cartesian(xlim = c(-10,-4)) +
  ggtitle("Estimated Intercept")+  xlab("") + ylab("")

p1 + p2 + plot_layout(guides = "collect") 


#' ## Parameters for the SPDE field
#+


# SPDE params-----------------------------------------------------------


model_list = c("models0","models0a", "models1")
covariate = c("broad", "broad_smooth", "temp")
spde_param = c()
for(k in 1:length(model_list))
{
  print(k)
  assign("mm", get(model_list[[k]]))

  aa = purrr::map(mm,function(x) x$summary.hyperpar )
  names = c("range", "sd")
  res = data.frame(aa[[1]][,c(1,3,5)], mesh = 1, names = names)
  colnames(res) = c("mean","0.025quant", "0.975quant", "mesh", "names")
  for(i in 2:n_meshes)
  {
    if(!is.null(dim(aa[[i]])[1]))
    {
      bb = aa[[i]][,c(1,3,5)]
      res = rbind(res, as.data.frame(cbind(bb, mesh = i, names = names)))
    }
  }
  res$model = covariate[[k]]
  spde_param = rbind(spde_param, res)
}

p1 = spde_param %>% dplyr::filter(names =="range") %>%
  ggplot() + geom_point(aes(x = mean, y = mesh)) +
  geom_errorbar(aes(y = mesh, xmin = `0.025quant`, xmax = `0.975quant`,
                    color = factor(mesh))) +
  facet_grid(model~.) +
  geom_vline(xintercept = true_range) +
  ggtitle("Estimated range")+  xlab("") + ylab("")

p2 = spde_param %>% dplyr::filter(names =="sd") %>%
  ggplot() + geom_point(aes(x = mean, y = mesh)) +
  geom_errorbar(aes(y = mesh, xmin = `0.025quant`, xmax = `0.975quant`,
                    color = factor(mesh))) +
  facet_grid(model~.) +
  geom_vline(xintercept = true_sd) +
  coord_cartesian(xlim = c(0,6)) +
  ggtitle("Estimated sd")+  xlab("") + ylab("")

p1 + p2 + plot_layout(guides = "collect")

#' #' ## Predicted fields
#' #+

# estimated field ---------------------------------------------------------

#' 
#' mesh <- fm_mesh_2d_inla(boundary = outline,
#'                         max.edge = c(40,150),
#'                         cutoff = 5,
#'                         crs= crs(outline))
#' pix =  fm_pixels(mesh, mask = outline)
#' preds = list()
#' for(k in 1:length(model_list))
#' {
#'   print(k)
#'   assign("mm", get(model_list[[k]]))
#'   preds[[k]] = list()
#'   for(i in 1:n_meshes)
#'   {
#'     preds[[k]][[i]] <- predict(
#'       mm[[i]], pix,
#'       ~ data.frame(
#'         loglambda = (maternPO +
#'                        InterceptPO +
#'                         cov
#'                      ),
#'         spde = maternPO,
#'         cov = cov
#'         ))
#'   }
#' }
#' 
#' #' ### Posterior mean and sd of log-lambda
#' #+
#' pmean = pix
#' for(k in 1:length(model_list))
#'   for(i in c(1,n_meshes))
#'   {
#'     colname= paste(model_list[k],".mesh",i,sep = "")
#'     pmean = pmean %>% add_column(!! (colname) := preds[[k]][[i]]$loglambda$mean)
#'   }
#' #pmean = pmean %>% pivot_longer(-geometry) %>%
#' #  separate(name, c("model", "mesh"))
#' 
#' #pmean = pmean %>%  dplyr::filter(!(model=="models1" &
#' #                                     mesh=="mesh5"))
#' 
#'   ggplot(data = pmean) + geom_sf(aes(color  = models0.mesh1)) + 
#'   theme_maps + ggtitle("models0.mesh1") + 
#'   ggplot(data = pmean) + geom_sf(aes(color  = models0.mesh4)) + 
#'   theme_maps + ggtitle("models0.mesh4") + 
#'   ggplot(data = pmean) + geom_sf(aes(color  = models1.mesh1)) + 
#'     theme_maps + ggtitle("models1.mesh1") + 
#'     ggplot(data = pmean) + geom_sf(aes(color  = models1.mesh4)) + 
#'     theme_maps + ggtitle("models1.mesh4") 
#'     
#' 
#' 
#' 
#' p1 = pmean %>%
#'   ggplot() + geom_sf(aes(color  = models0.mesh1)) + 
#'   theme_maps + ggtitle("Posterior mean of log-lambda")
#' 
#' 
#' ptrue_temp = ggplot() +
#'   geom_spatraster(data = log_lambda1)  +
#'   theme_maps +
#'   ggtitle("true lambda temp")
#' 
#' ptrue_broad = ggplot() +
#'   geom_spatraster(data =log_lambda0)  +
#'   theme_maps +
#'   ggtitle("true lambda broad")
#' p2 = ptrue_broad / ptrue_temp
#' p1 + p2
#' 
#' 
#' 
#' psd = pix
#' for(k in 1:length(model_list))
#'   for(i in c(1,n_meshes))
#'   {
#'     colname= paste(model_list[k],".mesh",i,sep = "")
#'     psd = psd %>% add_column(!! (colname) := preds[[k]][[i]]$loglambda$sd)
#'   }
#' p1 = psd %>% pivot_longer(-geometry) %>%
#'   separate(name, c("model", "mesh")) %>%
#'   ggplot() + geom_sf(aes(color  = value)) + facet_grid(mesh~model) +
#'   theme_maps + ggtitle("sd of log lambda")
#' p1
#' 
#' #' ### Posterior mean and sd of matern field
#' #+
#' 
#' pmean = pix
#' for(k in 1:length(model_list))
#'   for(i in c(1,n_meshes))
#'   {
#'     colname= paste(model_list[k],".mesh",i,sep = "")
#'     pmean = pmean %>% add_column(!! (colname) := preds[[k]][[i]]$spde$mean)
#'   }
#' pmean = pmean %>% pivot_longer(-geometry) %>%
#'   separate(name, c("model", "mesh"))
#' p1 = pmean %>%
#'   ggplot() + geom_sf(aes(color  = value)) + facet_grid(mesh~model) +
#'   theme_maps + ggtitle("Posterior mean of matern field")
#' 
#' if(0)
#' {
#'   pmean %>% dplyr::filter(mesh=="mesh5") %>%
#'     ggplot() + geom_sf(aes(color  = value)) + facet_grid(mesh~model) +
#'     theme_maps + ggtitle("Posterior mean of matern field")
#' 
#' }
#' 
#' ptrue_temp = ggplot() +
#'   geom_spatraster(data = true_GF)  +
#'   theme_maps +
#'   geom_sf(data = outline, alpha = 0) +
#'   ggtitle("true matern field")
#' 
#' 
#' p1 + ptrue_temp
#' 
#' 
#' 
#' psd = pix
#' for(k in 1:length(model_list))
#'   for(i in c(1,n_meshes))
#'   {
#'     colname= paste(model_list[k],".mesh",i,sep = "")
#'     psd = psd %>% add_column(!! (colname) := preds[[k]][[i]]$spde$sd)
#'   }
#' p1 = psd %>% pivot_longer(-geometry) %>%
#'   separate(name, c("model", "mesh")) %>%
#'   ggplot() + geom_sf(aes(color  = value)) + facet_grid(mesh~model) +
#'   theme_maps + ggtitle("sd of matern field")
#' p1
#' 
#' #' # Different meshes and covariate resolution
#' #+ fig.width=10, fig.height=10
#' # plot meshes and covariates via meshes -----------------------------------
#' 
#' for(i in 1:4)
#' {
#'   mesh = meshes[[i]]
#'   mesh_loc = st_as_sf(data.frame(x = mesh$loc[,1],
#'                                  y = mesh$loc[,2]),
#'                       coords = c("x","y"))
#'   st_crs(mesh_loc) = st_crs(outline)
#' 
#'   dual = st_as_sf(book.mesh.dual(mesh))
#'   st_crs(dual)= st_crs(outline)
#' 
#'   p0 =  ggplot() + gg(mesh) + coord_equal() + theme_maps
#' 
#'   val = eval_spatial(Broadlead, where = mesh_loc)
#'   val = bru_fill_missing(Broadlead, where = mesh_loc, values = val)
#'   p1 = dual %>% mutate(val = val) %>%
#'     ggplot() + geom_sf(aes(fill = val), color=alpha("black",0)) +
#'    theme_maps + ggtitle("Broadlead")+
#'     geom_sf(data = outline, alpha = 0)
#' 
#'   val = eval_spatial(Broad_smooth, where = mesh_loc)
#'   val = bru_fill_missing(Broad_smooth, where = mesh_loc, values = val)
#'   p2 = dual %>% mutate(val = val) %>%
#'     ggplot() + geom_sf(aes(fill = val), color=alpha("black",0)) +
#'     theme_maps + ggtitle("Broad smooth") +
#'     geom_sf(data = outline, alpha = 0)
#' 
#'   val = eval_spatial(Temp, where = mesh_loc)
#'   val = bru_fill_missing(Temp, where = mesh_loc, values = val)
#'   p3 = dual %>% mutate(val = val) %>%
#'     ggplot() + geom_sf(aes(fill = val), color=alpha("black",0)) +
#'     theme_maps + ggtitle("temp")+
#'     geom_sf(data = outline, alpha = 0)
#' 
#'   plot = p0 + p1 + p2 + p3 + plot_layout(nrow = 1)
#'   print(plot)
#' }
