## ----setup, include=FALSE, cache=FALSE----------------------------------------
library(knitr)
# set global chunk options
# opts_chunk$set(fig.path='figure/minimal-', fig.align='center', fig.show='hold')
# options(formatR.arrow=TRUE,width=90)
knitr::opts_chunk$set(dpi=100)

## ----eval=TRUE,message=FALSE,warning=FALSE------------------------------------
library("FRK")       # for carrying out FRK       
library("sp")        # for defining points/polygons
library("dplyr")     # for easy data manipulation
library("ggplot2")   # for plotting

## ----eval=TRUE,message=FALSE,warning=FALSE------------------------------------
m <- 250                                                   # Sample size
RNGversion("3.6.0"); set.seed(1)                           # Fix seed
zdf <- data.frame(x = runif(m), y= runif(m))               # Generate random locs
zdf$Y <- 3 + sin(7 * zdf$x) + cos(9 * zdf$y)               # Latent process
zdf$z <- rpois(m, lambda = exp(zdf$Y))                     # Simulate data
coordinates(zdf) = ~x+y                                    # Turn into sp object

## ----eval=TRUE,message=FALSE,warning=FALSE,results='hide'---------------------
S <- FRK(f = z ~ 1,               # Formula to FRK
         list(zdf),               # All datasets are supplied in list
         nres = 2,                # Low-rank model to reduce run-time
         response = "poisson",    # data model
         link = "log",            # link function 
         nonconvex_hull = FALSE)  # convex hull                     
pred <- predict(S)                # prediction stage

## ----eval=TRUE,message=FALSE,warning=FALSE,results='hide'---------------------
plot_list <- plot(S, pred$newdata)
plot_list <- c(plot_list, plot_spatial_or_ST(zdf, "z"))

## ----echo=FALSE, warning=FALSE,fig.align='center',fig.cap="(Left) Simulated Poisson spatial data. (Centre) Prediction of the mean process. (Right) Uncertainty quantification of predictions; specifically the width of the 90\\% posterior predictive interval.\\label{fig:example1}",fig.subcap=c("",""), fig.width = 10, fig.height = 4----
ggpubr::ggarrange(plot_list$z + labs(fill = "data"),
          plot_list$p_mu + labs(fill = "pred."),
          plot_list$interval90_mu + labs(fill = "pred.\nuncertainty"),
          nrow = 1, legend = "top")

## ----eval=TRUE,message=FALSE,warning=FALSE------------------------------------
data("Am_data")
coordinates(Am_data) = ~ Easting + Northing # convert to sp object
GZ_df <- data.frame("Easting" = 219868.09, "Northing" = 285320.8)

## ----echo=FALSE, warning=FALSE, results='hide'--------------------------------

# centre, width, and height
makeRectangle <- function(centre, w, h) {
  vertices <- rbind(c(centre[, 1] - w/2, centre[, 2] - h/2),
                    c(centre[, 1] - w/2, centre[, 2] + h/2),
                    c(centre[, 1] + w/2, centre[, 2] + h/2),
                    c(centre[, 1] + w/2, centre[, 2] - h/2),
                    c(centre[, 1] - w/2, centre[, 2] - h/2))
  Polygon(vertices)
}


## Following Paul and Cressie (2011),
## we predict over a series of concentric square blocks centred at Ground Zero 
## (GZ), as well as a series of concentric square blocks away from GZ.
  n_schemes <- 1
  n_block <- 5
construct_block_scheme <- function() {
  
  ratio <- 1.03  # width to height ratio of the blocks
  w     <- seq(43, 250, length.out = n_block)
  h     <- w / ratio
  
  ## Treat GZ as the centre, and expand relative to GZ.
  blocks <- list()
  for(i in 1:n_block) {
    blocks[[i]] <- makeRectangle(centre = GZ_df, w = w[i], h = h[i])
    blocks[[i]] <- Polygons(list(blocks[[i]]), paste0("block", i))
  }
  
  if (n_schemes == 2) {
    ## Now shift away from GZ
    centre <- GZ_df
    centre[, 1] <- GZ_df[, 1] - 153
    centre[, 2] <- GZ_df[, 2] + 125
    for(i in (n_block + 1):(2 * n_block)) {
      blocks[[i]] <- makeRectangle(centre = centre, w = w[i - n_block], h = h[i- n_block])
      blocks[[i]] <- Polygons(list(blocks[[i]]), paste0("block", i))
    }
  }
  
  ## (set the plotting order from largest to smallest)
  pred_polygons <- SpatialPolygons(blocks, (n_schemes * n_block):1)
  coordnames(pred_polygons) <- c("Easting", "Northing")
  
  pred_polygons$Scheme <- rep(as.character(n_schemes:1), each = length(pred_polygons)/n_schemes) 
  
  return(pred_polygons)
}

blocks <- construct_block_scheme()

## ----echo=FALSE, warning=FALSE, results='hide', fig.align='center',fig.cap="Americium soil data and blocking scheme.\\label{fig:Am_data}", fig.width = 13.6, fig.height = 4.5----
nasa_palette <- c("#03006d","#02008f","#0000b6","#0001ef","#0000f6","#0428f6","#0b53f7","#0f81f3",
                  "#18b1f5","#1ff0f7","#27fada","#3efaa3","#5dfc7b","#85fd4e","#aefc2a","#e9fc0d","#f6da0c","#f5a009",
                  "#f6780a","#f34a09","#f2210a","#f50008","#d90009","#a80109","#730005")

lab1 <- xlab(as.expression(bquote("Easting /" ~ 10^5 ~ "m")))
lab2 <- ylab(as.expression(bquote("Northing /" ~ 10^5 ~ "m")))

formatter <- function(x){ 
    x/10^5 
}
x_scale <- scale_x_continuous(breaks = 10^5 *c(2.197, 2.199, 2.201), labels = formatter)  
y_scale <- scale_y_continuous(breaks = 10^5 * c(2.852, 2.854, 2.856), labels = formatter)

## Basic plot to reduce code repetition
p_basic <- ggplot(data = as.data.frame(Am_data), 
                  aes(x = Easting, y = Northing)) +
  lab1 + lab2 + x_scale + y_scale + theme_bw() + coord_fixed()

## Data on the original scale
p_data <- p_basic +
  geom_point(aes(colour = Am), size = 1)  +
  geom_point(data = GZ_df, shape = 4, size = 5) +
  scale_colour_gradientn(colours = nasa_palette,
                         labels = scales::scientific, 
                         breaks = c(250000, 750000))

## Data on the log scale
p_data_log_scale <- p_basic +
  geom_point(aes(colour = log(Am)), size = 1) +
  geom_point(data = GZ_df, shape = 4, size = 5) +
  scale_colour_gradientn(colours = nasa_palette,
                         name = "Log-Americium", 
                         breaks = c(9, 11, 13))

## Blocking scheme
p_Scheme_1_2 <- p_basic +
  geom_point(size = 0.3) +
  geom_point(data = GZ_df, shape = 4, size = 5) +
  geom_polygon(data = FRK::SpatialPolygonsDataFrame_to_df(blocks), 
               aes(group = id, colour = Scheme), alpha = 0) +
  labs(colour = "Blocking Scheme")

ggpubr::ggarrange(p_data + theme(legend.text=element_text(angle = 20)) + 
            theme(text = element_text(size=17)), 
          p_data_log_scale + theme(text = element_text(size=17)), 
          p_Scheme_1_2 + theme(text = element_text(size=17)), 
          nrow = 1, align = "hv", legend = "top")

## ----eval=TRUE,message=FALSE,warning=FALSE------------------------------------
BAUs <- auto_BAUs(manifold = plane(), 
                  type = "grid",            
                  data = Am_data,           
                  nonconvex_hull = FALSE) 

## Add covariates to the BAUs
d_cutoff <- 30.48  
d_BAU <- distR(coordinates(BAUs), GZ_df)
BAUs$x1 <- as.numeric(d_BAU < d_cutoff)
BAUs$x2 <- d_BAU * BAUs$x1
BAUs$x3 <- as.numeric(d_BAU >= d_cutoff)
BAUs$x4 <- d_BAU * (BAUs$x3)

## ----eval=TRUE,message=FALSE,warning=FALSE,results='hide'---------------------
BAUs$fs     <- 1     
Am_data$std <- 1

S <- FRK(f = Am ~ -1 + x1 + x2 + x3 + x4, data = list(Am_data),
         response = "gaussian", 
         link = "log",
         BAUs = BAUs, 
         nres = 2, 
         est_error = FALSE)

## ----eval=TRUE,message=FALSE,warning=FALSE,results='hide'---------------------
pred <- predict(S, type = c("link", "mean"))
plot_list <- plot(S, pred$newdata)

## ----echo=FALSE, message=FALSE,warning=FALSE,results='hide',fig.align='center',fig.cap="Americium point predictions.\\label{fig:Am_BAU_predictions}",fig.subcap=c("",""),fig.width = 8, fig.height = 5----
plot_list <- lapply(
  plot_list, 
  function(gg) gg + lab1 + lab2 + x_scale + y_scale)

ggpubr::ggarrange(
  plot_list$p_Y + labs(fill = "Y pred.") +   
    scale_fill_gradientn(colours = nasa_palette, labels = scales::scientific),
  plot_list$RMSPE_Y + labs(fill = "RMSPE Y"),
  plot_list$p_mu + labs(fill = "mu pred.") +    
    scale_fill_gradientn(colours = nasa_palette, labels = scales::scientific),
  plot_list$RMSPE_mu + labs(fill = "RMSPE mu"), 
  align = "hv", nrow = 2, ncol =2) 

## ----eval=TRUE,message=FALSE,warning=FALSE,results='hide'---------------------
pred <- predict(S, newdata = blocks) 

## ----echo=FALSE, warning=FALSE,fig.align='center',fig.cap="Americium block level predictions.\\label{fig:Am_block_predictions_spatial}",fig.subcap=c("",""),fig.width = 9.5, fig.height = 3.7----

block_pred_df <- FRK::SpatialPolygonsDataFrame_to_df(pred$newdata)

## Change level order to reverse the order that the blocks are plotted.
block_pred_df$block <- factor(block_pred_df$id,
                              levels = paste0("block", (n_schemes * n_block):1))
p_block <- p_basic +
  geom_point(size = 0.3) + 
  geom_polygon(data = block_pred_df, aes(fill = p_mu, group = block),
               alpha = 1, colour = "black") + 
  labs(fill = "mu pred.") +   
  scale_fill_gradientn(colours = nasa_palette, labels = scales::scientific)
  

p_block_RMSPE <- p_basic +
  geom_point(size = 0.3) +
  geom_polygon(data = block_pred_df, aes(fill = RMSPE_mu, group = block),
               alpha = 1, colour = "black") +
  scale_fill_distiller(palette = "BrBG",
                       labels = scales::scientific) + 
  labs(fill = "RMSPE mu")

ggpubr::ggarrange(p_block, p_block_RMSPE, nrow = 1, align = "hv")

