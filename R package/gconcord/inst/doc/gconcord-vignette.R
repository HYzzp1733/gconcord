## ---- setup, include=FALSE----------------------------------------------------
if(!require("MASS")){ stop("You need to install the package MASS")}
if(!require("knitr")){ stop("You need to install the package knitr")}
if(!require("gridExtra")){ stop("You need to install the package gridExtra")}
if(!require("kableExtra")){ stop("You need to install the package kableExtra")}
if(!require("MASS")){ stop("You need to install the package MASS")}
if(!require("dplyr")){ stop("You need to install the package dplyr")}
if(!require("knitr")){ stop("You need to install the package knitr")}
knitr::opts_chunk$set(cache = FALSE)
indent1 = '    '
indent2 = paste(rep(indent1, 2), collapse='')
indent3 = paste(rep(indent1, 3), collapse='')

## ----message = FALSE, warning = FALSE, echo = FALSE---------------------------
library(gconcord)

## ----message = FALSE----------------------------------------------------------
data = get.data( start = "2017-12-01", end = "2017-12-31", type = "return")
dim(data)

## ----message = FALSE----------------------------------------------------------
set.seed(1)
res1 = cv.gconcord(data, K = 5)
res1$lam1.optimal  ## optimal lambda1
res1$lam2.optimal  ## optimal lambda2

## -----------------------------------------------------------------------------
negL <- function(omega, data){
  n <- nrow(data)
  loss <- -n * log(max(det(omega), 1e-15)) + ## add a threshold to avoid log(0)
    sum(diag(omega %*% t(data) %*% data))
  return(loss)
}

## ----message = FALSE----------------------------------------------------------
set.seed(1)
res2 = cv.gconcord(data, K = 5, FUN = negL)
res2$lam1.optimal   ## optimal lambda1
res2$lam2.optimal   ## optimal lambda2

## ----fig.align="center", fig.width=7, fig.height=3.5, fig.cap = "Predictive risk loss function"----
p1 <- cvplot(res1$val.error, main = "Loss values")
p2 <- cvplot(res1$val.error.quantile, main = "Loss quantiles")
library(gridExtra)
grid.arrange(p1, p2, ncol = 2)

## ----fig.align="center", fig.width=7, fig.height=3.5, fig.cap = "Negative log-likelihood loss function"----
p1 <- cvplot(res2$val.error, main = "Loss values")
p2 <- cvplot(res2$val.error.quantile, main = "Loss quantiles")
grid.arrange(p1, p2, ncol = 2)

## ----fig.align="center", fig.width=7, fig.height=3.5, message = FALSE---------
omega1 <- gconcord(data = data, lambda1 = res1$lam1.optimal, 
                   lambda2 = res1$lam2.optimal)
omega2 <- gconcord(data = data, lambda1 = res2$lam1.optimal, 
                   lambda2 = res2$lam2.optimal)
p1 <- graphplot(omega1, edge.width = 0.5, varnames = colnames(data), 
                main = "Predictive risk")
p2 <- graphplot(omega2, edge.width = 0.5, varnames = colnames(data), 
                main = "Negative log-likelihood")
grid.arrange(p1, p2, ncol = 2)

## -----------------------------------------------------------------------------
# Define categories
All <- colnames(data)
Tech <- c("AAPL","CSCO","IBM","INTC","MSFT","V") ## Technology
Food <- c("HD","KO","MCD","NKE","PG","WMT")      ## Food
Pham <- c("JNJ","MRK","PFE","UNH")               ## Pharmaceuticals
Cons <- c("BA","CAT","MMM","UTX")                ## Construction
Fina <- c("AXP","GS","JPM","TRV")                ## Finance
Engy <- c("CVX","GE","XOM")                      ## Energy
Chem <- c("DWDP")                                ## Chemistry
Ettm <- c("DIS")                                 ## Entertainment
Tele <- c("VZ")                                  ## Telecommunication

## -----------------------------------------------------------------------------
# Construct prior penalty matrix Lambda
Lam <- matrix(100, ncol(data), ncol(data), dimnames = list(All, All))
Lam[Tech, Tech] <- 0.01
Lam[Food, Food] <- 0.02
Lam[Pham, Pham] <- 0.05
Lam[Cons, Cons] <- 0.01
Lam[Fina, Fina] <- 0.01
Lam[Engy, Engy] <- 0.03
Lam[Chem, Chem] <- 0.04
Lam[Ettm, Ettm] <- 0.05
Lam[Tele, Tele] <- 0.04

## ----fig.align="center", fig.width=7, fig.height=3.5, fig.cap="Predictive risk loss function"----
# Predictive risk loss function
res3 <- cv.gconcord(data = data, lam1.vec = Lam, K = 5)  
par(mfrow=c(1,2))
cvplot(res3$val.error, ylab = "Loss values")
cvplot(res3$val.error.quantile, ylab = "Loss quantiles")

## ----echo = TRUE, fig.align="center", fig.width=7, fig.height=3.5, fig.cap="Negative log-likelihood loss function"----
# Negative log-likelihood loss function
res4 <- cv.gconcord(data = data, lam1.vec = Lam, K = 5, FUN = negL)
par(mfrow=c(1,2))
cvplot(res4$val.error, ylab = "Loss values")
cvplot(res4$val.error.quantile, ylab = "Loss quantiles")

## ----echo = FALSE, fig.align="center", fig.width=7, fig.height=3.5------------
omega3 <- gconcord(data = data, lambda1 = res3$lam1.optimal, lambda2 = res3$lam2.optimal)
omega4 <- gconcord(data = data, lambda1 = res4$lam1.optimal, lambda2 = res4$lam2.optimal)
p1 <- graphplot(omega3, edge.width = 0.5, varnames = colnames(data), 
                main = "Predictive risk")
p2 <- graphplot(omega4, edge.width = 0.5, varnames = colnames(data), 
                main = "Negative log-likelihood")
grid.arrange(p1, p2, ncol = 2)

## ----eval = FALSE-------------------------------------------------------------
#  data2 = get.data( start = "2008-03-01", end = "2008-03-31", na.rm = FALSE)
#  
#  # Assume data are not available and only S is available
#  S = cov(data2, use = "complete.obs")   ## assume only S is available
#  omega <- gconcord(S = S, lambda1 = 0.2, lambda2 = 0.1)

## ----echo = FALSE, eval = FALSE-----------------------------------------------
#  # Sparse matrix generator
#  Generate.sparse.mat <- function(k, sparsity = 0.3, seed = 1){
#    library(pracma)
#    while (TRUE) {
#      # generate the symmetric sparsity mask
#      set.seed(seed)
#      mask = rand(k)
#      mask = mask * (mask < sparsity)
#      mask[lower.tri(mask, diag = TRUE)] = 0
#      mask = mask + t(mask) + eye(k)
#      mask[mask > 0] = 1
#  
#      # generate the symmetric precision matrix
#      set.seed(seed)
#      omega = matrix(rnorm(k^2), k)
#      omega[lower.tri(omega, diag = TRUE)] = 0
#      omega = omega + t(omega) + eye(k)
#  
#      # apply the reqired sparsity
#      omega = omega * mask
#      omega = omega - (min(eig(omega))-.1) * eye(k)
#  
#      if(sum(eigen(omega)$values > 0) == k) {
#        break
#      } else {
#        print('Theta is not positive definite!')
#      }
#    }
#    return(omega)
#  }

## ----message = FALSE, warning = FALSE, eval = FALSE, echo = FALSE-------------
#  library(foreach)
#  library(iterators)
#  library(doParallel)
#  library(tcltk)
#  
#  Generate.running.time <- function(p, gamma){
#  
#    # generate sample covariance matrix
#    omega <- Generate.sparse.mat(p)
#    set.seed(p)
#    data <- MASS::mvrnorm(n = gamma * p, mu = rep(0, p), Sigma = solve(omega))
#    S = cov(data)
#  
#    # parallel computing gconcord under 3 methods
#    methods = c("coordinatewise", "ista", "fista")
#    n = length(methods)
#    a <- detectCores() - 1
#    cl <- makeCluster(a)
#    registerDoParallel(cl)
#    time3 <- system.time({
#      # clusterExport(cl, c("n"))
#      k <- foreach(i = icount(n), .packages = "tcltk", .combine = rbind.data.frame) %dopar% {
#        if(!exists("pb")) pb <- tkProgressBar("Parallel task", min = 1, max = n)
#        setTkProgressBar(pb, i)
#        Sys.sleep(0.05)
#        # Functions start here
#        library(gconcord)
#        start <- Sys.time()
#        omega <- gconcord(S = S, method = methods[i], lambda1 = 0.1, lambda2 = 0.2, maxit = 300)
#        end <- Sys.time()
#        ans <- data.frame(start = start, end = end, nitr = omega$nitr)
#        rownames(ans) <- methods[i]
#        ans
#      }
#    })
#    stopCluster(cl)
#    return(k)
#  }
#  

## ----eval = FALSE, echo = FALSE-----------------------------------------------
#  address = "C:/Users/Aaron Zhou/Desktop/Github/zhipuPhdResearch/JSS/manuscript/JSS manuscript"
#  for(gamma in c(0.5, 1, 2)){
#    p <- 500
#    res <- Generate.running.time(p = p, gamma = gamma)
#    save.image(paste(address, "/p", p, "g", gamma,".RData", sep = ""))
#  }

## ----echo = FALSE, warning = FALSE, message = FALSE, eval = FALSE-------------
#  address = "C:/Users/Aaron Zhou/Desktop/Github/zhipuPhdResearch/JSS/manuscript/JSS manuscript"
#  ans <- NULL
#  
#  for(p in c(500, 800, 1000, 1500, 2000)){
#    for(gamma in c(0.5, 1, 2)){
#      load(paste(address, "/p", p, "g", gamma, ".RData", sep = ""))
#      if( p < 400 ){
#        dif <- as.vector(round( difftime(res$end, res$start, units = "secs"), 2))
#      }else{
#        dif <- as.vector(round( difftime(res$end, res$start, units = "mins"), 2))
#      }
#      nitr <- res$nitr
#      tmp <- c(p, gamma, paste(dif[1], "/", nitr[1]), paste(dif[2], "/", nitr[2]), paste(dif[3], "/", nitr[3]))
#      ans <- rbind(ans, tmp)
#    }
#  }
#  colnames(ans) <- c("p", "gamma", "cdws", "ista", "fista")
#  rownames(ans) <- NULL
#  
#  library(dplyr)
#  library(knitr)
#  library(kableExtra)

## ----echo = FALSE-------------------------------------------------------------
ans <- data.frame(p = rep(c(500, 800, 1000, 1500, 2000), each = 3),
                  gamma = rep(c(0.5, 1, 2), 5),
                  coordinate.wise = c("0.39/76", "0.35/75", "0.35/74", 
                                      "5.82/78", "6.68/82", "5.78/83",
                                      "10.35/74", "14.33/102", "14.27/100", 
                                      "52.52/92", "52.62/92", "56.25/100", 
                                      "119.04/87", "104.35/76", "119.29/87"),
                  ista = c("0.16/52", "0.14/60", "0.11/50",
                           "0.46/57", "0.63/58", "0.49/49",
                           "0.77/52", "0.77/59", "0.87/65",
                           "2.22/69", "2.32/78", "2.04/70",
                           "4.17/59", "3.98/66", "3.94/64"),
                  fista = c("0.97/175", "1.42/268", "0.71/144",
                            "3.27/152", "4.30/174", "2.46/129", 
                            "7.46/300", "5.86/248", "4.17/142", 
                            "12.06/173", "11.53/163", "13.98/244",
                            "21.49/153",  "23.84/179", "25.14/219"))

## ----echo = FALSE-------------------------------------------------------------
library(dplyr)
library(knitr)
library(kableExtra)
ans %>% 
  kable(caption = "Running time comparison. For each case of p and gamma, running time in minutes / number of iterations are shown for cdws (coordinatewise descent), ISTA, and FISTA.") %>% 
  kable_styling(bootstrap_options = "striped") 

