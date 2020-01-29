---
title: 'Gconcord: A Python and an R package for graphical Concord'
tags:
  - Python
  - R
  - Sparse graphical model
  - Sparse estimation for the precision matrix
authors:
  - name: Sang-Yun Oh
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Zhipu Zhou
    affiliation: 2
affiliations:
 - name: Lawrence Berkeley National Laboratory
   index: 1
 - name: University of California - Santa Barbara
   index: 2
date: 26 January 2020
bibliography: paper.bib
---


# Summary

Graphs have played an important role in the modern statistical analysis due to their ability to intuitively capture and describe relationships among a group of variables. The sparse precision matrix is a popular tool to characterize the underlying dependency relationships of variables. There have existed several methods for the sparse precision matrix estimation. Graphical lasso is formulated for multivariate Gaussian distributed data when observations were limited [@Friedman2008]. Its estimator can be obtained by the R package ``glasso`` or the Python module ``graphical_lasso`` in scikit-learn. To extend the application to non-Gaussian data, the SPACE (Sparse PArtial Correlation Estimation) method is proposed [@Peng2009] and can be obtained by the R package ``space``.  However, the convergence is not guaranteed in SPACE method. As an improvement, graphical Concord (CONvex CORrelation selection methoD) is proposed as a sparse penalized maximum pseudo-likelihood estimator for the precision matrix [@Khare2014, @Ali2017]. Multiple approaches have been introduced to solve the Concord optimization problem [@Oh2014, @Ali2017], but no package implements and unifies all the approaches. The purpose of this software is to fill this gap.

We develop the software packages that use Concord to provide a sparse estimation for a precision matrix. They unify methods that solve the Concord optimization problem, including coordinate-wise descent algorithm, ISTA (Iterative Soft-thresholding Algorithm) method, and FISTA (Fast Iterative Soft-thresholding Algorithm) method. We provide both the R package ``gconcord`` and the Python package ``gconcord.graphical_concord_``, both of which share the same core file implemented using C++ languages. The packages are flexible, user-friendly, and easy to maintain. The structure of the packages can be visualized in Figure 1.


![structure](https://github.com/HYzzp1733/gconcord/blob/master/JOSS%20paper/Struct.JPG)

The R package contains the Concord solver function ``gconcord``, the cross-validation function ``cv.gconcord``, and other supporting functions. The Python package contains the Concord solver class ``GraphicalConcord`` and the cross-validation function ``GraphicalConcordCV``. 


# Statement of need

The software provides a sparse estimation for the precision matrix by solving pseudo-likelihood based objective function with a combination of ![l-1](http://latex.codecogs.com/gif.latex?l_1) and Frobenius norm penalties. The objective function does not depend on the Gaussian functional form and the convergence of the algorithm is guaranteed.


# The optimization problem of the package

Suppose that ![q1](http://latex.codecogs.com/gif.latex?y_1%2C%20y_2%2C%20%5Ccdots%2C%20y_n) are i.i.d. (independent and identically distributed) samples taking values in ![](http://latex.codecogs.com/gif.latex?%5Cmathbb%7BR%7D%5Ep%20%28p%5Cgeq%202%29) from some distribution with unknown true covariance matrix ![](http://latex.codecogs.com/gif.latex?%5Cmathbf%7B%5CSigma%7D). The sample covariance matrix is computed as

![eq1](http://latex.codecogs.com/gif.latex?%5Cmathbf%7BS%7D%20%3D%20%5Cfrac%7B1%7D%7Bn-1%7D%5Csum_%7Bi%20%3D%201%7D%5En%28y_i%20-%20%5Cbar%7By%7D%29%28y_i%20-%20%5Cbar%7By%7D%29%5E%7B%5Ctop%7D%5Cquad%5Ctext%7Bwhere%7D%5Cquad%20%5Cbar%7By%7D%20%3D%20%5Cfrac%7B1%7D%7Bn%7D%5Csum_%7Bi%20%3D%201%7D%5Eny_i)

Let ![](http://latex.codecogs.com/gif.latex?%5Cmathbf%7B%5COmega%7D_D) be a diagonal matrix whose diagonal elements are the diagonal elements of the matrix ![](http://latex.codecogs.com/gif.latex?%5Cmathbf%7B%5COmega%7D). Let ![](http://latex.codecogs.com/gif.latex?%5Cmathbf%7B%5COmega%7D_%7B%5Cbackslash%20D%7D%20%3D%20%5Cmathbf%7B%5COmega%7D%20-%20%5Cmathbf%7B%5COmega%7D_D), and ![](http://latex.codecogs.com/gif.latex?%5C%7C%5Cmathbf%7B%5COmega%7D%5C%7C_1) and ![](http://latex.codecogs.com/gif.latex?%5C%7C%5Cmathbf%7B%5COmega%7D%5C%7C_F) be the ![l-1](http://latex.codecogs.com/gif.latex?l_1) norm and the Frobenius norm of the matrix ![](http://latex.codecogs.com/gif.latex?%5Cmathbf%7B%5COmega%7D). The graphical concord estimator is the solution of the optimization problem

![eq2](http://latex.codecogs.com/gif.latex?%5Cmin_%7B%5Cmathbf%7B%5COmega%7D%5Cin%5Cmathbb%7BR%7D%5E%7Bp%5Ctimes%20p%7D%7D%5CBig%28%20-%5Cfrac%7B1%7D%7B2%7D%5Clog%5Cdet%28%28%5Cmathbf%7B%5COmega%7D_D%29%5E2%29%20&plus;%20%5Cfrac%7B1%7D%7B2%7D%5Ctext%7Btr%7D%28%5Cmathbf%7BS%7D%5Cmathbf%7B%5COmega%7D%5E2%29%20&plus;%20%5C%7C%5Cmathbf%7B%5CLambda%7D_1%5Ccirc%5Cmathbf%7B%5COmega%7D_%7B%5Cbackslash%20D%7D%5C%7C_1%20&plus;%20%5Clambda_2%5C%7C%5Cmathbf%7B%5COmega%7D%5C%7C_F%5E2%20%5CBig%29)

where ![](http://latex.codecogs.com/gif.latex?%5Cmathbf%7B%5CLambda%7D_1%20%3D%20%28%28%5Clambda_%7Bij%7D%29%29_%7B1%5Cleq%20i%2Cj%5Cleq%20p%7D) is a ![](http://latex.codecogs.com/gif.latex?p%5Ctimes%20p) matrix parameter for the ![l-1](http://latex.codecogs.com/gif.latex?l_1) penalty, and ![](http://latex.codecogs.com/gif.latex?%5Clambda_2%5Cgeq0) is a scalar parameter for the Frobenius norm penalty. The symbol ![](http://latex.codecogs.com/gif.latex?%5Ccirc) denotes an element-wise multiplication between two matrices. If all non-diagonal elements of ![](http://latex.codecogs.com/gif.latex?%5Cmathbf%7B%5CLambda%7D_1) are the same (i.e. ![](http://latex.codecogs.com/gif.latex?%5Clambda_%7Bij%7D%20%3D%20%5Clambda_1%5Cgeq0)) for all ![](http://latex.codecogs.com/gif.latex?1%5Cleq%20i%5Cneq%20j%5Cleq%20p), then the ![l-1](http://latex.codecogs.com/gif.latex?l_1) term is reduced to the term ![](http://latex.codecogs.com/gif.latex?%5C%7C%5Cmathbf%7B%5CLambda%7D_1%5Ccirc%5Cmathbf%7B%5COmega%7D_%7B%5Cbackslash%20D%7D%5C%7C_1%20%3D%20%5Clambda_1%5C%7C%5Cmathbf%7B%5COmega%7D_%7B%5Cbackslash%20D%7D%5C%7C_1).



# An example for R package

Suppose we hope to estimate a sparse precision matrix for the stock returns of Dow Jones Industrial Average (DJIA) component stocks from December 1, 2017 to December 31, 2017. We extract the required data, conduct the cross-validation to find out the optimal ![](http://latex.codecogs.com/gif.latex?%5Clambda_1)and ![](http://latex.codecogs.com/gif.latex?%5Clambda_2). In the cross-validation procedure, the loss values and the quantiles of loss values over the grid of ![](http://latex.codecogs.com/gif.latex?%28%5Clambda_1%2C%20%5Clambda_2%29) can be visualized in Figure 2.

```r
library(gconcord)
data = get.data(start = "2017-12-01", end = "2017-12-31", type = "return")
set.seed(1)
cv = cv.gconcord(data, K = 5)
p1 <- cvplot(cv$val.error, main = "Loss values")
p2 <- cvplot(cv$val.error.quantile, main = "Loss quantiles")
grid.arrange(p1, p2, ncol = 2)
```

![cvresults](https://github.com/HYzzp1733/gconcord/blob/master/JOSS%20paper/fig1.jpeg)


The sparse estimation of the precision matrix is computed as

```r
omega <- gconcord(data = data, lambda1 = cv$lam1.optimal, lambda2 = cv$lam2.optimal)
```

# An example for Python package

In this example, we generate multivariate normal random numbers, conduct the cross-validation to find out the optimal penalty parameters and compute the Concord estimator. The code is shown as follows.

```Python
import numpy as np
from gconcord.graphical_concord_ import GraphicalConcord, GraphicalConcordCV

p = 50    ## data dimension
mean = [0 for i in range(p)]  ## zero-mean
cov = np.diag(np.random.uniform(1, 2, p))  ## diagonal covariance matrix
x = np.random.multivariate_normal(mean, cov, 40)

cv = GraphicalConcordCV()     ## set cross-validation 
cvans = cv.fit(x)             ## fit the data
model = GraphicalConcord(lam1 = cv.lam1, lam2 = cv.lam2)  ## set Concord solver
ans = model.fit(x)            ## fit the data
ans.omega                     ## print out the estimator
```

# Acknowledgements

Not available.

# Reference

