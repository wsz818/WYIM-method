# WYIM method in Optimal Linear Combination of Biomarkers by Weighted Youden Index Maximization

We introduce a novel method for medical diagnosis modeling. Weighted Youden index serves as an important and flexible evaluation metric of diagnostic tests and we proposed a method to construct an optimal diagnosis score and determine the best cut-off point at the same time based on weighted Youden index maximization. To deal with the neither continuous nor smooth objective function, we adopt an iterative marginal optimization algorithm. It updates each component of the parameter vector marginally at each iteration step and guarantees that the objective function at each iteration step is monotonically increasing and hence is computationally stable yet reasonably fast. For high dimensional data, we also provide a forward selection algoritm for WYIM. 

Code for proposed WYIM method and WYIM forward selection is in file "WYIM_code.R". We present the usage of proposed method and compare their performances with other methods by two real examples: Duke Cardiac Catheterization Diagnostic Dataset (file "example_acath.R") and Diagnosis Alzheimer With Handwriting Dataset (file "example_DARWIN.R") with according data.

## WYIM_code.R
This file containing three functions:
### youden_w
#### Inputs
- 'y': A vector (length n) of the binary response.
- 'score': A vector (length n) of the diagnosis scores.
- 'w': Value of weight, and value range of w is greater than 0 and less than 1 (neither 0 nor 1 is allowed).

#### Outputs
- 'WYI': Estimation of weighted Youden index (WYI).
- 'cutoff': Estimation of cut-off points.

### WYIM
#### Inputs
- 'x': A matrix (n x p) of non-zero covariates.
- 'y': A vector (length n) of the binary response.
- 'w': Value of weight, and value range of w is greater than 0 and less than 1 (neither 0 nor 1 is allowed).
- 'beta_start': A vector (length p) of initial value for beta. The beta_start can be set as maximum likelihood estimator of logistic regression $\hat{\beta}\_{logistic} / \|\hat{\beta}\_{logistic}\|$, or just assign equal coefficients to each covariate $(1/\sqrt(p),\dots,1/\sqrt(p))$.
- 'maxiter': Number of maximum rounds of update iteration. If number of iterations exceeds maxiter, stop the iteration.
- 'tol': Value of threshold to stop the update iteration. If the dot product of the updated and previous normalized beta vectors differs from 1 by less than tol, stop the loop.​

#### Outputs
- 'bhat': A vector (length p) of estimated coefficient vector. The vector is normalized to unit length.
- 'iter': Number of iterations when the loop stop.
- 'WYI': Estimation of weighted Youden index (WYI) with estimated coefficient vector 'bhat'.
- 'cutoff': Estimation of cut-off points with estimated coefficient vector 'bhat'.

### WYIM_forward
#### Inputs
- 'x': A matrix (n x p) of non-zero covariates.
- 'y': A vector (length n) of the binary response.
- 'w': Value of weight, and value range of w is greater than 0 and less than 1 (neither 0 nor 1 is allowed).
- 'maxiter': Number of maximum rounds of update iteration. If number of iterations exceeds maxiter, stop the iteration.
- 'tol': Value of threshold to stop the update iteration. If the dot product of the updated and previous normalized beta vectors differs from 1 by less than tol, stop the loop.​

#### Outputs
- 'd': Number of covariates selected in final model.
- 'bhat': A vector (length p) of estimated coefficient vector in final model. The vector is normalized to unit length.
- 'WYI': Estimation of weighted Youden index (WYI) with estimated coefficient vector 'bhat' in final model.
- 'cutoff': Estimation of cut-off points with estimated coefficient vector 'bhat' in final model.
- 'bhat_forward': A matrix (d x p) of estimated coefficient vector during the forward selection.
- 'WYI_forward': A vector (length d) of estimated weighted Youden index (WYI) with according estimated coefficient vector 'bhat_forward' during the forward selection.
- 'cutoff_forward': A vector (length d) of estimated cut-off points with according estimated coefficient vector 'bhat_forward' during the forward selection.



## example_acath.R
This file provides code for example Duke Cardiac Catheterization Diagnostic Dataset. We apply four methods: logistic regression, single index, nonparametric minmax and WYIM to male group over the age of 60 of dataset. We record computing time, linear combination coefficient vector, cutoff point and according WYI for each method. Code for conditional density curves of of linear combinations given y=0 and y=1 from different methods is also provided.

## example_DARWIN.R 
This file probides code for example Diagnosis Alzheimer With Handwriting Dataset. We apply three forward selection methods: logistic regression, nonparametric minmax and WYIM to check the performances of these methods while confronted with high dimensional data. We record final linear combination coefficient vector, cutoff point and according WYI for each method in final model. Code for conditional density curves of of linear combinations given y=0 and y=1 from different methods is also provided.


## acath_male60.csv
This data is male group over the age of 60 from original data publicly available at { https://hbiostat.org/data/}. The data contains three covariates: age (Age), cad.dur (Duration of Symptoms of Coronary Artery Disease) and choleste (Cholesterol) and a binary response sigdz.
The missing data has been dealt with random forest imputation and all continuous variables has been normalized and the duration variable is log transformed. We set sigdz=1 notes for significant coronary disease.

## DARWIN_30.csv
This data includes top 30 handwriting covariates selected with individual Youden index and a response y. The original data is pulicly available at {https://archive.ics.uci.edu/dataset/732/darwin} with description of covariates. All continuous variables in this data has been normalized and we set y=1 notes for Alzheimer's disease patients.

