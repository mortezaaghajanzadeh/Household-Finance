# Structural Estimation in Household Finance

## Description

his project aims to estimate three parameters, time discount factor, risk aversion, and cost of stock market participation in a life cycle model of portfolio choices. In this regard, Section 2 we estimate the model, Section 3 do the counterfactual analysis for having a tax on capital income for different scenarios. Additionally, in sections 4 and 5 examine model robustness and subsample heterogeneity. 

## Table of Contents

- [Estimation](#estimation)
- [Counterfactual Analysis](#counterfactual-analysis)
- [Robustness](#robustness)
- [Subsample Heterogeneity](#subsample-heterogeneity)

## Estimation
We are using the "approximate" SMM method, so first we need to generate the simulated moments by choosing some values for the three parameters. We will do this by running the following code:
    
    ```Train.m```

we will save the selected parameters in the file ```Allp.mat``` and the simulated moments in the file ```Allm.mat```. Then we will run the following code to estimate the parameters:

    ```Estimation.m```
The estimated parameters will be saved in the file ```estimated_point.mat```.


## Counterfactual Analysis
We modify the simul function to include the tax on capital income. Then we run the following code to generate the simulated moments for different scenarios:

    ```Counterfactual.m```


## Robustness
We run the following code to estimate the parameters for different scenarios that we have in the project:

    ```Robustness.m```



## Subsample Heterogeneity
We run the following code to estimate the parameters for different subsamples that we have in the project:

    ```Subsample.m```


