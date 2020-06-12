# Ensemble Methods of Coordinate Descent Algorithms

## Group member:  Ninghui Li, Chenghan Sun, Han Chen 

## Summary of Project 
We summarise the paper [Coordinate Decent Algorithms ](https://arxiv.org/pdf/1502.04759.pdf) implement the algorithms in the repository.    
Coordinate descent algorithms are iterative methods in which each iterate is obtained by fixing most components of the variable vector x at their values from the current iteration, and approximately minimizing the objective with respect to the remaining components.  They have been used in applications for many years, and their popularity continues to grow because of their usefulness in data analysis, machine learning, and other areas of current interest. 

## Repository Description 
### [/codebase](/codebase)
[Accelerated_RCD.R](/codebase/Accelerate_RCD.R): Accelerated Randomized Coordinate Descent (Nesterov 2012)   
[Pathwise_CD.R](/codebase/Pathwise_CD.R): Pathwise Coordinate descent for the lasso    
[RCD.R](/codebase/RCD.R):  Coordinate Descent method with randomized / cyclic rules and fixed step size for solving quadratic form objective function and Gradient Descent as baseline model.   
[Separable_RCD.R](/codebase/Separable_RCD.R):  Separable Coordinate Descent Algorithm for solving quadratic form objective function with L1 penalty(LASSO)   

### [/notebook](/notebook)
[STA243_final.pdf](/notebook/STA243_final.pdf)   Final report

### [/test](/test)
[Experiment.R](/test/Experiment.R) Experiment file

### [/fig](/fig)  
Some figs used in the report.

