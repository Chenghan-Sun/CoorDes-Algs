# Ensemble Methods of Coordinate Descent Algorithms

## Numerical Experiments Detail

### Part I: Setup and Implementation <a class="anchor" id="sub1"></a>
We consider solving synthetic instances of the linear regression model with least-squares objective function: 

$$
f^{*}:=\min _{\beta \in \mathbb{R}^{p}} f(\beta)=\|y-X \beta\|_{2}^{2}
$$

using a baseline method Gradient Descent(GD), Randomized Coordinate Descent (RCD) Method (both were implemented codebase/RCD.R), and Accelerated Randomized Coordinate Descent (ARCD \ref{alg.acc}) Method (implemented in the submitted .R code named "Accelerated\_RCD.R"). The mechanism for generating the data(y, X) are described in the supplementary materials; functions for generating the simulation data could be found in the submitted .R code named "Experiment.R". 

We also consider solving the linear regression model with $L_1$ penalty: 
$$
f^{*}:=\min _{\beta \in \mathbb{R}^{p}} f(\beta)=\|y-X \beta\|_{2}^{2} + \lambda \| \beta\|_1
$$
using the Separable Coordinate Descent (SpCD \ref{alg.sep}) Method (implemented in the submitted .R code named "Separable\_RCD.R"). 

### Part I: Results and Analysis <a class="anchor" id="sub1"></a>

Figure \ref{fig.gap} shows the optimality gap versus numbers of iteration for solving different instances of linear regression with different conditinoal numbers of the matrix $X^TX$ using RCD, ARCD and GD algorithms. In each plot, the vertical axis is the objective value optimality gap $f(\beta^k) - f^*$ in log scale. The horizontal axis is the numbers of iteration. Each column corresponds to an instance with the prescribed condition number $\kappa$ of $X^TX$. 
