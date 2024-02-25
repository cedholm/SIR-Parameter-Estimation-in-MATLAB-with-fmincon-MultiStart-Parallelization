# SIR-Parameter-Estimation-in-MATLAB-with-fmincon-MultiStart-Parallelization
Base code to do Parameter Estimation for a SIR ODE model in MATLAB using fmincon MultiStart Parallelization

SIR_MultiStart.m calls SIR_Model.m to run parameter estimation analysis for an Ordinary Differential Equation (ODE) model, using fmincon to find a local minimum and MultiStart to do multiple runs of fmincon in the given parameter ranges. The code is currently set up to use parallelization for computation with MultiStart.

You can edit:

NoStartPoints - how many interactions of fmincon will be run in the parameter ranges using MultiStart

LowerBounds/UpperBounds - changes the parameter range, here for just beta

This code outputs for each run: the parameters, minimized difference function value, and the exit flags on convergence. The code also plots the "best fit".

SIR_Model.m is a base ODE model for a Susceptible Infected Recovered (SIR) system.

$$\begin{align*}
\frac{dS}{dt} &= -\beta S I\\
\frac{dI}{dt} &= \beta S I - \gamma I\\
\frac{dR}{dt} &= \gamma I\,
\end{align*}$$

We add an additional state equation to keep track of cumulative infections (CI): 

$$\begin{align*}
\frac{dCI}{dt} &= \beta S I
\end{align*}$$

The Math185 Class Scripps College Spring 2024 with Christina Edholm.pdf are class slides for Prof. Edholm's Math 285 (Methods in Modern Modeling) course which includes using the code here, class activity, and associated worksheet/homework problems in the To Do.
