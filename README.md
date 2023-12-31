# SIR-Parameter-Estimation-in-MATLAB-with-fmincon-MultiStart-Parallelization
Base code to do Parameter Estimation for an SIR ODE model in MATLAB using fmincon MultiStart Parallelization

SIR_MultiStart.m calls SIR_Model.m to run parameter estimation analysis for an Ordinary Differential Equation (ODE) model, using fmincon to find a local minimum and MultiStart to do multiple runs of fmincon in the give parameter ranges. The code is currently set-up to use paralellization for computation with MultiStart.

You can edits:

NoStartPoints - how many interations of fmincon will be run in the parameter ranges using MultiStart

LowerBounds/UpperBounds - chnages the paramter range, here for just beta

This code outputs for each run: the parameters, minimized difference function value, and the exit flags on convergence. The code also plots the "best fit".

SIR_Model.m is a base ODE model for a Susceptible Infected Recovered (SIR) system.

$$\begin{align*}
\frac{dS}{dt} &= -\beta S I\\
\frac{dI}{dt} &= \beta S I - \gamma I\\
\frac{dR}{dt} &= \gamma I\,
\end{align*}$$

We add in an additional state equation to keep track of cumulative infections (CI): 

$$\begin{align*}
\frac{dCI}{dt} &= \beta S I
\end{align*}$$
