# WaveEquationBox

Program that solves the wave equation on a box, given a suitable initial condition. 

The domain chosen is [0,1].

The discretization method is piecewise linear continuous functions.

The evolution is described via finite differences.

An example solution, at different time steps:

![Wavesolution](https://user-images.githubusercontent.com/29484930/75716215-be904180-5c9c-11ea-88b9-38fb5a8eb3b1.png)

Stability of the Energy, for one particular solution, at different time steps:

![Stability](https://user-images.githubusercontent.com/29484930/75716455-352d3f00-5c9d-11ea-9b3a-30fd4e69a95e.png)

Convergence of the energy, for t=10, at different resolutions:

![Convergence](https://user-images.githubusercontent.com/29484930/75716577-7b829e00-5c9d-11ea-8145-3fa6f7d2074e.png)
