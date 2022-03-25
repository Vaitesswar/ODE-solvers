# ODE-solvers

Reaction-diffusion in a spherical porous catalyst can be modelled by the following ordinary differential equation (ODE).

![Image 1](https://user-images.githubusercontent.com/81757215/160081557-ee363a1a-5b98-4be2-a5d2-7c59f23659b3.JPG)

At present, this ODE can be solved exactly only for reaction orders 1 and 5. Their exact solutions are as follows.

![Image 2](https://user-images.githubusercontent.com/81757215/160081602-52c286c6-fd5d-4afc-80cd-343cd2f973c2.JPG)

For other reaction orders, this ODE has to be still solved numerically in order to obtain the concentration profile. Hence, the objective of this project is to explore some promising numerical methods and determine the numerical method which is able to provide the most accurate results for a wide range of reaction and diffusion conditions. Particularly, since exact solutions are already available for reaction orders 1 and 5, the numerical method which is able to reproduce the exact solutions most accurately will be concluded as the best method for solving the ODE. The 3 numerical methods that will be explored in this project are namely method of finite difference, shooting method and Adam-Bashforth method.

The MATLAB codes for the different numerical methods are provided as well as the final report.
