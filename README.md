# Homogenization-FFT
Homogenization of Elasticty at Finite Strains with FFT-based methods

This repository contains the files that accompany the paper entitled

  "Efficient MATLAB implementation of FFT-based solvers for finite-strain computational homogenization".

## Description

This paper revisits the variational FFT-based methods for homogenization of finite-strain elasticity introduced by "de Geus et. al., Finite strain FFT-based non-linear solvers made simple, CMAME 318 (2017) 412-430". We improvement their implementation in various aspects. Such improvements and reasoning for these are described in the paper abstract:

We report a MATLAB-based numerical implementation of FFT-based computational homogenization with superior robustness and efficiency, as demonstrated by the application to problems in finite-strain elasticity. The microscopic boundary value problem of classical periodic homogenization is solved by Newton-Raphson iteration, combined with a conjugate-gradient solver applied to the resulting linearized system. As compared to the open-source Python code of "de Geus et al., 2017", our implementation has four major advantages: (i) it reduces the number of floating-point operations and (ii) it reduces the memory used, while (iii) simultaneously enhancing the numerical stability, (iv) it easily implements  multiple phases with distinct constitutive laws. Aspects (i) and (ii) are achieved by exploiting the symmetry properties of the Green projection operator and the local elastic tangent stiffness tensor to eliminate multiplications by zero and repeated computations of the symmetric parts. In addition, we introduce modified data structures for storing the deformation gradient field as well as the use of function handles to permit modeling multiple phases with distinct functional relations representing different constitutive laws. We specifically demonstrate the applicability by simulating composites with multiple phases governed by different hyperelastic models. Despite such generality, the computational costs have been optimized to outperform the existing, efficient open-source code of "de Geus et al., 2017", while maintaining the readability and compactness of the implementation. The provided MATLAB code provides an accurate and efficient tool for single- and two-scale homogenization problems with applications beyond the chosen example systems. Lastly, the provided code also allows easy conversion to the version that supports GPU computing, which yields massive reduction of computation time.


The following three files are minimal:
  1. main3D.m implements the entire solution procedure for solving microscopic boundary value problem, 
  2. neo_hookeam.m and saint_venant.m implement functions representing the Neo-Hookean and Saint-Venant material laws.

Besides the minimal files, several more files/scripts are provided to serve the purpose of reproducing the numerical results presented in this paper.
