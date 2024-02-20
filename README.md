# reflected_FBM

This project contains Fortran 90 code to simualate discrete-time reflected fractional Brownian motion (FBM) in one, two, and three dimension. FBM is characterized by the decay exponent gamma (gamma is related to the Hurst exponent via gamma = 2-2H). Geometries inlude:
- semi-infinite one-dimensional interval
- finite one-dimensional interval
- two-dimensional square
- two-dimensional circular disk
- two-dimensional circular sector
- three-dimensional sphere
- three-dimensional spherical cone
  
The programs compute general characteristics such as the mean-square displacement. They also compute the spatial probability density of the random walkers. Several ways of implementing the reflection condition are implemented (hard and soft walls), see Phys. Rev. E 102, 032108 (2020), [https://doi.org/10.1103/PhysRevE.102.032108](https://doi.org/10.1103/PhysRevE.102.032108) for a detailed discussion. 

The increments (steps) of the discrete-time FBM process, i.e., the fractional Gaussian noise, are created using the effective Fourier-filtering technique, allowing one to simulate long trajectories up to 2^27 (about 134 million) time steps. The code uses Ooura's FFT package, [https://github.com/biotrump/OouraFFT](https://github.com/biotrump/OouraFFT), to perform the Fourier transformations. Conditional compilation using preprocessor directives is used to allow the same code to run serially (non-parallel) or in parallel using MPI.


If you use these codes, or programs building on these codes in academic work, please cite the following publications:
- Alexander H. O. Wada and Thomas Vojta, Fractional Brownian motion with a reflecting wall, Phys. Rev. E 97, 020102(R) (2018), [https://doi.org/10.1103/PhysRevE.97.020102](https://doi.org/10.1103/PhysRevE.97.020102)
- Alexander H. O. Wada, Alex Warhover, and Thomas Vojta: Non-Gaussian behavior of reflected fractional Brownian motion, J. Stat. Mech. 2019 033209 (2019), [https://doi.org/10.1088/1742-5468/ab02f1](https://doi.org/10.1088/1742-5468/ab02f1)
- Thomas Vojta, Samuel Halladay, Sarah Skinner, Skirmantas Janušonis, Tobias Guggenberger, and Ralf Metzler, Reflected fractional Brownian motion in one and higher dimensions, Phys. Rev. E 102, 032108 (2020), [https://doi.org/10.1103/PhysRevE.102.032108](https://doi.org/10.1103/PhysRevE.102.032108)  


