# Notebooks

!!! warning
    The notebooks were the original documentation of the project. Most of the material of the notebooks has been incorporated into the Manual section and is reproduced in the Tutorial. The notebooks are linked here for backward compatibility and to cover some content that is not yet included in the Manual section (e.g. the dual tandem DMA) or the FORTRAN integration.


The project is documented in a journal manuscript. The original submission included
12 Supplementary Jupyter Notebooks. The links open the notebooks in viewer mode via 
NBViewer. The notebooks have been edited for compatibility and readability. 
A virtual machine with working copies of the originally submitted notebooks are archived
on [zenodo](https://doi.org/10.5281/zenodo.2652893). Instructions for setting up the 
machine are in the supporting information of the manuscript.

[Manuscript (Open Access)](https://www.tandfonline.com/doi/full/10.1080/02786826.2018.1530724) 

[Notebook S1. Differential Mobility Analyzer](https://nbviewer.jupyter.org/github/mdpetters/DifferentialMobilityAnalyzers.jl/blob/master/docs/Notebook%20S01.%20Differential%20Mobility%20Analyzer.ipynb) 

This notebook introduces the Differential Mobility Analyzer (DMA) and demonstrates basic 
functions. The notebook includes Figures of the schematic of the DMA, the size dependence 
of the Cunningham slip flow correction factor, particle diffusion coefficient, penetration 
efficiency through the DMA, and the fractional charging efficiency of the bipolar charger. 
It also includes examples of the normalized DMA transfer functions.

[Notebook S2. Fredholm Integral Equation](https://nbviewer.jupyter.org/github/mdpetters/DifferentialMobilityAnalyzers.jl/blob/master/docs/Notebook%20S02.%20Fredholm%20Integral%20Equation.ipynb) 

This notebook introduces the Fredholm integral equation and derives the discretized 
solution via the forward convolution matrix. The notebook demonstrates how the convolution matrix is computed for any set of transmission functions.

[Notebook S3. Size Distribution Arithmetic](https://nbviewer.jupyter.org/github/mdpetters/DifferentialMobilityAnalyzers.jl/blob/master/docs/Notebook%20S03.%20Size%20Distribution%20Arithmetic.ipynb) 

This notebook introduces the SizeDistribution type. Seven unique operations are defined and showcased: (a) Multiplication of scalar and size distribution, (b) Multiplication of vector and size distribution, (c) Multiplication of matrix and size distribution, (d) Multiplication of size distribution and size distribution, (e) Division of size distribution and size distribution (f) Dot product of scalar and size distribution, (g) Dot product of vector and size distribution, and (h) Addition of two size distributions.

[Notebook S4. Single Mobility Classification](https://nbviewer.jupyter.org/github/mdpetters/DifferentialMobilityAnalyzers.jl/blob/master/docs/Notebook%20S04.%20Single%20Mobility%20Classification.ipynb) 

This notebook demonstrates how the software can be used to find the true size distribution of monodisperse mobility selected particles. It also demonstrates how to compute the selected number, surface area, and volume concentration.

[Notebook S5. Size Distribution Inversion Using Regularization](https://nbviewer.jupyter.org/github/mdpetters/DifferentialMobilityAnalyzers.jl/blob/master/docs/Notebook%20S05.%20Size%20Distribution%20Inversion%20Using%20Regularization.ipynb) 

This notebook demonstrates how to invert a size distribution from a measured noisy response function. Application of the convolution matrix together with a Twomey inverse and the L-curve algorithm are used to invert a synthetic dataset.

[Notebook S6. Size Distribution Inversion of Ambient Data](https://nbviewer.jupyter.org/github/mdpetters/DifferentialMobilityAnalyzers.jl/blob/master/docs/Notebook%20S06.%20Size%20Distribution%20Inversion%20of%20Ambient%20Data.ipynb) 

This notebook applies the regularized inverse to a published dataset. Results are compared to the inversion output of the manufacturer supplied software and evaluates the degree of agreement between this package and the manufacturer software.

[Notebook S7. Size resolved CCN measurements](https://nbviewer.jupyter.org/github/mdpetters/DifferentialMobilityAnalyzers.jl/blob/master/docs/Notebook%20S07.%20Size%20resolved%20CCN%20measurements.ipynb) 

This notebook demonstrates calculations related to the the configuration where a single differential mobility analyzer is used together with a condensation particle counter and cloud condensation nuclei counter. It is demonstrated how to express the response function in terms of the language and how to fit the response function to infer the cloud droplet activation diameter.

[Notebook S8. Hygroscopicity Tandem DMA](https://nbviewer.jupyter.org/github/mdpetters/DifferentialMobilityAnalyzers.jl/blob/master/docs/Notebook%20S08.%20Hygroscopicity%20Tandem%20DMA.ipynb) 

This notebook demonstrates calculations related to the configuration where the first differential mobility is used as classifier. The output is conditioned in a humidifier and the resulting size distribution is measured using a second DMA in scanning mode together with a condensation particle counter (or other instrument) as detector. It is demonstrated how the language can be used to express the response function of the second differential mobility analyzer used either with or without bipolar charger.

[Notebook S9. Volatility Tandem DMA](https://nbviewer.jupyter.org/github/mdpetters/DifferentialMobilityAnalyzers.jl/blob/master/docs/Notebook%20S09.%20Volatility%20Tandem%20DMA.ipynb) 

This notebook demonstrates calculations related to the configuration where the first differential mobility is used as classifier. The output is conditioned in an evaporator or condenser and the resulting size distribution is measured using a second differential mobility analyzer in scanning mode together with a condensation particle counter (or other instrument) as detector. It is demonstrated how the language can be used to express the response function of the second differential mobility analyzer used either with or without bipolar charger.

[Notebook S10. Dimer Coagulation and Isolation](https://nbviewer.jupyter.org/github/mdpetters/DifferentialMobilityAnalyzers.jl/blob/master/docs/Notebook%20S10.%20Dimer%20Coagulation%20and%20Isolation.ipynb) 

This notebook demonstrates calculations related to the configuration where two DMAs are used to size select particles of opposite charge. The two populations are merged and allowed to coagulate. Coagulated dimers are isolated using an electrostatic filter. The dimers are charge neutralized and the size distribution is measured using a DMA operated in stepping or scanning mode. It is demonstrated how the language can be used to express the response function of the third differential mobility analyzer.

[Notebook S11. PartMC Simulations](https://nbviewer.jupyter.org/github/mdpetters/DifferentialMobilityAnalyzers.jl/blob/master/docs/Notebook%20S11.%20PartMC%20Simulations.ipynb)

This notebook presents an overview over coagulation theory. It is shown how to predict the coagulated distribution using Size Distribution Arithmetic (Notebook S3). The approach is compared to model predictions with the Particle-resolved Monte Carlo code for atmospheric aerosol simulation [(PartMC)](http://lagrange.mechse.illinois.edu/partmc/)

[Notebook S12. FORTRAN API](https://nbviewer.jupyter.org/github/mdpetters/DifferentialMobilityAnalyzers.jl/blob/master/docs/Notebook%20S12.%20FORTRAN%20API.ipynb) 

This notebook demonstrates how to construct a convolution matrix using a DMA transfer function defined in a FORTRAN routine. The notebook explains how to compile the routine to a shared library and setup ```ccall``` to pass variables to and from the FORTRAN routine.

[Link to Virtual Machine](https://doi.org/10.5281/zenodo.2652893)
