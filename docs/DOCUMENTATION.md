## Documentation

### Overview
The main documentation of the package is described in a manuscript submitted ([preprint](https://)) for [Open Access](https://en.wikipedia.org/wiki/Open_access) publication to the journal [Aerosol Science & Technology](https://www.tandfonline.com/action/journalInformation?journalCode=uast20). In addition, a series of [Jupyter](http://jupyter.org/) notebooks provides example usage and additional detailed documentation. The Jupyter notebooks also serve as formal supplement to the manuscript. This remainder of this document is structured as follows.

1. Why Julia?
2. Managing Expectations (for Newcomers to Julia)
3. Installation & Getting Started
4. Supplementary Notebooks

### 1. Why Julia?
This software started as a MATLAB script to model the response function of a specific [instrument](https://www.tandfonline.com/doi/abs/10.1080/02786826.2016.1221050) developed and used in [my laboratory](http://www4.ncsu.edu/~mdpetter/). Over time the project gradually morphed from a script performing a single task to a more formal language that could be used to express differential mobility analyzer response functions in a unified way. Implementation of this language required formalism that MATLAB did not adequately support, including recursion over inline functions, support of the constructs map and mapreduce, full Unicode support within the actual code, and the ability to create and share code and visualizations that non-expert programmers can interact with in a meaningful way. I also wanted the entire software suite to be free software.

[Julia](https://julialang.org/) is a relatively young programming language that is a promising to bridge general and scientific computing. The language bears some resemblance to MATLAB syntax. In many ways Julia code is intuitive and ideally suited for scientific programming. Julia is a dynamically-typed programming language that can rival C or FORTRAN execution speeds. It supports all of the above mentioned requirements to design the language. In addition, Julia has an excellent package managing system that allows straightforward installation of the new packages and it's dependencies, which is a great advantage to entrain new users to the language community.

### 2. Managing Expectations (for Newcomers to Julia)
Although Julia is promising, the above advantages come at some cost.

#### (a) Just in Time Compiler
Julia uses a [Just in Time (JIT) compiler](https://en.wikipedia.org/wiki/Just-in-time_compilation). When Julia first starts and packages are loaded, there can be significant JIT lag. For example, this results in slow (to very slow) time to first plot. The reason for this is that the entire [Plots.jl](http://docs.juliaplots.org/latest/) package is compiled each time the interpreter starts and the Plots.jl package is first invoked. On the flip side, once a package is loaded, execution is generally almost as fast as achievable with the hardware. Because of this issue, the execution speed is slow the first time the cells in the Jupyter notebooks are called. Once the notebook has executed once, cells can be modified to interact with the code and fast results are obtained. When working with the [Julia IDE](http://junolab.org/), the Julia interpreter is usually loaded once and the issue is mitigated.

#### (b) Ease of Use
Although Julia is dynamically-typed and in syntax apparently similar to MATLAB, the language is more abstract. For example, to optimize code types are specified explicitly. The [formal type system](https://docs.julialang.org/en/stable/manual/types/) may require a more in-depth engagement by the user with programming fundamentals relative to other common scientific programming languages such as IDL/GDL or MATLAB/GNU Octave.

#### (c) Continuing Julia Development
The Julia language is undergoing rapid development and not yet stable. This means that some language constructs can change from version to version. The current version of the Package is written in Julia 0.6 and will not run on older versions. Julia 1.0 is likely being released in 2018 after which language deprecations are expected to become less of a concern. I will be using the Package for continuing data analysis and maintain the DifferentialMobilityAnalyzers.jl package to remain compatible with future versions of the language. Users that simply apply the package to solve a problem will only be affected should they opt to upgrade to the newest Julia version in the future.

### 3. Installation & Getting Started
To interactively use the Notebooks and/or use the software follow the steps below. Skip this step if you just want to view the Notebooks.

#### (a) Julia
Julia binaries for Window, MacOS, or Linux can be downloaded [here](https://julialang.org/downloads/).
#### (a) DifferentialMobilityAnalyzers.jl
The package  <b> DifferentialMobiliyAnalyzer </b> can be installed from the Julia REPL prompt with
```julia
julia> Pkg.add("DifferentialMobilityAnalyzers")
```
This installs the package and any missing dependencies.

#### (c) Notebooks
At Julia REPL prompt invoke the Jupyter Notebook server
```julia
julia> using IJulia
julia> notebook(detached=true)
```
This opens the Jupyter tree view in your default web browser.

#### (d) Load Notebooks
Copy the Jupyter notebooks in the <b> docs/ </b> folder to a location of your choice. Load a Notebook. The notebooks contain interactive figures in Javascript. Click on the "Not Trusted" button in the top right corner to enable the interactive graphics. Alternatively, execute all cells in notebook. Each code block can be executed independently with Shift-Enter.

### Supplementary Notebooks
The main documentation of the package is described in a submitted manuscript ([preprint](https://)). There are 12 Supplementary Notebooks. The links can be followed to view the Notebooks using the [online Jupyter NBViewer tool](https://nbviewer.jupyter.org/)

[Notebook S1. Differential Mobility Analyzer]() <br>
This notebook introduces the Differential Mobility Analyzer (DMA) and demonstrates basic functions embedded in the package <b> DifferentialMobilityAnalyzers.jl </b>. The notebook includes Figures of the schematic of the DMA, the size dependence of the Cunningham slip flow correction factor, particle diffusion coefficient, penetration efficiency through the DMA, and the fractional charging efficiency of the through the bipolar charger. It also includes examples of the normalized DMA transfer functions.

[Notebook S2. Fredholm Integral Equation]() <br>
This notebook introduces the Fredholm integral equation and derives the discretized solution via the forward convolution matrix. The notebook demonstrates how the convolution matrix is computed for any set of transmission functions.

[Notebook S3. Size Distribution Arithmetic]() <br>
This notebook introduces the SizeDistribution type. Seven unique operations are defined and showcased: (a) Multiplication of scalar and size distribution,
(b) Multiplication of vector and size distribution, (c) Multiplication of matrix and size distribution, (d) Multiplication of size distribution and size distribution, (e) Division of size distribution and size distribution (f) Dot product of scalar and size distribution, (g) Dot product of vector and size distribution, and (h) Addition of two size distributions.  

[Notebook S4. Single Mobility Classification]() <br>
This notebook demonstrates how to the software can be used to find the true size distribution of monodisperse mobility selected particles. It also demonstrates how to compute the selected number, surface area, and volume concentration.

[Notebook S5. Size Distribution Inversion Using Regularization]() <br>
This notebook demonstrates how to invert a size distribution from a measured noisy response function. Application of the convolution matrix together with a Twomey inverse and the L-curve algorithm are used to invert a synthetic dataset.

[Notebook S6. Size Distribution Inversion of Ambient Data]() <br>
This notebook applies the regularized inverse to a published dataset. Results are compared to the inversion output of the manufacturer supplied software and evaluates the degree of agreement between this package and the manufacturer software.

[Notebook S7. Size resolved CCN measurements]() <br>
This notebook demonstrates calculations related to the the configuration where a single differential mobility analyzer is used together with a condensation particle counter and cloud condensation nuclei counter. It is demonstrated how to express the response function in terms of the language and how to fit the response function to infer the cloud droplet activation diameter.

[Notebook S8. Hygroscopicity Tandem DMA]() <br>
This notebook demonstrates calculations related to the configuration where the first differential mobility is used as classifier. The output is conditioned in a humidifier and the resulting size distribution is measured using a second DMA in scanning mode together with a condensation particle counter (or other instrument) as detector. It is demonstrated how to the language can be used to express the response function of the second differential mobility analyzer used either with or without bipolar charger.

[Notebook S9. Volatility Tandem DMA]() <br>
This notebook demonstrates calculations related to the configuration where the first differential mobility is used as classifier. The output is conditioned in an evaporator or condenser and the resulting size distribution is measured using a second differential mobility analyzer in scanning mode together with a condensation particle counter (or other instrument) as detector. It is demonstrated how to the language can be used to express the response function of the second differential mobility analyzer used either with or without bipolar charger.

[Notebook S10. Dimer Coagulation and Isolation]() <br>
This notebook demonstrates calculations related to the configuration where two DMAs are used to size select particles of opposite charge. The two populations are merged and allowed to coagulate. Coagulated dimers are isolated using an electrostatic filter. The dimers are charge neutralized and the size distribution is measured using a DMA operated in stepping or scanning mode. It is demonstrated how to the language can be used to express the response function of the third differential mobility analyzer.

[Notebook S11. PartMC Simulations]() <br>
This notebook presents an overview over coagulation theory. It is shown how to predict the coagulated distribution using Size Distribution Arithmetic (Notebook S3). The approach is compared to model predictions with the Particle-resolved Monte Carlo code for atmospheric aerosol simulation [(PartMC)](http://lagrange.mechse.illinois.edu/partmc/)

[Notebook S12. FORTRAN API]() <br>
This notebook demonstrates how to construct a convolution matrix using a DMA transfer function defined in a FORTRAN routine. The notebook explains how to compile the routine to a shared library and setup the ccall to pass variables to and from the FORTRAN routine.
