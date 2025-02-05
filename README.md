# CorrelatorAnalysis.jl
This is an analysis package to compute and propagate Markov Chain Monte Carlo error of correlation functions. It is written in Julia and relies on the package [ADerrors.jl](https://igit.ific.uv.es/alramos/aderrors.jl).

## Installation
To install this package [ADerrors.jl](https://igit.ific.uv.es/alramos/aderrors.jl) has to be installed first. Then, you can install it the following two ways in the julia Pkg mode.

Using HTTPS
```
pkg> add https://github.com/ResStump/CorrelatorAnalysis.jl
```

Using SSH
```
pkg> add git@github.com:ResStump/CorrelatorAnalysis.jl
```

## Documentation
Documentation is provided in the doc strings of the functions.


## Examples
The folder `examples` contains Jupyter notebooks which show the basic functionalities of `CorrelatorAnalysis`.

The notebook `correlator_analysis_random.ipynb` shows the error analysis of a randomly generated meson like correlator with autocorrelation.

The notebooks `correlator_analysis.ipynb` and `correlator_analysis_2.ipynb` repeats the same analysis for actual data of a pion correlator (two different ensembles with different number of configurations).

In `GEVP_analysis.ipynb` a sample analysis of a correlator matrix by solving the generalized eigenvalue problem (GEVP) is shown.
