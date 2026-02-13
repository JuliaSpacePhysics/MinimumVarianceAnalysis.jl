# MinimumVarianceAnalysis

[![DOI](https://zenodo.org/badge/1046063112.svg)](https://doi.org/10.5281/zenodo.18635364)
[![version](https://juliahub.com/docs/General/MinimumVarianceAnalysis/stable/version.svg)](https://juliahub.com/ui/Packages/General/MinimumVarianceAnalysis)

A Julia package for minimum or maximum variance analysis (MVA).

- [x] Maximum Variance Analysis on Magnetic Field (MVAB)
- [x] Maximum Variance Analysis on Electric Field (MVAE)

## API Reference

```@autodocs
Modules = [MinimumVarianceAnalysis]
```


Error estimates for MVA:

```@docs; canonical=false
MinimumVarianceAnalysis.Δφij
MinimumVarianceAnalysis.B_x3_error
MinimumVarianceAnalysis.E_x1_error
```

## Validation with PySPEDAS

References: [`mva_eigen`](@ref), [test_minvar.py - PySPEDAS](https://github.com/spedas/pyspedas/blob/master/pyspedas/cotrans_tools/tests/test_minvar.py)

```@example pyspedas
using MinimumVarianceAnalysis
using PySPEDAS
using PySPEDAS.PythonCall
@py import pyspedas.cotrans_tools.tests.test_minvar: TestMinvar
@py import pyspedas.cotrans_tools.minvar_matrix_make: minvar_matrix_make

isapprox_eigenvector(v1, v2) = isapprox(v1, v2) || isapprox(v1, -v2)

pytest = TestMinvar()
pytest.setUpClass()

thb_fgs_gsm = get_data("idl_thb_fgs_gsm_mvaclipped1")
jl_mva_eigen = mva_eigen(thb_fgs_gsm)
jl_mva_mat = jl_mva_eigen.vectors
jl_mva_vals = jl_mva_eigen.values

py_mva_vals = PyArray(pytest.vals.y) |> vec
py_mva_mat = PyArray(pytest.mat.y[0])'
@assert isapprox(jl_mva_vals, py_mva_vals)
@assert all(isapprox_eigenvector.(eachcol(jl_mva_mat), eachcol(py_mva_mat)))
```

Since eigenvectors are only unique up to sign; therefore, the test checks if each Julia eigenvector is approximately equal to the corresponding Python eigenvector or its negative.

### Benchmark

```@example pyspedas
using Chairmarks
@b mva_eigen(thb_fgs_gsm), minvar_matrix_make("idl_thb_fgs_gsm_mvaclipped1")
```

Julia demonstrates a performance advantage of approximately 1000 times over Python, with significantly reduced memory allocations. 
Moreover, Julia's implementation is generalized for N-dimensional data.
