# Smg2s.jl: Sparse Matrix Generator with Given Spectrum

__Smg2s.jl__ is a Julia implementation of the __Sparse Matrix Generator with Given Spectrum (Smg2s)__. Smg2s was initially implemented based MPI (Message Passing Interface) and C++, which is able to generate very large-scale non-Hermitian/Symmetric __Sparse__ matrices in parallel on modern supercomputers. For more details about this parallel implementation of Smg2s, click [here](https://github.com/Smg2s/Smg2s).  

The idea of creating a sparse matrix generator came from the fact that the spectrum of matrix have large impacts on the convergence behaviour of iterative linear solvers, such as the [Krylov subspace method](https://en.wikipedia.org/wiki/Krylov_subspace). Generating very large sparse with given spectrum would be beneficial both for the study/research on the numerical methods and benchmarking of the parallel performance of existing iterative solvers on supercomputers.

The main consideration of the MPI/C++ implementation was its parallel efficiency on supercomputers, thus a lot constraints were artificially posed, which leave very few space to be customized by the users. This results in a very limited types of sparsity patterns of generated matrices.

However, the sparsity patterns are in fact really important for the evaluation the performance of algorithms of sparse matrices, especially their parallel efficiency. Different sparsity patterns can results in very different parallel performance which requires specific implementation and optimization.     

Therefore, different with the MPI/C++ version, Smg2s.jl is introduced which is try to give as much as possible room to the users to generate the matrices with a diverse types of sparsity patterns.

👉  [[Gallery of some generated sparsity patterns]](gallery.md)

## Features

- Support of generating of both **Non-Hermitian** and **Non-Symmetric** sparse matrix
- The generated matrices are naturally **Sparse** with **non-trivial** sparsity pattern
- **Given Spectrum**: the spectrum of generated matrix is the same as the one specified by the users
- Sparsity Patterns is **Controllable** (in some sense)
- it is able to very **efficiently** generate high dimension matrices

![](assets/example.png)

## Installation

Currently, Smg2s.jl is only able to be installed from its github repository.

The command to install Smg2s.jl within Julia REPL (terminal) is

```julia
using Pkg;Pkg.add(PackageSpec(url="https://github.com/Smg2s/Smg2s.jl", rev="main"))
```

## Quick Example

```julia
using Smg2s

#Define size of matrix to be generated
size = 100

#maximum number of continous `1` in the diagonal of generated nilpotent matrix
nbOne = 4
#offsets of lower diagonals between which the matrix will be initially filled
diag_l = -20
diag_u = -5

#a function to generate the user-provided spectrum, `idx` is the indexing of eigenvalues
function f(idx::Integer, size = size)
    return 10 * cos(((idx-1) / 2) * 2 * π / size) + 5
end

#create a empty vector to store the generated eigenvalues
#Attention that this vector should always in complex scalar type
spec = zeros(ComplexF64, size)

#Setting up the spectrum with function `f`
Spectrum!(spec, f, size)

#Generating a non-symmetric sparse matrix and store it into `genMat`
genMat = nonsym(nbOne, size, diag_l, diag_u, spec)
```

## Credits

The following people are involved in the development of Smg2s:
- [Xinzhe Wu](https://www.fz-juelich.de/SharedDocs/Personen/IAS/JSC/EN/staff/wu_x.html?nn=362224) (main development and algorithms)
- Serge G. Petiton (maths and algorithms)
- Hervé Gachlier (maths)
- ...

If this project is useful for your work please consider
* [Citing](citing.md) the relevant paper
* Leaving a star on the [GitHub repository](https://github.com/Smg2s/Smg2s.jl)

## Contributing

We always would like to re-factor the MPI/C++ version following the scheme in Smg2s.jl. However, we lack of man power, because
the current contributors can only work on that during very limited part time. If you have intention to take part in this project, feel
free to contact [us](mailto:xin.wu@fz-juelich.de).

## Contact

Any questions? Just send an email to [us](mailto:xin.wu@fz-juelich.de).
