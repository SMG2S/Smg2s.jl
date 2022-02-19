# SMG2S
Sparse Matrix Generator with Given Spectrum

SMG2S.jl is a Julia implementation of the Sparse Matrix Generation with Given Spectrum (SMG2S).

## Features
- Both Non-Hermitian and Non-Symmetric
- Sparse and non-trivial
- Efficiency: it is able to efficiently generate a very high dimension matrices
- Given Spectrum: their spectra must be known and can be customized
- Controllable Sparsity Patterns

![Comparison of generated spectrum with given spectrum](fig/example.png)

## Installation

Currently, SMG2S.jl is only able to be installed from its github repository.

The command to install SMG2S.jl within Julia REPL (terminal) is

```julia
using Pkg;Pkg.add(PackageSpec(url="https://github.com/SMG2S/SMG2S.jl", rev="main"))
```

## Documentation
The documentation is still under development, click here for the [dev](https://smg2s.github.io/SMG2S.jl/dev/) version.


## Quick Example

```julia
using SMG2S

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

The following people are involved in the development of SMG2S:
- [Xinzhe Wu](https://github.com/brunowu) (main development and algorithms)
- [Serge Petiton](https://www.cristal.univ-lille.fr/~petiton/) (math and algorithms)
- Hervé Gachlier (maths)
- ...

If this project is useful for your work please consider
* [Citing](citing.md) the relevant paper
* Leaving a star on the [GitHub repository](https://github.com/SMG2S/SMG2S.jl)

## Contact

Any questions? Just send a email to [us](mailto:xin.wu@fz-juelich.de).


## Citing SMG2S

If you find SMG2S useful in your project, we kindly request that you cite the following paper:
```
@article{wu2020parallel,
  title={A parallel generator of non-Hermitian matrices computed from given spectra},
  author={Wu, Xinzhe and Petiton, Serge G and Lu, Yutong},
  journal={Concurrency and Computation: Practice and Experience},
  volume={32},
  number={20},
  pages={e5710},
  year={2020},
  publisher={Wiley Online Library},
  doi={https://doi.org/10.1002/cpe.5710}
}
```
A preprint can be downloaded [here](https://hal.archives-ouvertes.fr/hal-02469027/document).

The very initial idea of SMG2S is coming from [here](http://www.vecpar.org/posters/vecpar2014_submission_41.pdf).

## Licence

SMG2S.jl is licensed under the MIT License. For more details click [here](https://github.com/SMG2S/SMG2S.jl/blob/main/LICENSE).
