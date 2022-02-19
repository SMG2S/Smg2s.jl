__SMG2S.jl__ is a Julia implementation of the __Scalable Matrix Generation with Given Spectrum (SMG2S)__. SMG2S was initially implemented based MPI (Message Passing Interface) and C++, which is able to generate very large-scale non-Hermitian/Symmetric matrices on modern supercomputers. For more details about this parallel implementation of SMG2S, click [here](https://github.com/SMG2S/SMG2S). The main target of the MPI/C++ implementation is its parallel performance on supercomputers, thus a lot constraints
were artificially posed, which leave very few space to be able to be customized by the users.

Different with the MPI/C++ version, this Julia implementation of SMG2S is introduced which is try to give as much as possible room to the users, who generate the matrices with different sparsity patterns.  

## Features

- **Non-Hermitian** and **Non-Symmetric**:
- **Efficiency**:
- **Given Spectrum**:
- **Controllable Sparsity Patterns**:

## Installation

Currently, SMG2S.jl is only able to be installed from its github repository.

The command to install SMG2S.jl within Julia REPL (terminal) is

```julia
using Pkg;Pkg.add(PackageSpec(url="https://github.com/SMG2S/SMG2S.jl", rev="main"))
```

## Quick Example

```julia
using SMG2S

size = 100
diag_l = -20
diag_u = -5
function f(idx::Integer, size = size)
    return 10 * cos(((idx-1) / 2) * 2 * π / size) + 5
end

spec = zeros(ComplexF64, size)
Spectrum!(spec, f, size)

genMat = nonsym(size, diag_l, diag_u, spec)
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
