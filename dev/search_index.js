var documenterSearchIndex = {"docs":
[{"location":"api/#API-Reference","page":"API Reference","title":"API Reference","text":"","category":"section"},{"location":"api/#Nilpotent-Matrix","page":"API Reference","title":"Nilpotent Matrix","text":"","category":"section"},{"location":"api/","page":"API Reference","title":"API Reference","text":"Nilp(nbOne::Ti, size::Ti) where {Ti<:Integer}","category":"page"},{"location":"api/#SMG2S.Nilp-Union{Tuple{Ti}, Tuple{Ti, Ti}} where Ti<:Integer","page":"API Reference","title":"SMG2S.Nilp","text":"Nilp(nbOne::Ti, size::Ti) where {Ti<:Integer}\n\n\n\n\n\n","category":"method"},{"location":"api/","page":"API Reference","title":"API Reference","text":"Nilp(vector::AbstractVector, size::Ti) where {Ti<:Integer}","category":"page"},{"location":"api/#SMG2S.Nilp-Union{Tuple{Ti}, Tuple{AbstractVector{T} where T, Ti}} where Ti<:Integer","page":"API Reference","title":"SMG2S.Nilp","text":"Nilp(vector::AbstractVector, size::Ti) where {Ti<:Integer}\n\n\n\n\n\n","category":"method"},{"location":"api/","page":"API Reference","title":"API Reference","text":"Nilp(matrix::SparseMatrixCSC{Tv, Ti} ,size::Ti; maxdegree::Ti=80) where{Tv <: Real, Ti <: Integer}","category":"page"},{"location":"api/#SMG2S.Nilp-Union{Tuple{Ti}, Tuple{Tv}, Tuple{SparseMatrixCSC{Tv, Ti}, Ti}} where {Tv<:Real, Ti<:Integer}","page":"API Reference","title":"SMG2S.Nilp","text":"Nilp(matrix::SparseMatrixCSC{Tv, Ti} ,size::Ti; maxdegree::Ti=80) where{Tv <: Real, Ti <: Integer}\n\n\n\n\n\n","category":"method"},{"location":"api/","page":"API Reference","title":"API Reference","text":"Nilp(vec::AbstractVector, diag::Ti, size::Ti) where{Ti <: Integer}","category":"page"},{"location":"api/#SMG2S.Nilp-Union{Tuple{Ti}, Tuple{AbstractVector{T} where T, Ti, Ti}} where Ti<:Integer","page":"API Reference","title":"SMG2S.Nilp","text":"Nilp(vec::AbstractVector, diag::Ti, size::Ti) where{Ti <: Integer}\n\n\n\n\n\n","category":"method"},{"location":"api/","page":"API Reference","title":"API Reference","text":"Nilp(nbOne::Ti, diag::Ti, size::Ti) where{Ti <: Integer}","category":"page"},{"location":"api/#SMG2S.Nilp-Union{Tuple{Ti}, Tuple{Ti, Ti, Ti}} where Ti<:Integer","page":"API Reference","title":"SMG2S.Nilp","text":"Nilp(nbOne::Ti, diag::Ti, size::Ti) where{Ti <: Integer}\n\n\n\n\n\n","category":"method"},{"location":"api/#Matrix-Initialization","page":"API Reference","title":"Matrix Initialization","text":"","category":"section"},{"location":"api/","page":"API Reference","title":"API Reference","text":"initMat!(matrix::SparseMatrixCSC{Tv, Ti}, diag_l::Ti, diag_u::Ti, size::Ti; scale::Real = 1.0, shift::Real = 0.0, sparsity::Real = 0.9) where {Tv<:Complex, Ti<:Integer}","category":"page"},{"location":"api/#SMG2S.initMat!-Union{Tuple{Ti}, Tuple{Tv}, Tuple{SparseMatrixCSC{Tv, Ti}, Ti, Ti, Ti}} where {Tv<:Complex, Ti<:Integer}","page":"API Reference","title":"SMG2S.initMat!","text":"initMat!(matrix::SparseMatrixCSC{Tv, Ti}, diag_l::Ti, diag_u::Ti, size::Ti; scale::Real = 1.0, shift::Real = 0.0, sparsity::Real = 0.9) where {Tv<:Complex, Ti<:Integer}\n\n\n\n\n\n","category":"method"},{"location":"api/","page":"API Reference","title":"API Reference","text":"initMat!(matrix::SparseMatrixCSC{Tv, Ti}, diag_l::Ti, diag_u::Ti, size::Ti; scale::Real = 1.0, shift::Real = 0.0, sparsity::Real = 0.9) where {Tv<:Real, Ti<:Integer}","category":"page"},{"location":"api/#SMG2S.initMat!-Union{Tuple{Ti}, Tuple{Tv}, Tuple{SparseMatrixCSC{Tv, Ti}, Ti, Ti, Ti}} where {Tv<:Real, Ti<:Integer}","page":"API Reference","title":"SMG2S.initMat!","text":"initMat!(matrix::SparseMatrixCSC{Tv, Ti}, diag_l::Ti, diag_u::Ti, size::Ti; scale::Real = 1.0, shift::Real = 0.0, sparsity::Real = 0.9) where {Tv<:Real, Ti<:Integer}\n\n\n\n\n\n","category":"method"},{"location":"api/#Set-Spectrum","page":"API Reference","title":"Set Spectrum","text":"","category":"section"},{"location":"api/","page":"API Reference","title":"API Reference","text":"Spectrum!(spec::AbstractVector{Tv}, f::Function, size::Ti) where {Tv<:Complex, Ti<:Integer}","category":"page"},{"location":"api/#SMG2S.Spectrum!-Union{Tuple{Ti}, Tuple{Tv}, Tuple{AbstractVector{Tv}, Function, Ti}} where {Tv<:Complex, Ti<:Integer}","page":"API Reference","title":"SMG2S.Spectrum!","text":"Spectrum!(spec::AbstractVector{Tv}, f::Function, size::Ti) where {Tv<:Complex, Ti<:Integer}\n\n\n\n\n\n","category":"method"},{"location":"api/","page":"API Reference","title":"API Reference","text":"Spectrum!(spec::AbstractVector{Tv}, vec::AbstractVector{Tv}, size::Ti) where {Tv<:Complex, Ti<:Integer}","category":"page"},{"location":"api/#SMG2S.Spectrum!-Union{Tuple{Ti}, Tuple{Tv}, Tuple{AbstractVector{Tv}, AbstractVector{Tv}, Ti}} where {Tv<:Complex, Ti<:Integer}","page":"API Reference","title":"SMG2S.Spectrum!","text":"Spectrum!(spec::AbstractVector{Tv}, vec::AbstractVector{Tv}, size::Ti) where {Tv<:Complex, Ti<:Integer}\n\n\n\n\n\n","category":"method"},{"location":"api/#Generation","page":"API Reference","title":"Generation","text":"","category":"section"},{"location":"api/#Non-Hermitian-Matrix","page":"API Reference","title":"Non-Hermitian Matrix","text":"","category":"section"},{"location":"api/","page":"API Reference","title":"API Reference","text":"nonherm(nbOne::Ti, size::Ti, diag_l::Ti, diag_u::Ti, spectrum::AbstractVector{Tv}) where {Tv<:Complex, Ti<:Integer}","category":"page"},{"location":"api/#SMG2S.nonherm-Union{Tuple{Ti}, Tuple{Tv}, Tuple{Ti, Ti, Ti, Ti, AbstractVector{Tv}}} where {Tv<:Complex, Ti<:Integer}","page":"API Reference","title":"SMG2S.nonherm","text":"nonherm(nbOne::Ti, size::Ti, diag_l::Ti, diag_u::Ti, spectrum::AbstractVector{Tv}) where {Tv<:Complex, Ti<:Integer}\n\n\n\n\n\n","category":"method"},{"location":"api/","page":"API Reference","title":"API Reference","text":"nonherm(nbOne::Ti, size::Ti, diag_l::Ti, diag_u::Ti, spectrum::AbstractVector{Tv}, Am::SparseMatrixCSC{Tv, Ti}) where {Tv<:Complex, Ti<:Integer}","category":"page"},{"location":"api/#SMG2S.nonherm-Union{Tuple{Ti}, Tuple{Tv}, Tuple{Ti, Ti, Ti, Ti, AbstractVector{Tv}, SparseMatrixCSC{Tv, Ti}}} where {Tv<:Complex, Ti<:Integer}","page":"API Reference","title":"SMG2S.nonherm","text":"nonherm(nbOne::Ti, size::Ti, diag_l::Ti, diag_u::Ti, spectrum::AbstractVector{Tv}, Am::SparseMatrixCSC{Tv, Ti}) where {Tv<:Complex, Ti<:Integer}\n\n\n\n\n\n","category":"method"},{"location":"api/","page":"API Reference","title":"API Reference","text":"nonherm(size::Ti, diag_l::Ti, diag_u::Ti, spectrum::AbstractVector{Tv}, Am::SparseMatrixCSC{Tv, Ti}, nilp::Nilpotent{Ti}) where {Tv<:Complex, Ti<:Integer}","category":"page"},{"location":"api/#SMG2S.nonherm-Union{Tuple{Ti}, Tuple{Tv}, Tuple{Ti, Ti, Ti, AbstractVector{Tv}, SparseMatrixCSC{Tv, Ti}, Nilpotent{Ti}}} where {Tv<:Complex, Ti<:Integer}","page":"API Reference","title":"SMG2S.nonherm","text":"nonherm(size::Ti, diag_l::Ti, diag_u::Ti, spectrum::AbstractVector{Tv}, Am::SparseMatrixCSC{Tv, Ti}, nilp::Nilpotent{Ti}) where {Tv<:Complex, Ti<:Integer}\n\n\n\n\n\n","category":"method"},{"location":"api/#Non-Symmetric-Matrix","page":"API Reference","title":"Non-Symmetric Matrix","text":"","category":"section"},{"location":"api/","page":"API Reference","title":"API Reference","text":"nonsym(nbOne::Ti, size::Ti, diag_l::Ti, diag_u::Ti, spectrum::AbstractVector{Tv}) where {Tv<:Complex, Ti<:Integer}","category":"page"},{"location":"api/#SMG2S.nonsym-Union{Tuple{Ti}, Tuple{Tv}, Tuple{Ti, Ti, Ti, Ti, AbstractVector{Tv}}} where {Tv<:Complex, Ti<:Integer}","page":"API Reference","title":"SMG2S.nonsym","text":"nonsym(nbOne::Ti, size::Ti, diag_l::Ti, diag_u::Ti, spectrum::AbstractVector{Tv}) where {Tv<:Complex, Ti<:Integer}\n\n\n\n\n\n","category":"method"},{"location":"api/","page":"API Reference","title":"API Reference","text":"nonsym(nbOne::Ti, size::Ti, diag_l::Ti, diag_u::Ti, spectrum::AbstractVector{Tv}, Am::SparseMatrixCSC{Tv2, Ti},) where {Tv<:Complex, Tv2<:Real, Ti<:Integer}","category":"page"},{"location":"api/#SMG2S.nonsym-Union{Tuple{Ti}, Tuple{Tv2}, Tuple{Tv}, Tuple{Ti, Ti, Ti, Ti, AbstractVector{Tv}, SparseMatrixCSC{Tv2, Ti}}} where {Tv<:Complex, Tv2<:Real, Ti<:Integer}","page":"API Reference","title":"SMG2S.nonsym","text":"nonsym(nbOne::Ti, size::Ti, diag_l::Ti, diag_u::Ti, spectrum::AbstractVector{Tv}, Am::SparseMatrixCSC{Tv2, Ti},) where {Tv<:Complex, Tv2<:Real, Ti<:Integer}\n\n\n\n\n\n","category":"method"},{"location":"api/","page":"API Reference","title":"API Reference","text":"nonsym(size::Ti, diag_l::Ti, diag_u::Ti, spectrum::AbstractVector{Tv}, Am::SparseMatrixCSC{Tv2, Ti}, nilp::Nilpotent{Ti}) where {Tv<:Complex, Tv2<:Real, Ti<:Integer}","category":"page"},{"location":"api/#SMG2S.nonsym-Union{Tuple{Ti}, Tuple{Tv2}, Tuple{Tv}, Tuple{Ti, Ti, Ti, AbstractVector{Tv}, SparseMatrixCSC{Tv2, Ti}, Nilpotent{Ti}}} where {Tv<:Complex, Tv2<:Real, Ti<:Integer}","page":"API Reference","title":"SMG2S.nonsym","text":"nonsym(size::Ti, diag_l::Ti, diag_u::Ti, spectrum::AbstractVector{Tv}, Am::SparseMatrixCSC{Tv2, Ti}, nilp::Nilpotent{Ti}) where {Tv<:Complex, Tv2<:Real, Ti<:Integer}\n\n\n\n\n\n","category":"method"},{"location":"custorm/#Customization","page":"Customization","title":"Customization","text":"","category":"section"},{"location":"custorm/#Nilpotent-Matrix","page":"Customization","title":"Nilpotent Matrix","text":"","category":"section"},{"location":"custorm/#Initialization-of-Matrix","page":"Customization","title":"Initialization of Matrix","text":"","category":"section"},{"location":"custorm/#Assembling-the-Customizations","page":"Customization","title":"Assembling the Customizations","text":"","category":"section"},{"location":"getting_started/#Getting-Started","page":"Getting Started","title":"Getting Started","text":"","category":"section"},{"location":"gallery/#Gallery","page":"Gallery: Some Sparsity Patterns","title":"Gallery","text":"","category":"section"},{"location":"gallery/","page":"Gallery: Some Sparsity Patterns","title":"Gallery: Some Sparsity Patterns","text":"size: 10000, diag_l: -7800, diag_u: -7790, nbOne: 10, nilpdiag: 1000, sp: 0.0008","category":"page"},{"location":"gallery/","page":"Gallery: Some Sparsity Patterns","title":"Gallery: Some Sparsity Patterns","text":"(Image: Pattern 1)","category":"page"},{"location":"gallery/","page":"Gallery: Some Sparsity Patterns","title":"Gallery: Some Sparsity Patterns","text":"size: 10000, diag_l: -7800, diag_u: -7790, nbOne: 10, nilpdiag: 100, sp: 0.05","category":"page"},{"location":"gallery/","page":"Gallery: Some Sparsity Patterns","title":"Gallery: Some Sparsity Patterns","text":"(Image: Pattern 2)","category":"page"},{"location":"gallery/","page":"Gallery: Some Sparsity Patterns","title":"Gallery: Some Sparsity Patterns","text":"size: 10000, diag_l: -800, diag_u: -790, nbOne: 10, nilpdiag: 100, sp: 0.05","category":"page"},{"location":"gallery/","page":"Gallery: Some Sparsity Patterns","title":"Gallery: Some Sparsity Patterns","text":"(Image: Pattern 3)","category":"page"},{"location":"gallery/","page":"Gallery: Some Sparsity Patterns","title":"Gallery: Some Sparsity Patterns","text":"size: 10000, diag_l: -800, diag_u: -790, nbOne: 10, nilpdiag: 10, sp: 0.05","category":"page"},{"location":"gallery/","page":"Gallery: Some Sparsity Patterns","title":"Gallery: Some Sparsity Patterns","text":"(Image: Pattern 4)","category":"page"},{"location":"gallery/","page":"Gallery: Some Sparsity Patterns","title":"Gallery: Some Sparsity Patterns","text":"size: 10000, diag_l: -7800, diag_u: -6000, nbOne: 10, nilpdiag: 7000, sp: 0.0008","category":"page"},{"location":"gallery/","page":"Gallery: Some Sparsity Patterns","title":"Gallery: Some Sparsity Patterns","text":"(Image: Pattern 5)","category":"page"},{"location":"gallery/","page":"Gallery: Some Sparsity Patterns","title":"Gallery: Some Sparsity Patterns","text":"size: 10000, diag_l: -7800, diag_u: -7790, nbOne: 10, nilpdiag: 1000, sp: 0.8","category":"page"},{"location":"gallery/","page":"Gallery: Some Sparsity Patterns","title":"Gallery: Some Sparsity Patterns","text":"(Image: Pattern 6)","category":"page"},{"location":"gallery/","page":"Gallery: Some Sparsity Patterns","title":"Gallery: Some Sparsity Patterns","text":"size: 10000, diag_l: -7800, diag_u: -7000, nbOne: 20, nilpdiag: 1000, sp: 0.005","category":"page"},{"location":"gallery/","page":"Gallery: Some Sparsity Patterns","title":"Gallery: Some Sparsity Patterns","text":"(Image: Pattern 7)","category":"page"},{"location":"gallery/","page":"Gallery: Some Sparsity Patterns","title":"Gallery: Some Sparsity Patterns","text":"size: 10000, diag_l: -1800, diag_u: -10, nbOne: 1000, nilpdiag: 1000, sp: 0.0005","category":"page"},{"location":"gallery/","page":"Gallery: Some Sparsity Patterns","title":"Gallery: Some Sparsity Patterns","text":"(Image: Pattern 8)","category":"page"},{"location":"gallery/","page":"Gallery: Some Sparsity Patterns","title":"Gallery: Some Sparsity Patterns","text":"size: 10000, diag_l: -9800, diag_u: -6800, nbOne: 10, nilpdiag: 7000, sp: 0.0005","category":"page"},{"location":"gallery/","page":"Gallery: Some Sparsity Patterns","title":"Gallery: Some Sparsity Patterns","text":"(Image: Pattern 9)","category":"page"},{"location":"gallery/","page":"Gallery: Some Sparsity Patterns","title":"Gallery: Some Sparsity Patterns","text":"size: 10000, diag_l: -9800, diag_u: -6800, nbOne: 20, nilpdiag: 1500, sp: 0.004","category":"page"},{"location":"gallery/","page":"Gallery: Some Sparsity Patterns","title":"Gallery: Some Sparsity Patterns","text":"(Image: Pattern 10)","category":"page"},{"location":"gallery/","page":"Gallery: Some Sparsity Patterns","title":"Gallery: Some Sparsity Patterns","text":"size: 10000, diag_l: -9800, diag_u: -6800, nbOne: 70, nilpdiag: 1500, sp: 0.004","category":"page"},{"location":"gallery/","page":"Gallery: Some Sparsity Patterns","title":"Gallery: Some Sparsity Patterns","text":"(Image: Pattern 11)","category":"page"},{"location":"gallery/","page":"Gallery: Some Sparsity Patterns","title":"Gallery: Some Sparsity Patterns","text":"size: 10000, diag_l: -9900, diag_u: -800, nbOne: 70, nilpdiag: 3500, sp: 0.0008","category":"page"},{"location":"gallery/","page":"Gallery: Some Sparsity Patterns","title":"Gallery: Some Sparsity Patterns","text":"(Image: Pattern 12)","category":"page"},{"location":"","page":"Home","title":"Home","text":"SMG2S.jl is a Julia implementation of the Scalable Matrix Generation with Given Spectrum (SMG2S). SMG2S was initially implemented based MPI (Message Passing Interface) and C++, which is able to generate very large-scale non-Hermitian/Symmetric matrices on modern supercomputers. For more details about this parallel implementation of SMG2S, click here. The main target of the MPI/C++ implementation is its parallel performance on supercomputers, thus a lot constraints were artificially posed, which leave very few space to be able to be customized by the users.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Different with the MPI/C++ version, this Julia implementation of SMG2S is introduced which is try to give as much as possible room to the users, who generate the matrices with different sparsity patterns.  ","category":"page"},{"location":"#Features","page":"Home","title":"Features","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Non-Hermitian and Non-Symmetric:\nEfficiency:\nGiven Spectrum:\nControllable Sparsity Patterns:","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Currently, SMG2S.jl is only able to be installed from its github repository.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The command to install SMG2S.jl within Julia REPL (terminal) is","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg;Pkg.add(PackageSpec(url=\"https://github.com/SMG2S/SMG2S.jl\", rev=\"main\"))","category":"page"},{"location":"#Quick-Example","page":"Home","title":"Quick Example","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"using SMG2S\n\nsize = 100\ndiag_l = -20\ndiag_u = -5\nfunction f(idx::Integer, size = size)\n    return 10 * cos(((idx-1) / 2) * 2 * π / size) + 5\nend\n\nspec = zeros(ComplexF64, size)\nSpectrum!(spec, f, size)\n\ngenMat = nonsym(size, diag_l, diag_u, spec)","category":"page"},{"location":"#Credits","page":"Home","title":"Credits","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The following people are involved in the development of SMG2S:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Xinzhe Wu (main development and algorithms)\nSerge Petiton (math and algorithms)\nHervé Gachlier (maths)\n...","category":"page"},{"location":"","page":"Home","title":"Home","text":"If this project is useful for your work please consider","category":"page"},{"location":"","page":"Home","title":"Home","text":"Citing the relevant paper\nLeaving a star on the GitHub repository","category":"page"},{"location":"#Licence","page":"Home","title":"Licence","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"SMG2S.jl is licensed under the MIT License. For more details click here.","category":"page"},{"location":"citing/#Citing-SMG2S","page":"Citing SMG2S","title":"Citing SMG2S","text":"","category":"section"},{"location":"citing/","page":"Citing SMG2S","title":"Citing SMG2S","text":"If you find SMG2S useful in your project, we kindly request that you cite the following paper:","category":"page"},{"location":"citing/","page":"Citing SMG2S","title":"Citing SMG2S","text":"@article{wu2020parallel,\n  title={A parallel generator of non-Hermitian matrices computed from given spectra},\n  author={Wu, Xinzhe and Petiton, Serge G and Lu, Yutong},\n  journal={Concurrency and Computation: Practice and Experience},\n  volume={32},\n  number={20},\n  pages={e5710},\n  year={2020},\n  publisher={Wiley Online Library},\n  doi={https://doi.org/10.1002/cpe.5710}\n}","category":"page"},{"location":"citing/","page":"Citing SMG2S","title":"Citing SMG2S","text":"A preprint can be downloaded here.","category":"page"},{"location":"citing/","page":"Citing SMG2S","title":"Citing SMG2S","text":"The very initial idea of SMG2S is coming from here.","category":"page"}]
}
