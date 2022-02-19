# Getting Started

## Matrix generation framework

Let's consider a collection of matrices ``M(t) \in \mathbb{C}^{n \times n}``, ``n \in \mathbb{N}^*``. If ``M(t)`` verifies:

```math
\begin{aligned}
&\frac{dM(t)}{dt} = AM(t) - M(t)A, \\
&M(t=0) = M_0, \\
\end{aligned}
```

Then ``M(t)`` and ``M_0`` are similar. ``M(t)`` has the same eigenvalues as ``M_0``.

A simple proof is available in the the relevant [paper](citing.md).

## Matrix Generation Method

The idea may seem simple, but many parameters need to be determined to achieve our objective.

Denote a linear operator of matrix ``M`` determined by matrix ``A`` as ``\widetilde{A_A} = AM-MA``, ``\forall A \in \mathbb{C}^{n \times n}``, ``M \in \mathbb{C}^{n \times n}``, ``n \in \mathbb{N^*}``. Here ``AM`` and ``MA`` are the matrix-matrix multiplication operation of matrices ``A`` and ``M``. By solving the differential equation introduced in the previous section, we can firstly get the formula of ``M(t)`` with the exponential operator and then extend it by the __Taylor series formula__:

```math
\begin{aligned}
&M(t) = e^{\widetilde{A_A}(M_0)t}, \\
&M(t) = \sum_{k=0}^{\infty}\frac{t^k}{k!}(\widetilde{A_A})^k (M_0). \\
\end{aligned}
```

Through the loop ``M_{i+1}=M_i+\frac{1}{i!}(\widetilde{A_A})^i(M_0), i\in(0,+\infty)``, a very simple initial matrix ``M_0 \in \mathbb{C}^{n \times n}`` can be converted into a new sparse, non-trivial and non-Hermitian matrix ``M_{+\infty} \in \mathbb{C}^{n \times n}``, which has the same spectra but different eigenvectors with ``M_0``.

A good selection of matrix ``A`` can make ``\widetilde{(A_A)}^i`` tend to __``{0}``__ in limited steps. We select ``A`` as a [nilpotent matrix](https://en.wikipedia.org/wiki/Nilpotent_matrix), such that there exists an integer ``k`` which leads to ``A^i=0`` for all ``i \ge k``. Such ``k`` is called the nilpotency of ``A``.

The initial matrix ``M_0`` is selected as a sparse low triangular matrices, whose eigenvalues are exactly the same as the given spectrum.

!!! note
    For generating non-Hermitian and non-Symmetric sparse matrices, the initializations of ``M₀`` are different.


## Building Blocks

There are three important building blocks, which should be taken care of (in the other words, customized) by the users.

### Spectrum specified by Users

These are different constraints for the spectrum provided by users to generate targeting matrices, especially for the non-symmetric case.

#### Non-Hermitian matrix

For the Hermitian matrices, the spectrum can be specified as a `n`-element vector, whose entries (eigenvalues) can be any real or complex number, e.g.,

```math
spectrum = [1,2+i,3-3i,4,5+i, 6-21i, 7, 8 ]
```

#### Non-Symmetric matrix

For the non-Symmetric case, the spectrum should also be specified as a `n`-element vector with real and/or complex entries (eigenvalues). However, for the entries with
complex number, there are should also exist the conjugate pairs. For example, if the first complex entry is ``a+bi``, either its left or right neighboring entry should
be its conjugate ``a-bi``. This constraint is posed because of the fact, the eigenvalues of non-symmetric matrices with real values always come out as conjugate pairs, if they are complex. Voilà an example:

```math
spectrum = [1, 2-i, 2+i, 3-2i, 3+2i, -5, 4+9i, 4-9i]
```

### Initialization of ``M_0``

The initialization of ``M_0`` for non-Hermitian and non-symmetric cases are slightly difference, which is caused by their different ways to set the given spectrum.

#### Non-Hermitian matrix

For non-Hermitian matrix, the diagonal of ``M_0`` are set to be eigenvalues given by users. The strict upper triangular part is empty, and the strict lower triangular part can be
customized by the users. An example of ``M_0`` is given as below, which takes the example of spectrum we gave in previous subsection. In this example, the entries marked as ``\times`` can be either ``0`` or any other numbers specified by the users.

```math
M_0=\begin{bmatrix}
1 &  &  &  &  &  &  & \\
\times & 2+i &  &  &  &  &  & \\
\times & \times & 3-3i & & &  &  & \\
\times & \times & \times & 4 &  &  &  & \\
\times & \times & \times & \times & 5+i &  &  & \\
\times & \times & \times & \times & \times & 6-21i & &\\
\times & \times & \times & \times & \times & \times & 7 & \\
\times & \times & \times & \times & \times & \times & \times & 8\\
\end{bmatrix}
```

#### Non-Symmetric matrix

For the non-symmetric case, three diagonals of offsets ``(-1, 0, 1)`` are reserved for the given spectrum, thus, users can only customize ``M_0`` from its diagonal of offset ``-2``. Here is an example, which uses also the spectrum introduced in previous section.


```math
M_0=\begin{bmatrix}
1 &  &  &  &  &  &  & \\
 & \textcolor{green}{2} & \textcolor{green}{1} &  &  &  &  & \\
\times & \textcolor{green}{-1} & \textcolor{green}{2} &  & &  &  & \\
\times & \times &  & \textcolor{red}{3} & \textcolor{red}{2} &  &  & \\
\times & \times & \times & \textcolor{red}{-2} & \textcolor{red}{3} &  &  & \\
\times & \times & \times & \times &  & -5 &  &\\
\times & \times & \times & \times & \times &  & \textcolor{blue}{4} &  \textcolor{blue}{9}\\
\times & \times & \times & \times & \times & \times & \textcolor{blue}{-9} & \textcolor{blue}{4}\\
\end{bmatrix}
```

Cleary, for each conjugate pair of eigenvalues ``a+bi`` and ``a-bi``, three diagonals of ``M_0``
are filled by a small block of matrix

```math
\begin{vmatrix}
a & |b| \\
-|b| & a
\end{vmatrix},
```
whose eigenvalues are exactly ``a-|b|i`` and ``a+|b|i``.

!!! note
    For the initialization of ``M_0``, SMG2S.jl will take care of setting the given spectrum. The user only needs to customize the others entries of its lower triangular part.

### Generation of Nilpotent Matrix

A simple nilpotent matrix ``A`` can be a sparse matrix whose only one upper diagonal of indexing ``diag`` are non-zeros, which is filled with continuous ``1`` and ``0``. The maximum of continuous can be
fixed as a user-specific value ``nbOne``, this will leads to ``A^i=0`` in limited number of steps. Here are an example of nilpotent matrix with ``diag=1`` and ``nbOne=3``:

```math
M_0=\begin{bmatrix}
  &  & 0 &  &  &  &  &  &  &  &  & \\
  &  &  & 1 &  &  &  &  &  &  &  & \\
  &  &  &  & 1 &  &  & &  &  &  &  \\
  &  &  &  &  & 1 &  &  &  &  &  & \\
  &  &  &  &  &  & 0 &  &  &  &  & \\
  &  &  &  &  &  &  & 1 &  &  &  & \\
  &  &  &  &  &  &  &  & 1 &  &  & \\
  &  &  &  &  &  &  & &  & 0 &   & \\
  &  &  &  &  &  &  & &  &  & 1 & \\
  &  &  &  &  &  &  & &  &  &  & 1\\
  &  &  &  &  &  &  & &  &  &  & \\
  &  &  &  &  &  &  & &  &  &  &   
\end{bmatrix}
```



!!! note
     - Initialization of ``M_0`` and generation of different nilpotent matrices will result in sparse matrices with different sparsity patterns (see [Gallery: Sparsity Patterns](gallery.md)),
     - Therefore SMG2S.jl leaves as much as possible room to the users for the customization.
     - Meanwhile, we provides also some default routines for them, which are free to be used by the users if they don't want to customize by themselves.
