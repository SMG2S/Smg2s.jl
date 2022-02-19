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

The initial matrix ``M_0`` is selected as a sparse low triangular matrices, whose up triangular parts are zeros, and its diagonal is set to be the given spectrum.

!!! note
    - This is the general idea for generating non-Hermitian matrix.
    - If non-symmetric matrices are required whose entries are real number, but with possible conjugated eigenvalues, the main three diagonals are to be filled in.


## 3 Building Blocks

- Initialization of lower triangular part

- Generation of Nilpotent Matrix

- Spectrum specified by Users
