# Neumann.jl

[![Documentation (stable)][docs-stable-img]][docs-stable-url] [![Documentation (development)][docs-dev-img]][docs-dev-url] [![Build status][ci-status-img]][ci-status-url] [![Coverage][coverage-img]][coverage-url]

This package provides a single function `neumann` which applies Neumann's principle to determine the forbidden and allowed components of a response tensor of arbitrary order, subject to a finite set of point group symmetry constraints.

## Installation
The package can be installed via Julia's package prompt (entered by typing `]` at the REPL) and subsequently imported by calling:
```jl
pkg> add https://github.com/thchr/Neumann.jl
julia> using Neumann
```

## Theory

Neumann's principle states that any macroscopic response tensor must transform into itself under all microscopic (isogonal) point group symmetries $P$ of the associated physical system.<sup>[1]</sup> Equivalently, the response tensor must be invariant under the elements of the group $P$.

> <sup>[1]</sup> See e.g. the International Tables of Crystallography, Volume A (2016), [Section 3.2.2.1.1](https://onlinelibrary.wiley.com/iucr/itc/Ac/ch3o2v0001/sec3o2o2o1o1/).

We consider response tensors $A_{ij\ldots k}$ of order $N$ connecting an induced response $w_i$ to the product of $N$ perturbations $\{v_j^{(1)}, \ldots, v_k^{(N)}\}$, i.e.,

$$w_i = B_{ij\ldots k}v_j^{(1)}\cdots v_k^{(N)}.$$

Under a symmetry operation $g$, the response tensor transforms according to:

$$
gA_{ij\ldots k} = g_{il} g_{mj}^{-1} \cdots g_{nk}^{-1} A_{lm\ldots n}.
$$

If $g$ is provided in an orthonormal basis we have $g^{-1} = g^{\mathrm{T}}$, and the transformation rule simplifies to $gA_{ij\ldots k} = g_{il}g_{jm}\cdots g_{kn} A_{lm\ldots n}$.

Neumann's principle simply requires that $gA_{ij\ldots k} = A_{ij\ldots k}$ for every element $g$ in the isogonal point group $P$ of the underlying system.

## Examples

Neumann.jl exports a single function, `neumann`, which satisfies the linear relations imposed by Neumann's principle by solving for the null space of the associated equation system. The null space imposes relations between certain components of the response tensor and forbids (i.e., requires vanishing value of) other components.

As example, any odd-rank tensor vanishes completely under inversion symmetry:
```jl
julia> using Neumann
julia> inversion = [-1 0 0; 0 -1 0; 0 0 -1]
julia> N = 3 # tensor order (e.g., corresponding to second-harmonic generation)
julia> neumann(inversion, N)
1-element Vector{String}:
 "xxx = yxx = zxx = xyx = yyx = z" ⋯ 102 bytes ⋯ "yz = zyz = xzz = yzz = zzz = 0"
```

At a more advanced level, we can consider e.g. the constraints imposed by 4-fold rotation symmetry ($C_4$ in Schoenflies notation; 4 in Hermann-Mauguinn notation). For convenience, we can load the relevant generators of the group from [Crystalline.jl](https://github.com/thchr/Crystalline.jl):

```jl
julia> using Crystalline
julia> ops_C₄ = generators("4", PointGroup{3}) # generators of the group C₄ (4)
julia> neumann(ops_C₄, N)
8-element Vector{String}:
 "zxx = zyy"
 "zyx = -zxy"
 "xzx = yzy"
 "yzx = -xzy"
 "xxz = yyz"
 "yxz = -xyz"
 "zzz"
 "xxx = yxx = xyx = yyx = zzx = x" ⋯ 23 bytes ⋯ "zzy = zxz = zyz = xzz = yzz = 0"
```

In the above example, the symmetry operations are already returned in a Cartesian basis. For several symmetry settings of interest, this is not usually the case. In such cases, we suggest that the generators returned by Crystalline first be converted to a Cartesian setting. As example, we may consider the case of 3-fold rotation symmetry ($C_3$ in Schoenflies notation; 3 in Hermann-Mauguinn notation):
```jl
julia> ops_C₃ = generators("3", PointGroup{3}) # generators of the group C₃ (3)
julia> Rs = crystal(1,1,1,π/2,π/2,2π/3)        # a conventional coordinate system for hexagonal systems
julia> ops_C₃′ = cartesianize.(ops_C₃, Ref(Rs))
julia> neumann(ops_C₃′, N)
10-element Vector{String}:
 "xxx = -yyx = -yxy = -xyy"
 "yxx = xyx = xxy = -yyy"
 "zxx = zyy"
 "zyx = -zxy"
 "xzx = yzy"
 "yzx = -xzy"
 "xxz = yyz"
 "yxz = -xyz"
 "zzz"
 "zzx = zzy = zxz = zyz = xzz = yzz = 0"
```

Additional information is available in the documentation of `neumann` (accessible by typing `?neumann` at the Julia REPL).

### Kleinmann symmetry

For low-frequency harmonic generation, a response tensor may additionally exhibit [Kleinmann symmetry](https://en.wikipedia.org/wiki/Kleinman_symmetry). For e.g., second-harmonic generation, this implies that the response tensor exhibits the index permutation symmetry $A_{ijk}(\omega_3; \omega_1+\omega_2) = A_{ikj}(\omega_3; \omega_1+\omega_2)$ with $\omega_3 = \omega_1+\omega_2$ and $\omega_{1,2,3}$ assumed small relative to any intrinsic frequency scales of the material.
More generally, we may consider Kleinmann-like permutation symmetries of the form $A_{ij\ldots k} = A_{i P(j\ldots k)}$ with $P(j\ldots k)$ denoting any permutation of the indices $j\ldots k$.

To incorporate Kleinmann symmetry, the `kleinmann = true` keyword argument can be passed to `neumann`.
For instance, in $C_4$, the addition of Kleinmann symmetry reduces the number of independent components from 7 to 4:
```jl
julia> neumann(ops_C₄, N; kleinmann = true)  # C₄ symmetry + Kleinmann symmetry
5-element Vector{String}:
 "zxx = zyy"
 "xzx = yzy = xxz = yyz"
 "yzx = -xzy = yxz = -xyz"
 "zzz"
 "xxx = yxx = xyx = yyx = zyx = z" ⋯ 35 bytes ⋯ "zzy = zxz = zyz = xzz = yzz = 0"
```
While, in $C_3$, Kleinmann symmetry reduces the number of independent components from 9 to 6:
```jl
julia> neumann(ops_C₃′, N; kleinmann = true) # C₃ symmetry + Kleinmann symmetry
7-element Vector{String}:
 "xxx = -yyx = -yxy = -xyy"
 "yxx = xyx = xxy = -yyy"
 "zxx = zyy"
 "xzx = yzy = xxz = yyz"
 "yzx = -xzy = yxz = -xyz"
 "zzz"
 "zyx = zzx = zxy = zzy = zxz = zyz = xzz = yzz = 0"
```


[ci-status-img]:   https://github.com/thchr/Neumann.jl/workflows/CI/badge.svg
[ci-status-url]:   https://github.com/thchr/Neumann.jl/actions
[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://thchr.github.io/Neumann.jl/stable
[docs-dev-img]:    https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]:    https://thchr.github.io/Neumann.jl/dev
[coverage-img]:    https://codecov.io/gh/thchr/Neumann.jl/branch/master/graph/badge.svg
[coverage-url]:    https://codecov.io/gh/thchr/Neumann.jl