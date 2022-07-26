# Neumann.jl

This package provides a single function `extract_neumann` which applies Neumann's principle to determine the forbidden and allowed components of a response tensor of arbitrary order, subject to a finite set of point group symmetry constraints.

## Installation
The package can be installed via Julia's package prompt (entered by typing `]` at the REPL):
```jl
pkg> add https://github.com/thchr/Neumann.jl
```
The package can subsequently be imported by calling:
```jl
julia> using Neumann
```

## Theory

Neumann's principle states that any macroscopic response tensor must transform into itself under all microscopic (isogonal) point group symmetries $P$ of the associated physical system [1]. Equivalently, the response tensor must be invariant under the elements of the group $P$.

[1] See e.g. the International Tables of Crystallography, Volume A (2016), [Section 3.2.2.1.1](https://onlinelibrary.wiley.com/iucr/itc/Ac/ch3o2v0001/sec3o2o2o1o1/).

We consider response tensors $A_{ij\ldots k}$ of order $n$ connecting an induced response $w_i$ to the product of $n$ perturbations $\{v_j^{(1)}, \ldots, v_k^{(n)}\}$, such that $v_i = B_{ij\ldots k}v_j^{(1)}\cdotsv_k^{(n)}$. Under a symmetry operation $g$, the response tensor transforms to:

$$
(gA)_{ij\ldots k} = g_{il}g^{-1}_{mj}\cdots g^{-1}_{nk} A_{ij\ldots k}
$$
If $g$ is provided in an orthonormal basis we have $g^{-1} = g^{\mathrm{T}}$, and the transformation rule simplifies to $A_{ij\ldots k}' = g_{il}g_{jm}\cdots g_{kn} A_{ij\ldots k}$.

Neumann's principle simply states that $(gA)_{ij\ldots k} = A_{ij\ldots k}$ for every element $g$ in the isogonal point group of the underlying system.

## Examples

Neumann.jl exports a single function, `extract_neumann`, which satisfies the linear relations imposed by Neumann's principle by solving for the null space of the associated equation system. The null space imposes relations between certain components of the response tensor and forbids (i.e., requires vanishing value of) other components.

As example, any odd-rank tensor vanishes completely under inversion symmetry:
```jl
julia> using Neumann
julia> inversion = [-1 0 0; 0 -1 0; 0 0 -1]
julia> TD = 3 # tensor order (e.g., corresponding to second-harmonic generation)
julia> extract_neumann(inversion, TD)
1-element Vector{String}:
 "xxx = yxx = zxx = xyx = yyx = z" ⋯ 102 bytes ⋯ "yz = zyz = xzz = yzz = zzz = 0"
```

At a more advanced level, we can consider e.g. the constraints imposed by 4-fold rotation symmetry ($C_4$ in Schoenflies notation; 4 in Hermann-Mauguinn notation). For convenience, we can load the relevant generators of the group from [Crystalline.jl](https://github.com/thchr/Crystalline.jl):

```jl
julia> using Crystalline
julia> ops = generators("4", PointGroup{3}) # the generators of the group C₄ (4 in Hermann-Mauguinn notation)
julia> neumann_relations(ops, TD)
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
julia> ops = generators("3", PointGroup{3})
julia> Rs = crystal(1,1,1,π/2,π/2,2π/3)# a conventional coordinate system for a hexagonal system
julia> ops′ = cartesianize.(ops, Ref(Rs))
julia> neumann_relations(ops′, TD)
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