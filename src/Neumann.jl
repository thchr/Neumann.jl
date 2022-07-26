module Neumann

using LinearAlgebra: checksquare, I, nullspace
using RowEchelon: rref, rref!

export neumann_relations

# ---------------------------------------------------------------------------------------- #

const CHARS = ('x','y','z','w')

# ---------------------------------------------------------------------------------------- #

"""
Build a matrix representation of the right-hand side of the constraint χᵢⱼₖ = gᵢₗgⱼₖgₗₘχₗₘₙ
(here, for third-rank tensor χ) associated with Neumann's principle.
"""
function neumann_constraint(g, 
                            #=tensor_order=#Nᵛ::Val{N},
                            #=operator_dim=#  ::Val{D}) where {N, D}

    isinteger(N) && N ≥ 1 || error("tensor order must be an integer greater than 0")
    isinteger(D)  && D ≥ 1  || error("operator dimension must be an integer greater than 0")
    minimum(size(g)) == D   || error("incompatible dimensions of `g` and `D`")
    
    g⁻¹ = inv(g)
    constraint_block = zeros(D^N, D^N)
    for (q,ijks) in enumerate(CartesianIndices(ntuple(_->1:D, Nᵛ)))
        for (p,lmns) in enumerate(CartesianIndices(ntuple(_->1:D, Nᵛ)))
            # compute gᵢₗ g⁻¹ₘⱼ g⁻¹ₙₖ … = gᵢₗ gⱼₘ gₖₙ …: we use the former variant in order to
            # be able to treat operators provided in a non-orthonormal basis (where g⁻¹≠gᵀ)
            v = g[ijks[1], lmns[1]]
            for idx = 2:N
                v *= g⁻¹[lmns[idx], ijks[idx]]
                # equivalent to `g[ijks[idx], lmns[idx]]` if g's basis is orthonormal
            end
            constraint_block[q,p] = v
        end
    end

    return constraint_block
end

# ---------------------------------------------------------------------------------------- #

function xyz_sorting(Nᵛ::Val{N}, ::Val{D}) where {N, D}
    D > 4 && error("does not support dimensions higher than 4")
    xyzs = Vector{String}(undef, D^N)
    for (q,ijks) in enumerate(CartesianIndices(ntuple(_->1:D, Nᵛ)))
        xyzs[q] = xyz_string(Tuple(ijks))
    end
    return xyzs
end

function xyz_string(ijks)
    io = IOBuffer()
    for i in ijks
        write(io, CHARS[i])
    end
    return String(take!(io))
end

# ---------------------------------------------------------------------------------------- #

"""
Return the matrix of relations used by `neumann_relations`.
"""
function neumann(ops,
                 Nᵛ::Val{N}=Val(3), Dᵛ::Val{D}=Val(3);
                 sparsify::Bool=true,
                 rref_tol::Union{Nothing,Float64}=1e-11,
                 nullspace_kws...) where {N, D}
    Id = I(D^N)
    constraint_eqs = mapfoldl(vcat, ops) do op
        Id - neumann_constraint(op, Nᵛ, Dᵛ)
    end
    A = nullspace(constraint_eqs; nullspace_kws...)

    return sparsify ? poormans_sparsification(A; rref_tol) : A
end
function neumann(op::AbstractMatrix{<:Real}, Nᵛ::Val=Val(3), Dᵛ::Val=Val(3); kws...)
    neumann(tuple(op), Nᵛ, Dᵛ; kws...)
end

function extract_relations(A::AbstractMatrix{<:Real},
                           Nᵛ::Val{N}=Val(3), Dᵛ::Val{D}=Val(3);
                           atol=1e-10) where {N, D}
    xyzs = xyz_sorting(Nᵛ, Dᵛ)

    strs = Vector{String}(undef, size(A,2)+1)
    for (i,col) in enumerate(eachcol(A))
        idxs = findall(x->abs(x)>atol, col)
        io = IOBuffer()
        for (j,idx) in enumerate(idxs)
            v = col[idx]
            v < 0 && print(io, '-')
            abs(abs(v) - 1) > atol && print(io, abs(v))
            print(io, xyzs[idx])
            j ≠ length(idxs) && print(io, " = ")
        end
        strs[i] = String(take!(io))
    end

    idxs_zerorows = findall(row->all(x->abs(x)<atol, row), eachrow(A))
    if !isempty(idxs_zerorows)
        strs[end] = join(xyzs[idxs_zerorows], " = ") * " = 0"
    else
        resize!(strs, size(A,2))
    end
    return strs
end

"""
    neumann_relations(ops, N::Integer; kws...)  -->  Vector{String}

Return the symmetry-constraints among components of a tensor imposed by the set of point
group symmetries in `ops`. The tensor order is provided via `N`. `ops` can be either a
single symmetry operation (provided as an `AbstractMatrix`) or any iterable of symmetry
operations.

The obtained relations are returned as a `Vector{String}`, with each element of this vector
giving either:

1. a free, nonzero component of the tensor (e.g., "xxx"),
2. free, nonzero - but interrelated - components of the tensor (e.g., "xxx = yyy"),
3. forbidden, zero components of the tensor (e.g., "xyz = xzy = 0").

Note that if the provided operators are not given in a Cartesian basis, the returned
relations among tensor components will match the basis of the provided operators: we
generally recommend supplying operators in a Cartesian basis for ease of interpretation.

## Keyword arguments `kws`

- `sparsify` (default, `true`): whether to attempt to sparsify relations between allowed
  nonzero components of the response tensor. A "poor man's" approximation of matrix
  sparsification is carried out using the reduced row echelon form.
- `rref_tol` (default, `1e-11`): the absolute tolerance used during row echelon reduction
  in sparsification. As the row echelon form is numerically unstable for larger matrices,
  it may be necessary to increase this tolerance for larger tensor orders. Set to `nothing`
  to use a tolerance of zero.
- `atol` (default, `1e-10`): the absolute tolerance used in assessing whether a term in a
  relation is considered vanishing or not (and similarly used to assess whether a term has
  integer coefficient). Must be smaller than `rref_tol` (if not nothing).
- `nullspace_kws` (default, empty): keyword arguments passed to `nullspace`.

## Examples

Second harmonic generation is forbidden in inversion symmetric materials
```
julia> using Neumann
julia> N = 3 # second-harmonic generation → third-rank tensor
julia> inversion = [-1 0 0; 0 -1 0; 0 0 -1]
julia> neumann_relations(inversion, N)
1-element Vector{String}:
 "xxx = yxx = zxx = xyx = yyx = z" ⋯ 102 bytes ⋯ "yz = zyz = xzz = yzz = zzz = 0"
```

But several components of the second-harmonic response are allowed under e.g. the point
group symmetry of D₃ (321):
```
julia> using Crystalline
julia> ops = generators("321", PointGroup{3}) # obtain generators of D₃ (321) from Crystalline
julia> ops .= cartesianize.(ops, Ref(crystal(1,1,1,π/2,π/2,2π/3))); # convert to a Cartesian basis
julia> neumann_relations(ops, N)
5-element Vector{String}:
 "xxx = -yyx = -yxy = -xyy"
 "zyx = -zxy"
 "yzx = -xzy"
 "yxz = -xyz"
 "yxx = zxx = xyx = xzx = zzx = x" ⋯ 41 bytes ⋯ "yyz = zyz = xzz = yzz = zzz = 0"
```
"""
function neumann_relations(ops, Nᵛ::Val{N}, Dᵛ::Val{D};
                           atol::Real=1e-10, neumann_kws...) where {N, D}
    if haskey(neumann_kws, :rref_tol)
        if atol ≤ something(neumann_kws[:rref_tol], 0.0)
            error("`atol` must exceed `rref_tol`: increase `atol` or reduce `rref_tol`")
        end
    end

    A = neumann(ops, Nᵛ, Dᵛ; neumann_kws...)
    return extract_relations(A, Nᵛ, Dᵛ; atol)
end
function neumann_relations(ops, N::Integer=3; atol::Real=1e-10, neumann_kws...)
    # determine dimension of provided operators
    D = if ops isa AbstractMatrix{<:Real} # just a single operator provided
        minimum(size(ops))
    else                                  # iterator/vector/tuple of provided operators
        D′ = minimum(size(first(ops)))
        if any(op->minimum(size(op)) != D′, ops[2:end])
            error("mismatched dimensions of provided operators")
        end
        D′
    end

    return neumann_relations(ops, Val(N), Val(D); atol, neumann_kws...)
end

# ---------------------------------------------------------------------------------------- #

"""
Poor man's "sparsification" via the reduced row echelon form.
"""
function poormans_sparsification(A; rref_tol::Union{Nothing,Float64}=1e-11)
    # following appendix E of the Qsymm paper (https://arxiv.org/abs/1806.08363)
    if !isnothing(rref_tol)
        # use a relatively low tolerance in `rref` to avoid explosions of errors
        # NB: this optional tolerance argument of `rref!` is undocumented :(
        return transpose(rref!(copy(transpose(A)), rref_tol))
    end
    return return transpose(rref(transpose(A)))
end

end # module