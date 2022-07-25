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
                            #=tensor_order=#TDⱽ::Val{TD}=Val(3),
                            #=operator_dim=#   ::Val{D} =Val(3)) where {TD, D}

    isinteger(TD) && TD ≥ 1 || error("tensor order must be an integer greater than 0")
    isinteger(D)  && D ≥ 1  || error("operator dimension must be an integer greater than 0")
    minimum(size(g)) == D   || error("incompatible dimensions of `g` and `D`")
    
    constraint_block = zeros(D^TD,D^TD)
    for (q,ijks) in enumerate(CartesianIndices(ntuple(_->1:D, TDⱽ)))
        for (p,lmns) in enumerate(CartesianIndices(ntuple(_->1:D, TDⱽ)))
            # compute g[i,l]*g[j,m]*g[k,n]...
            v = g[ijks[1], lmns[1]]
            for idx = 2:TD
                v *= g[ijks[idx], lmns[idx]]
            end
            constraint_block[q,p] = v
        end
    end

    return constraint_block
end

# ---------------------------------------------------------------------------------------- #

function xyz_sorting(TDᵛ::Val{TD}=Val(3), ::Val{D}=Val(3)) where {TD, D}
    D > 4 && error("does not support dimensions higher than 4")
    xyzs = Vector{String}(undef, D^TD)
    for (q,ijks) in enumerate(CartesianIndices(ntuple(_->1:D, TDᵛ)))
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
                 TDᵛ::Val{TD}=Val(3), Dᵛ::Val{D}=Val(3);
                 sparsify::Bool=true,
                 rref_tol::Union{Nothing,Float64}=1e-11,
                 nullspace_kws...) where {TD, D}
    Id = I(D^TD)
    constraint_eqs = mapfoldl(vcat, ops) do op
        Id - neumann_constraint(op, TDᵛ, Dᵛ)
    end
    A = nullspace(constraint_eqs; nullspace_kws...)

    return sparsify ? poormans_sparsification(A; rref_tol) : A
end
function neumann(op::AbstractMatrix{<:Real}, TDᵛ::Val=Val(3), Dᵛ::Val=Val(3); kws...)
    neumann(tuple(op), TDᵛ, Dᵛ; kws...)
end

function extract_relations(A::AbstractMatrix{<:Real},
                           TDᵛ::Val{TD}=Val(3), Dᵛ::Val{D}=Val(3);
                           atol=1e-10) where {TD, D}
    xyzs = xyz_sorting(TDᵛ, Dᵛ)

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
    neumann_relations(ops, ::Val{TD}, ::Val{D})

Return the free components of a tensor in Cartesian indexing, subject to the constraints of
point group symmetries in `ops`. The tensor order is provided via `Val{TD}` (default,
`TD=3`).
The provided operators must be given in a Cartesian basis.

The assumed spatial dimensionality - which must be consistent with the dimension of the
provided operators `ops` - must be supplied as a `Val(D)` (default, `D=3`).

## Examples

Second harmonic generation is forbidden in inversion symmetric materials
```
julia> using Neumann, Crystalline

julia> TD = 3 # second-harmonic generation → third-rank tensor
julia> inversion = S"-x,-y,-z"
julia> neumann_relations(rotation(inversion), Val(TD))
1-element Vector{String}:
 "xxx = yxx = zxx = xyx = yyx = z" ⋯ 102 bytes ⋯ "yz = zyz = xzz = yzz = zzz = 0"
```

But several components of the second-harmonic response are allowed under e.g. the point
group symmetry of D₃ (321):
```
julia> ops = generators("321", PointGroup{3}) # get relevant operators from Crystalline
julia> ops .= cartesianize.(ops, Ref(crystal(1,1,1,π/2,π/2,2π/3))); # convert to a Cartesian basis
julia> println.(neumann_relations(ops));
yxx = xyx = yyx = xxy = yxy = xyy
zyx = -zxy
yzx = -xzy
yxz = -xyz
xxx = zxx = xzx = zzx = yyy = zyy = yzy = zzy = xxz = zxz = yyz = zyz = xzz = yzz = zzz = 0
```
"""
function neumann_relations(ops, TDᵛ::Val{TD}=Val(3), Dᵛ::Val{D}=Val(3);
                           atol::Real=1e-10, neumann_kws...) where {TD, D}
    A = neumann(ops, TDᵛ, Dᵛ; neumann_kws...)
    return extract_relations(A, TDᵛ, Dᵛ; atol)
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