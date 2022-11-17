var documenterSearchIndex = {"docs":
[{"location":"","page":"API","title":"API","text":"CurrentModule = Neumann","category":"page"},{"location":"#Neumann.jl","page":"API","title":"Neumann.jl","text":"","category":"section"},{"location":"","page":"API","title":"API","text":"Neumann.jl exports a single method:","category":"page"},{"location":"","page":"API","title":"API","text":"Modules = [Neumann]\nPrivate = false\nOrder   = [:function]","category":"page"},{"location":"#Neumann.neumann-Union{Tuple{D}, Tuple{N}, Tuple{Any, Val{N}, Val{D}}} where {N, D}","page":"API","title":"Neumann.neumann","text":"neumann(ops, N::Integer; kws...)  -->  Vector{String}\n\nReturn the symmetry-constraints among components of a tensor imposed by the set of point group symmetries in ops. The tensor order is provided via N. ops can be either a single symmetry operation (provided as an AbstractMatrix) or any iterable of symmetry operations.\n\nThe obtained relations are returned as a Vector{String}, with each element of this vector giving either:\n\na free, nonzero component of the tensor (e.g., \"xxx\"),\nfree, nonzero - but interrelated - components of the tensor (e.g., \"xxx = yyy\"),\nforbidden, zero components of the tensor (e.g., \"xyz = xzy = 0\").\n\nNote that if the provided operators are not given in a Cartesian basis, the returned relations among tensor components will match the basis of the provided operators: we generally recommend supplying operators in a Cartesian basis for ease of interpretation.\n\nKeyword arguments kws\n\nkleinman (default, false): whether to incorporate Kleinman symmetry, i.e. whether to enforce that A[i,j,k,…] = A[i,perm(j,k,…)...] with perm(j,k,…) denoting all unique permutations of [j,k,…]. Relevant e.g., to low-frequency second-harmonic generation.\nsparsify (default, true): whether to attempt to sparsify relations between allowed nonzero components of the response tensor. A \"poor man's\" approximation of matrix sparsification is carried out using the reduced row echelon form.\nrref_tol (default, 1e-11): the absolute tolerance used during row echelon reduction in sparsification. As the row echelon form is numerically unstable for larger matrices, it may be necessary to increase this tolerance for larger tensor orders. Set to nothing to use the default tolerance of of RowEchelon.jl's rref.\natol (default, 1e-10): the absolute tolerance used in assessing whether a term in a relation is considered vanishing or not (and similarly used to assess whether a term has integer coefficient). Must be greater than rref_tol.\nnullspace_kws (default, empty): keyword arguments passed to nullspace.\n\nExamples\n\nSecond harmonic generation is forbidden in inversion symmetric materials\n\njulia> using Neumann\njulia> N = 3 # second-harmonic generation → third-rank tensor\njulia> inversion = [-1 0 0; 0 -1 0; 0 0 -1]\njulia> neumann(inversion, N)\n1-element Vector{String}:\n \"xxx = yxx = zxx = xyx = yyx = z\" ⋯ 102 bytes ⋯ \"yz = zyz = xzz = yzz = zzz = 0\"\n\nBut several components of the second-harmonic response are allowed under e.g. the point group symmetry of D₃ (321):\n\njulia> using Crystalline\njulia> ops = generators(\"321\", PointGroup{3}) # obtain generators of D₃ (321) from Crystalline\njulia> ops .= cartesianize.(ops, Ref(crystal(1,1,1,π/2,π/2,2π/3))); # convert to a Cartesian basis\njulia> neumann(ops, N)\n5-element Vector{String}:\n \"xxx = -yyx = -yxy = -xyy\"\n \"zyx = -zxy\"\n \"yzx = -xzy\"\n \"yxz = -xyz\"\n \"yxx = zxx = xyx = xzx = zzx = x\" ⋯ 41 bytes ⋯ \"yyz = zyz = xzz = yzz = zzz = 0\"\n\n\n\n\n\n","category":"method"}]
}
