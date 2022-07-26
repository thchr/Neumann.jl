using Neumann
using Test
using Crystalline
using LinearAlgebra

@testset "Neumann.jl" begin
    @testset "Centrosymmetry" begin
        N = 3 # second-harmonic generation → third-rank tensor
        inversion = [-1 0 0; 0 -1 0; 0 0 -1]
        for N in 1:7
            rels = neumann(inversion, N)
            if iseven(N)
                @test length(rels) == 3^N
                @test all(s->length(s)==N, rels)
            else # isodd(N)
                @test length(rels) == 1
                @test only(rels) == join(Neumann.xyz_sorting(Val(N), Val(3)), " = ")*" = 0"
            end
        end
    end

    @testset "Kleinmann" begin
        # second-rank tensor/matrix (Kleinmann symmetry has no effect)
        for D in 1:4
            rels₂ = neumann(I(D), 2, kleinmann=true)
            @test length(rels₂) == D^2
            @test rels₂ == Neumann.xyz_sorting(Val(2), Val(D))
        end
        
        # third-rank tensor (Kleinmann symmetry couples several components)
        rels₃ = neumann(I(3), 3, kleinmann=true)
        @test length(rels₃) == 18
    end

    @testset "README examples" begin
        # C₄
        ops_C₄ = generators("4", PointGroup{3}) # generators of the group C₄ (4)
        rels_C₄ = neumann(ops_C₄, 3)
        @test rels_C₄ == ["zxx = zyy", "zyx = -zxy", "xzx = yzy", "yzx = -xzy", "xxz = yyz",
                          "yxz = -xyz", "zzz", 
                          "xxx = yxx = xyx = yyx = zzx = xxy = yxy = xyy = yyy = zzy = "*
                          "zxz = zyz = xzz = yzz = 0"]
        # C₃
        ops_C₃ = generators("3", PointGroup{3})
        Rs = crystal(1,1,1,π/2,π/2,2π/3)
        ops_C₃′ = cartesianize.(ops_C₃, Ref(Rs))
        rels_C₃ = neumann(ops_C₃′, 3)
        @test rels_C₃ == ["xxx = -yyx = -yxy = -xyy", "yxx = xyx = xxy = -yyy",
                          "zxx = zyy", "zyx = -zxy", "xzx = yzy", "yzx = -xzy",
                          "xxz = yyz", "yxz = -xyz", "zzz",
                          "zzx = zzy = zxz = zyz = xzz = yzz = 0"]
        # C₄ + Kleinmann
        klein_rels_C₄ = neumann(ops_C₄, 3; kleinmann=true)
        @test klein_rels_C₄ == ["zxx = zyy", "xzx = yzy = xxz = yyz",
                                "yzx = -xzy = yxz = -xyz", "zzz",
                                "xxx = yxx = xyx = yyx = zyx = zzx = xxy = yxy = zxy = "*
                                "xyy = yyy = zzy = zxz = zyz = xzz = yzz = 0"]
        # C₃ + Kleinmann
        klein_rels_C₃ = neumann(ops_C₃′, 3; kleinmann=true)
        @test klein_rels_C₃ == ["xxx = -yyx = -yxy = -xyy", "yxx = xyx = xxy = -yyy",
                                "zxx = zyy", "xzx = yzy = xxz = yyz",
                                "yzx = -xzy = yxz = -xyz", "zzz",
                                "zyx = zzx = zxy = zzy = zxz = zyz = xzz = yzz = 0"]
    end

    @testset "Miscellaneous" begin
        ops = [[0 -1 0; 1 0 0; 0 0 1], [0 -1 0; -1 0 0; 0 0 1]]

        # input as tuple, iterator, and view
        @test (neumann(ops, 3) == neumann(Tuple(ops), 3)
                               == neumann((op for op in ops), 3)
                               == neumann((@view ops[1:end]), 3))

        # span identical column space with/without sparsification
        for N in 2:6
            A = Neumann.neumann_matrix(ops, Val(N), Val(3); sparsify=true)
            B = Neumann.neumann_matrix(ops, Val(N), Val(3); sparsify=false)
            A_to_B, B_to_A = A\B, B\A
            @test A*A_to_B ≈ B
            @test B*B_to_A ≈ A
        end

        # check: sparsification is essentially insensitive to reasonable `rref_tol` choices
        @test (neumann(ops, 3; rref_tol=1e-12, atol=1e-11) ==
               neumann(ops, 3; rref_tol=1e-5, atol=1e-4))

        # check that sparsification succeeds, even for quite large tensor orders
        for N in 1:7
            rels = neumann(ops, N)
            # check that we only get integer coefficients
            @test all(rel -> !occursin(r"[0-9]", rel), @views rels[1:end-1])
            @test endswith(rels[end], " = 0")
        end
    end

    @testset "Error handling" begin
        @test_throws ErrorException neumann(I(3), 3; atol=1e-6, rref_tol=1e-5)
        @test_throws ErrorException neumann([I(3), [0 1; 1 0]], 3)
    end
end
