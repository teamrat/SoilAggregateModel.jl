# test_particle_size.jl
# Tests for src/physics/particle_size.jl
#
# The Sauter mean is defined by preserving the volume-to-surface ratio, so the
# strongest tests are the ones that check that property directly rather than
# re-running the formula.

using Test

@testset "Sauter mean diameter" begin

    @testset "monodisperse returns the diameter itself" begin
        for d in (0.001, 0.01, 0.5, 2.0)
            @test sauter_mean_diameter([1.0], [d]) ≈ d
        end
    end

    @testset "preserves volume-to-surface ratio — the defining property" begin
        # For spheres, surface per unit solid volume is S_v = 6 Σ fᵢ/dᵢ.
        # A single sphere of diameter d₃₂ must have the same S_v = 6/d₃₂.
        f = [0.44, 0.45, 0.11]
        d = [0.316, 0.010, 0.001]
        S_v_mix = 6.0 * sum(f[i] / d[i] for i in 1:3)
        d32 = sauter_mean_diameter(f, d)
        @test 6.0 / d32 ≈ S_v_mix
    end

    @testset "bounded by the extreme diameters" begin
        f = [0.3, 0.4, 0.3]
        d = [0.316, 0.010, 0.001]
        d32 = sauter_mean_diameter(f, d)
        @test d32 > minimum(d)
        @test d32 < maximum(d)
    end

    @testset "dominated by the fine fraction, not the coarse" begin
        # This is the whole reason for preferring it to a geometric mean.
        # 90% coarse sand by MASS, 10% clay: d32 must still sit near the clay.
        d32 = sauter_mean_diameter([0.90, 0.10], [0.316, 0.001])
        @test d32 < 0.011           # nowhere near the 0.316 mass-dominant class
        # arithmetic mean would give 0.285 — two orders of magnitude larger
        @test d32 < 0.05 * (0.90 * 0.316 + 0.10 * 0.001)
    end

    @testset "adding fines always lowers d32" begin
        base = sauter_mean_diameter([1.0, 0.0 + eps()], [0.316, 0.001])
        prev = base
        for f_clay in (0.02, 0.05, 0.10, 0.25)
            d32 = sauter_mean_diameter([1.0 - f_clay, f_clay], [0.316, 0.001])
            @test d32 < prev
            prev = d32
        end
    end

    @testset "sauter_from_texture — De Gryze soils" begin
        # Hand-checked: d32 = 1/(sand/0.316 + silt/0.010 + clay/0.001)
        # soil 3: 1/(0.44/0.316 + 0.45/0.01 + 0.11/0.001) = 1/(1.392+45+110)
        @test sauter_from_texture(0.44, 0.45, 0.11) ≈
              1.0 / (0.44/0.316 + 0.45/0.010 + 0.11/0.001)
        @test sauter_from_texture(0.44, 0.45, 0.11) ≈ 0.006394 rtol=1e-3

        # The five soils must order finest-to-coarsest as 5 < 4 < 2 ≈ 3 < 1
        d = [sauter_from_texture(s, si, c) for (s, si, c) in
             ((0.53, 0.40, 0.07), (0.33, 0.57, 0.10), (0.44, 0.45, 0.11),
              (0.44, 0.43, 0.13), (0.21, 0.52, 0.27))]
        @test d[5] < d[4] < d[3] < d[1]
        @test d[5] < d[2] < d[1]
    end

    @testset "ordering is robust to the clay class midpoint" begin
        # The spread changes, the ranking must not — this is what lets the
        # texture comparison be reported at all (see the docstring).
        for d_clay in (0.0005, 0.001, 0.002)
            d1 = sauter_from_texture(0.53, 0.40, 0.07; d_clay=d_clay)
            d4 = sauter_from_texture(0.44, 0.43, 0.13; d_clay=d_clay)
            d5 = sauter_from_texture(0.21, 0.52, 0.27; d_clay=d_clay)
            @test d5 < d4
            @test d4 < d1
        end
    end

    @testset "input validation" begin
        @test_throws DimensionMismatch sauter_mean_diameter([0.5, 0.5], [0.1])
        @test_throws ArgumentError sauter_mean_diameter(Float64[], Float64[])
        @test_throws ArgumentError sauter_mean_diameter([0.5, 0.4], [0.1, 0.2])   # sums to 0.9
        @test_throws ArgumentError sauter_mean_diameter([1.5, -0.5], [0.1, 0.2])  # negative
        @test_throws ArgumentError sauter_mean_diameter([0.5, 0.5], [0.1, 0.0])   # zero diameter
        @test_throws ArgumentError sauter_from_texture(0.5, 0.3, 0.1)             # sums to 0.9
    end

    @testset "wired into SoilProperties and the stability threshold" begin
        # G_c = τ_w·d_32/κ_b must rise linearly with d_32: a coarser soil needs
        # more binder. This is the direction Albalasmeh & Ghezzehei report.
        fine   = SoilProperties(d_32 = 0.003)
        coarse = SoilProperties(d_32 = 0.009)
        @test coarse.d_32 ≈ 3.0 * fine.d_32      # ≈, not ==: 3.0*0.003 = 0.009000000000000001
        @test fine.κ_b == coarse.κ_b          # same binder, different texture
        @test fine.w_E > 0.0
    end
end
