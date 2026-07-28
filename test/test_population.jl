# test_population.jl
# Tests for src/postprocessing/population.jl
#
# These are analytic tests: every expected value below is computed by hand from
# the formulas in the docstrings, not copied from a run.

using Test

@testset "Population upscaling" begin

    @testset "sieve_class" begin
        s = [0.053, 0.25, 2.0]
        @test sieve_class(0.0, s)    == 1
        @test sieve_class(0.052, s)  == 1
        @test sieve_class(0.053, s)  == 2      # at the aperture -> retained
        @test sieve_class(0.24, s)   == 2
        @test sieve_class(0.25, s)   == 3
        @test sieve_class(1.999, s)  == 3
        @test sieve_class(2.0, s)    == 4
        @test sieve_class(9.0, s)    == 4

        # No sieves -> one class
        @test sieve_class(1.0, Float64[]) == 1
    end

    @testset "aggregate_mass" begin
        # Shell only: r_0 = 0.5, D = 2.0 -> r = 1.0
        #   V = (4/3)π(1³ - 0.5³) = (4/3)π(0.875)
        ρ_b = 1300.0
        V_shell = (4.0/3.0) * π * (1.0^3 - 0.5^3)
        @test aggregate_mass(2.0, 0.5, 0.0; ρ_b=ρ_b, f_C_POM=0.443) ≈ ρ_b * V_shell

        # Core only: D equal to the POM diameter -> zero shell
        @test aggregate_mass(1.0, 0.5, 44.3; ρ_b=ρ_b, f_C_POM=0.443) ≈ 100.0

        # Shell floors at zero when D_agg < 2·r_0
        @test aggregate_mass(0.5, 0.5, 0.0; ρ_b=ρ_b, f_C_POM=0.443) == 0.0

        @test_throws ArgumentError aggregate_mass(2.0, 0.5, 1.0; ρ_b=ρ_b, f_C_POM=0.0)
    end

    @testset "population_statistics — single class, hand-checked" begin
        # One size class, one time. r_0 = 0.5 mm, D_agg = 2.0 mm, no POM left.
        # n = 1 particle in a 1 mm³ reference volume with ρ_b = 1000 µg/mm³.
        D    = reshape([2.0], 1, 1)
        P    = reshape([0.0], 1, 1)
        CO2  = reshape([7.0], 1, 1)
        flux = reshape([0.5], 1, 1)

        m = 1000.0 * (4.0/3.0) * π * (1.0^3 - 0.5^3)

        st = population_statistics(D, P, CO2, flux, [0.5], [3.0];
                                   soil_volume_mm3=1.0e6, ρ_b=1000.0, f_C_POM=0.5)

        @test st.MWD_agg_only[1] ≈ 2.0             # only one diameter present
        @test st.POM_mass_frac[1] == 0.0           # no residue left
        @test st.f_agg[1] ≈ 3.0 * m / (1000.0 * 1.0e6)
        @test st.f_agg_vol[1] ≈ 3.0 * (4.0/3.0) * π * 1.0^3 / 1.0e6
        @test st.CO2_total[1] ≈ 21.0               # 3 particles × 7
        @test st.CO2_flux_total[1] ≈ 1.5
        @test isnan(st.MWD_fixed_weight[1])        # no nominals supplied
        @test size(st.class_pct) == (1, 1)         # no sieves -> one class
    end

    @testset "population_statistics — mass weighting is not volume weighting" begin
        # Two classes with the SAME aggregate diameter but different residue
        # left inside. Volume weighting cannot tell them apart; mass weighting
        # must, because the core term differs.
        D    = [2.0 2.0; 2.0 2.0]
        P    = [100.0 0.0; 100.0 100.0]     # class 1 loses its core, class 2 keeps it
        CO2  = zeros(2, 2)
        flux = zeros(2, 2)

        st = population_statistics(D, P, CO2, flux, [0.5, 0.5], [1.0, 1.0];
                                   soil_volume_mm3=1.0e6, ρ_b=1000.0, f_C_POM=0.5)

        @test st.f_agg_vol[1] ≈ st.f_agg_vol[2]    # geometry unchanged
        @test st.f_agg[2] < st.f_agg[1]            # mass fell as residue was respired
        @test st.POM_mass_frac[2] < st.POM_mass_frac[1]
        @test st.MWD_agg_only[1] ≈ 2.0             # one diameter, unaffected
    end

    @testset "class shares sum to 100 % with a mineral distribution" begin
        sieves  = [0.053, 0.25, 2.0]
        mineral = [0.56, 0.188, 0.252, 0.0]

        D    = [0.3 2.5; 1.0 3.0]
        P    = [10.0 5.0; 10.0 5.0]
        CO2  = zeros(2, 2)
        flux = zeros(2, 2)

        st = population_statistics(D, P, CO2, flux, [0.1, 0.2], [1.0e3, 5.0e2];
                                   sieve_sizes=sieves,
                                   mineral_class_fractions=mineral,
                                   class_nominal_mm=[0.0265, 0.1515, 1.125, 5.0],
                                   soil_volume_mm3=1.0e6, ρ_b=1370.0, f_C_POM=0.443)

        for t in 1:2
            @test sum(st.class_pct[:, t]) ≈ 100.0 atol=1e-8
        end

        # Aggregates moved from the 0.25-2 mm class into >2 mm
        @test st.class_pct[4, 2] > st.class_pct[4, 1]

        # Fixed-weight MWD is a convex combination of the nominals, so it is
        # bounded by them whatever the aggregates do.
        for t in 1:2
            @test st.MWD_fixed_weight[t] >= 0.0265
            @test st.MWD_fixed_weight[t] <= 5.0
        end
    end

    @testset "without a mineral distribution the classes sum to f_agg" begin
        sieves = [0.25]
        D    = reshape([1.0, 0.1], 2, 1)
        P    = zeros(2, 1)
        CO2  = zeros(2, 1)
        flux = zeros(2, 1)

        st = population_statistics(D, P, CO2, flux, [0.05, 0.02], [1.0e3, 1.0e3];
                                   sieve_sizes=sieves,
                                   soil_volume_mm3=1.0e6, ρ_b=1300.0, f_C_POM=0.443)

        @test sum(st.class_pct[:, 1]) ≈ 100.0 * st.f_agg[1] atol=1e-8
    end

    @testset "shell thickness — the continuum-validity diagnostic" begin
        # Two classes, r_0 = 0.25 and 0.5, retained at D = 1.0 and 2.0.
        # Shells are 0.5 - 0.25 = 0.25 and 1.0 - 0.5 = 0.5 mm.
        # Equal aggregate masses would give 0.375; they are not equal, so the
        # expected value is the mass-weighted mean, computed here by hand.
        D    = reshape([1.0, 2.0], 2, 1)
        P    = zeros(2, 1)
        CO2  = zeros(2, 1)
        flux = zeros(2, 1)
        r0   = [0.25, 0.5]
        nd   = [1.0, 1.0]
        ρ_b  = 1000.0

        m1 = ρ_b * (4/3) * π * (0.5^3 - 0.25^3)
        m2 = ρ_b * (4/3) * π * (1.0^3 - 0.5^3)
        expected = (m1 * 0.25 + m2 * 0.5) / (m1 + m2)

        st = population_statistics(D, P, CO2, flux, r0, nd;
                                   soil_volume_mm3=1.0e6, ρ_b=ρ_b, f_C_POM=0.5)
        @test st.shell_thickness_mm[1] ≈ expected

        # A shell can never be negative even if r_agg were to sit below r_0
        D2 = reshape([0.2, 0.2], 2, 1)
        st2 = population_statistics(D2, P, CO2, flux, r0, nd;
                                    soil_volume_mm3=1.0e6, ρ_b=ρ_b, f_C_POM=0.5)
        @test st2.shell_thickness_mm[1] >= 0.0
    end

    @testset "f_agg is reported unclamped when aggregates exceed the sample" begin
        # Absurdly dense population: aggregate mass must exceed sample mass.
        # f_agg must SHOW that (> 1), the mineral top-up must not go negative,
        # and the class shares must stay non-negative.
        D    = reshape([2.0], 1, 1)
        P    = reshape([0.0], 1, 1)
        CO2  = zeros(1, 1)
        flux = zeros(1, 1)

        st = (@test_logs (:warn,) match_mode=:any population_statistics(
                  D, P, CO2, flux, [0.5], [1.0e9];
                  sieve_sizes=[0.25],
                  mineral_class_fractions=[0.6, 0.4],
                  soil_volume_mm3=1.0e6, ρ_b=1000.0, f_C_POM=0.5))

        @test st.f_agg[1] > 1.0                      # visible, not hidden
        @test all(st.class_pct[:, 1] .>= 0.0)        # top-up floored at zero
        @test sum(st.class_pct[:, 1]) > 100.0        # and the excess is evident
    end

    @testset "cell occupancy — per-class share of owned soil" begin
        # Two classes. Aggregate volumes are (4/3)pi(D/2)^3 = 0.5236 and 4.18879.
        # Cells supplied as 2.0 and 4.18879 mm^3 -> occupancy 0.2618 and 1.0.
        D    = reshape([1.0, 2.0], 2, 1)
        P    = zeros(2, 1)
        CO2  = zeros(2, 1)
        flux = zeros(2, 1)
        cells = [2.0, (4/3)*π]

        st = population_statistics(D, P, CO2, flux, [0.25, 0.5], [1.0, 1.0];
                                   cell_volume_mm3=cells,
                                   soil_volume_mm3=1.0e6, ρ_b=1000.0, f_C_POM=0.5)

        @test st.cell_occupancy[1, 1] ≈ (4/3)*π*0.5^3 / 2.0
        @test st.cell_occupancy[2, 1] ≈ 1.0          # class 2 exactly fills its cell
        @test st.cell_occupancy[2, 1] >= st.cell_occupancy[1, 1]

        # Without cell volumes the diagnostic is NaN, not silently zero
        st2 = population_statistics(D, P, CO2, flux, [0.25, 0.5], [1.0, 1.0];
                                    soil_volume_mm3=1.0e6, ρ_b=1000.0, f_C_POM=0.5)
        @test all(isnan, st2.cell_occupancy)

        @test_throws DimensionMismatch population_statistics(
            D, P, CO2, flux, [0.25, 0.5], [1.0, 1.0];
            cell_volume_mm3=[1.0], soil_volume_mm3=1.0e6, ρ_b=1000.0, f_C_POM=0.5)
    end

    @testset "input validation" begin
        D    = zeros(2, 3)
        P    = zeros(2, 3)
        CO2  = zeros(2, 3)
        flux = zeros(2, 3)
        r0   = [0.1, 0.2]
        n    = [1.0, 1.0]
        kw   = (ρ_b=1300.0, f_C_POM=0.443)

        @test_throws DimensionMismatch population_statistics(D, zeros(3,3), CO2, flux, r0, n; kw...)
        @test_throws DimensionMismatch population_statistics(D, P, CO2, flux, [0.1], n; kw...)
        @test_throws DimensionMismatch population_statistics(D, P, CO2, flux, r0, [1.0]; kw...)
        @test_throws ArgumentError population_statistics(D, P, CO2, flux, r0, n; ρ_b=0.0, f_C_POM=0.443)
        @test_throws ArgumentError population_statistics(D, P, CO2, flux, r0, n; sieve_sizes=[2.0, 0.25], kw...)

        # mineral_class_fractions must match the class count and sum to 1
        @test_throws ArgumentError population_statistics(D, P, CO2, flux, r0, n;
            sieve_sizes=[0.25], mineral_class_fractions=[0.5, 0.4], kw...)
        @test_throws ArgumentError population_statistics(D, P, CO2, flux, r0, n;
            sieve_sizes=[0.25], mineral_class_fractions=[1.0], kw...)

        # nominals without a mineral distribution are not a whole-sample MWD
        @test_throws ArgumentError population_statistics(D, P, CO2, flux, r0, n;
            sieve_sizes=[0.25], class_nominal_mm=[0.1, 1.0], kw...)
    end
end
