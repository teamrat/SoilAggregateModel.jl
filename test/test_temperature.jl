"""
Tests for temperature/*.jl modules
"""

import SoilAggregateModel: arrhenius, arrhenius_ratio, q10_equivalent,
                           water_viscosity, D_DOC_water, D_O2_water, D_O2_air,
                           D_fungal_translocation, henry_vant_hoff, K_H_O2, O2_saturation

@testset "Arrhenius" begin
    @testset "Basic behavior" begin
        # At reference temperature, factor should be 1.0
        @test arrhenius(60_000.0, 293.15, 293.15) ≈ 1.0

        # Higher temperature → factor > 1
        @test arrhenius(60_000.0, 298.15, 293.15) > 1.0

        # Lower temperature → factor < 1
        @test arrhenius(60_000.0, 288.15, 293.15) < 1.0

        # Higher Ea → stronger temperature sensitivity
        f_high = arrhenius(80_000.0, 298.15, 293.15)
        f_low  = arrhenius(40_000.0, 298.15, 293.15)
        @test f_high > f_low
    end

    @testset "Typical biological values" begin
        # Bacterial metabolism: Ea = 60 kJ/mol
        f_bac = arrhenius(60_000.0, 298.15, 293.15)
        @test 1.3 < f_bac < 1.6  # ~40-50% increase over 5 K

        # Fungal metabolism: Ea = 55 kJ/mol
        f_fun = arrhenius(55_000.0, 298.15, 293.15)
        @test 1.25 < f_fun < 1.5

        # MAOC sorption: Ea = 25 kJ/mol (lower)
        f_maoc = arrhenius(25_000.0, 298.15, 293.15)
        @test 1.1 < f_maoc < 1.2  # weaker temperature effect
    end

    @testset "Arrhenius ratio" begin
        # Should equal ratio of individual factors
        Ea = 60_000.0
        T1 = 293.15
        T2 = 303.15

        ratio_direct = arrhenius_ratio(Ea, T1, T2)
        ratio_manual = arrhenius(Ea, T2, 293.15) / arrhenius(Ea, T1, 293.15)

        @test ratio_direct ≈ ratio_manual rtol=1e-10

        # Symmetric property
        @test arrhenius_ratio(Ea, T1, T2) ≈ 1.0 / arrhenius_ratio(Ea, T2, T1)
    end

    @testset "Q10 equivalent" begin
        # Typical biological Q10 should be 2-3
        q10_bac = q10_equivalent(60_000.0, 293.15)
        @test 2.0 < q10_bac < 3.0

        # Higher Ea → higher Q10
        q10_high = q10_equivalent(80_000.0, 293.15)
        q10_low  = q10_equivalent(40_000.0, 293.15)
        @test q10_high > q10_low

        # Q10 decreases slightly with temperature
        q10_cold = q10_equivalent(60_000.0, 283.15)
        q10_warm = q10_equivalent(60_000.0, 303.15)
        @test q10_cold > q10_warm
    end
end

@testset "Water viscosity (VFT)" begin
    @testset "Known values" begin
        # At 20°C (293.15 K): η ≈ 1.002 mPa·s
        η_20 = water_viscosity(293.15)
        @test η_20 ≈ 1.002 rtol=0.02

        # At 25°C (298.15 K): η ≈ 0.890 mPa·s
        η_25 = water_viscosity(298.15)
        @test η_25 ≈ 0.890 rtol=0.02

        # At 0°C (273.15 K): η ≈ 1.73 mPa·s (VFT gives ~1.73, NIST is 1.787)
        η_0 = water_viscosity(273.15)
        @test η_0 ≈ 1.73 rtol=0.03
    end

    @testset "Temperature trends" begin
        # Viscosity decreases with warming
        η_cold = water_viscosity(273.15)
        η_warm = water_viscosity(313.15)
        @test η_warm < η_cold

        # Should be monotonic
        temps = 273.15:5.0:313.15
        viscosities = [water_viscosity(T) for T in temps]
        @test all(diff(viscosities) .< 0)  # Strictly decreasing
    end
end

@testset "DOC diffusion (Stokes-Einstein)" begin
    @testset "Reference condition" begin
        # At T_ref, should return D_ref
        D_ref = 0.864  # mm²/day
        T_ref = 293.15
        @test D_DOC_water(T_ref, D_ref, T_ref) ≈ D_ref
    end

    @testset "Temperature trends" begin
        D_ref = 0.864
        T_ref = 293.15

        # Warming increases diffusion
        D_cold = D_DOC_water(283.15, D_ref, T_ref)
        D_warm = D_DOC_water(303.15, D_ref, T_ref)
        @test D_warm > D_ref > D_cold

        # Roughly 2-3% per degree near 20°C
        D_plus5 = D_DOC_water(298.15, D_ref, T_ref)
        percent_change = 100 * (D_plus5 / D_ref - 1.0)
        @test 10.0 < percent_change < 20.0  # ~10-20% over 5 K
    end
end

@testset "O₂ diffusion in water (Han-Bartels)" begin
    @testset "Known values" begin
        # At 20°C (293.15 K): D ≈ 150 mm²/day
        # (1.75e-5 cm²/s × 8.64e6 = 151.2 mm²/day)
        D_20 = D_O2_water(293.15)
        @test 145 < D_20 < 160

        # At 25°C (298.15 K): D ≈ 172 mm²/day
        D_25 = D_O2_water(298.15)
        @test 170 < D_25 < 180

        # Should be faster than DOC (smaller molecule)
        D_DOC = 0.864
        @test D_20 > D_DOC
    end

    @testset "Temperature trends" begin
        # Warming increases diffusion
        D_cold = D_O2_water(278.15)
        D_warm = D_O2_water(308.15)
        @test D_warm > D_cold

        # Should be monotonic
        temps = 273.15:5.0:313.15
        diffusivities = [D_O2_water(T) for T in temps]
        @test all(diff(diffusivities) .> 0)  # Strictly increasing
    end
end

@testset "O₂ diffusion in air (Chapman-Enskog)" begin
    @testset "Reference condition" begin
        D_ref = 1.73e6  # mm²/day
        T_ref = 293.15
        @test D_O2_air(T_ref, D_ref, T_ref) ≈ D_ref
    end

    @testset "Temperature scaling" begin
        D_ref = 1.73e6
        T_ref = 293.15

        # Should scale as T^1.75
        T1 = 280.0
        T2 = 320.0
        D1 = D_O2_air(T1, D_ref, T_ref)
        D2 = D_O2_air(T2, D_ref, T_ref)

        expected_ratio = (T2 / T1)^1.75
        actual_ratio = D2 / D1
        @test actual_ratio ≈ expected_ratio rtol=1e-10

        # Much faster than in water
        D_air = D_O2_air(293.15, D_ref, T_ref)
        D_water = D_O2_water(293.15)
        @test D_air > 1000 * D_water
    end
end

@testset "Fungal translocation" begin
    @testset "Arrhenius behavior" begin
        D_ref = 1.0  # mm²/day
        Ea_F = 55_000.0
        T_ref = 293.15

        # At reference temperature
        @test D_fungal_translocation(T_ref, D_ref, Ea_F, T_ref) ≈ D_ref

        # Warming increases translocation (biological process)
        D_warm = D_fungal_translocation(298.15, D_ref, Ea_F, T_ref)
        @test D_warm > D_ref

        # Should match Arrhenius factor
        f = arrhenius(Ea_F, 298.15, T_ref)
        @test D_fungal_translocation(298.15, D_ref, Ea_F, T_ref) ≈ D_ref * f
    end
end

@testset "Henry's law (van't Hoff)" begin
    @testset "Reference condition" begin
        K_H_ref = 31.25
        T_ref = 298.15
        @test henry_vant_hoff(K_H_ref, -12_000.0, T_ref, T_ref) ≈ K_H_ref
    end

    @testset "O₂ temperature trends" begin
        # ΔH_sol < 0 (exothermic dissolution)
        # Warming → K_H increases → C_aq = C_gas/K_H decreases (LESS soluble)
        # Cooling → K_H decreases → C_aq = C_gas/K_H increases (MORE soluble)
        K_H_cold = henry_vant_hoff(31.25, -12_000.0, 283.15, 298.15)
        K_H_warm = henry_vant_hoff(31.25, -12_000.0, 313.15, 298.15)

        # Cold water should have LOWER K_H (more dissolved O₂)
        @test K_H_cold < K_H_warm
        @test K_H_cold < 31.25  # Colder than reference
        @test K_H_warm > 31.25  # Warmer than reference
    end

    @testset "K_H_O2 convenience function" begin
        # Should use standard values
        K_H = K_H_O2(298.15)
        @test K_H ≈ 31.25

        # Cold water → lower K_H → higher solubility (correct physics)
        K_H_cold = K_H_O2(273.15)
        K_H_warm = K_H_O2(313.15)
        @test K_H_cold < K_H_warm

        # Typical range
        @test 20.0 < K_H_cold < 30.0
        @test 33.0 < K_H_warm < 40.0
    end

    @testset "Sign convention check" begin
        K_H_ref = 31.25
        T_ref = 298.15

        # Exothermic dissolution (ΔH_sol < 0): warming INCREASES K_H (less soluble)
        K_H_exo_cold = henry_vant_hoff(K_H_ref, -12_000.0, 283.15, T_ref)
        K_H_exo_warm = henry_vant_hoff(K_H_ref, -12_000.0, 313.15, T_ref)
        @test K_H_exo_warm > K_H_exo_cold  # Higher K_H = less dissolved O₂

        # Endothermic dissolution (ΔH_sol > 0): warming DECREASES K_H (more soluble)
        K_H_endo_cold = henry_vant_hoff(K_H_ref, 12_000.0, 283.15, T_ref)
        K_H_endo_warm = henry_vant_hoff(K_H_ref, 12_000.0, 313.15, T_ref)
        @test K_H_endo_warm < K_H_endo_cold  # Lower K_H = more dissolved
    end
end

@testset "O₂ saturation" begin
    @testset "Typical values" begin
        # At 20°C, sea level
        # C_gas ≈ 0.28 μg/mm³, K_H ≈ 29 → C_aq ≈ 0.28/29 ≈ 0.01 μg/mm³
        C_sat_20 = O2_saturation(293.15)
        @test 0.008 < C_sat_20 < 0.012

        # Cold water holds MORE O₂ (correct physics)
        C_sat_0 = O2_saturation(273.15)
        C_sat_30 = O2_saturation(303.15)
        @test C_sat_0 > C_sat_20 > C_sat_30
    end

    @testset "Altitude effects" begin
        # Higher altitude → lower partial pressure → lower saturation
        C_sea = O2_saturation(293.15, 0.21)
        C_high = O2_saturation(293.15, 0.15)
        @test C_high < C_sea

        # Should scale linearly with P_atm
        ratio = C_high / C_sea
        @test ratio ≈ 0.15 / 0.21 rtol=0.01
    end
end

@testset "Temperature function consistency" begin
    @testset "All positive and finite" begin
        T_range = 273.15:5.0:313.15

        for T in T_range
            # Arrhenius
            @test arrhenius(60_000.0, T, 293.15) > 0
            @test isfinite(arrhenius(60_000.0, T, 293.15))

            # Viscosity
            @test water_viscosity(T) > 0
            @test isfinite(water_viscosity(T))

            # Diffusivities
            @test D_DOC_water(T, 0.864, 293.15) > 0
            @test D_O2_water(T) > 0
            @test D_O2_air(T, 1.73e6, 293.15) > 0
            @test isfinite(D_O2_water(T))

            # Henry's law
            @test K_H_O2(T) > 0
            @test isfinite(K_H_O2(T))
        end
    end

    @testset "Type stability" begin
        # All functions should be type-stable
        @inferred arrhenius(60_000.0, 298.15, 293.15)
        @inferred water_viscosity(293.15)
        @inferred D_DOC_water(298.15, 0.864, 293.15)
        @inferred D_O2_water(293.15)
        @inferred D_O2_air(298.15, 1.73e6, 293.15)
        @inferred K_H_O2(293.15)
        @inferred O2_saturation(293.15)
    end
end
