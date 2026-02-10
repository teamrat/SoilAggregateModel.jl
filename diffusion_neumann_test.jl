using SoilAggregateModel
import SoilAggregateModel: crank_nicolson_step!, neumann_flux, neumann_zero, GridInfo

# Set up grid matching Theme 1
n = 100
r_0 = 0.25  # POM radius for 0.5mm diameter
r_max = 5.0
grid = GridInfo(n, r_0, r_max)

# Initial state: uniform C = 1.0
u = ones(n)

# Compute initial mass using conservation weights
mass_initial = sum(u[i] * grid.W[i] for i in 1:n)

# Uniform diffusion coefficient
D = fill(0.01, n)

# Known flux density at inner boundary
J_P = 10.0  # μg/mm²/day

# Workspace
lower = zeros(n-1)
diag_vec = zeros(n)
upper_vec = zeros(n-1)
rhs_vec = zeros(n)

# One CN half-step
dt_half = 0.5  # days

crank_nicolson_step!(u, D, dt_half, grid.r_grid, grid.h,
                     lower, diag_vec, upper_vec, rhs_vec,
                     neumann_flux, neumann_zero, J_P, 0.0)

# Compute final mass
mass_final = sum(u[i] * grid.W[i] for i in 1:n)

# Expected mass increase: J_P × surface_area × dt_half
surface_area = 4.0 * π * r_0^2
expected_increase = J_P * surface_area * dt_half

actual_increase = mass_final - mass_initial
error = actual_increase - expected_increase

println("Surface area: $(surface_area)")
println("Initial mass: $(mass_initial)")
println("Final mass:   $(mass_final)")
println("Actual ΔM:    $(actual_increase)")
println("Expected ΔM:  $(expected_increase)")
println("Error:        $(error)")
println("Rel error:    $(error / expected_increase)")
println()

# Now test 100 steps to check if error accumulates
u2 = ones(n)
for step in 1:100
    crank_nicolson_step!(u2, D, dt_half, grid.r_grid, grid.h,
                         lower, diag_vec, upper_vec, rhs_vec,
                         neumann_flux, neumann_zero, J_P, 0.0)
end
mass_100 = sum(u2[i] * grid.W[i] for i in 1:n)
expected_100 = mass_initial + 100 * expected_increase
error_100 = mass_100 - expected_100
println("After 100 steps:")
println("  Mass:       $(mass_100)")
println("  Expected:   $(expected_100)")
println("  Error:      $(error_100)")
println("  Rel error:  $(error_100 / (100 * expected_increase))")
