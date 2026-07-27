# reactions.jl — TWO CHANGES for C_aq O₂ state variable
#
# In compute_source_terms():
#
# CHANGE 1 — O_aq computation (STEP 1)
#
#   BEFORE:
#     O_aq = O * θ / (θ + temp_cache.K_H_O * θ_a)
#
#   AFTER:
#     O_aq = O    # State variable is already aqueous concentration
#
#
# CHANGE 2 — S_O computation (STEP 7)
#
#   BEFORE:
#     S_O_val = -bio.α_O * Resp_total_val
#
#   AFTER:
#     capacity_O = θ + temp_cache.K_H_O * θ_a
#     S_O_val = -bio.α_O * Resp_total_val / capacity_O
#
#
# CHANGE 3 — Update docstring for O argument
#
#   BEFORE:
#     - `O`: Total oxygen [μg/mm³]
#
#   AFTER:
#     - `O`: Aqueous oxygen concentration [μg/mm³]
