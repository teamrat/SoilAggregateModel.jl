# initial_conditions.jl — ONE CHANGE for C_aq O₂ state variable
#
# In create_initial_state(), Step 6:
#
#   BEFORE:
#     K_H = K_H_O2(ic.T_0)
#     O_total = ic.O2_gas * (θ_0 + K_H * θ_a_0) / K_H
#
#   AFTER:
#     K_H = K_H_O2(ic.T_0)
#     O_aq = ic.O2_gas / K_H
#
# And in Step 9 (Fill state):
#
#   BEFORE:
#     state.O .= O_total
#
#   AFTER:
#     state.O .= O_aq
