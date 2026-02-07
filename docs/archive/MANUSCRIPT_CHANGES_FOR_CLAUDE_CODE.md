# Manuscript Changes from Audit — Implementation Update

The manuscript has been updated in three places. If you have already written code based on the earlier version, these are the specific changes needed.

---

## 1. POM Dissolution: R_P vs J_P (new notation, same physics)

The manuscript now defines two quantities to eliminate a units ambiguity:

**J_P** [μg-C/mm²/day] — flux density at the POM surface:
```
J_P = R_P^max(T) · (P/P₀) · B₀/(K_{B,P}+B₀) · F_{n,0}/(K_{F,P}+F_{n,0})
      × θ₀/(θ_P+θ₀) · O_{aq,0}/(L_P+O_{aq,0})
```

**R_P** [μg-C/day] — total dissolution rate:
```
R_P = 4π r₀² · J_P
```

R_P^max has units [μg-C/mm²/day] (surface-specific rate).

**Where each is used:**
- POM scalar ODE: `dP/dt = −R_P = −4πr₀²·J_P`
- Neumann BC for C at r = r₀: `−D_C (∂C/∂r)|_{r₀} = J_P`
- CO₂ diagnostics: R_P is total carbon entering the domain per day

**What to check in code:** The ghost-node implementation for the inner Neumann BC must use J_P (without the 4πr₀² prefactor). The POM ODE must use R_P (with the 4πr₀² prefactor). If your code already computes R_P with the 4πr₀² included, just make sure the BC divides it back out, or compute J_P separately.

---

## 2. Resp_F^conv: absolute value added

**Old (broken during mobilization):**
```
Resp_F^conv = (1−η)·[(β_i·Π − α_i·Π^δ)·F_i + (β_n·Π − α_n·Π^δ)·F_n]
```

**New (corrected):**
```
Resp_F^conv = (1−η)·|[(β_i·Π − α_i·Π^δ)·F_i + (β_n·Π − α_n·Π^δ)·F_n]|
```

**Why:** When mobilization dominates (Π small, α·Π^δ > β·Π), the bracketed expression goes negative. Without abs(), Resp_F^conv becomes negative — nonphysical negative respiration that breaks carbon conservation.

**What to check in code:** Wrap the bracketed sum in `abs()`. This affects:
- `S_O` (oxygen consumption): `S_O = −α_O·(Resp_B + Resp_F + Resp_F^conv)`
- `S_Fm` (mobile fungi source): `S_Fm = Γ_F − η·(...)·F_i − η·(...)·F_n − Resp_F^conv`
- CO₂ cumulative diagnostic

**Note:** The η-weighted transfer terms in S_Fi, S_Fn, S_Fm are UNCHANGED — they correctly handle both directions as written. Only the Resp_F^conv respiration diagnostic needs abs().

---

## 3. EPS Degradation: C → C_aq

**Old:**
```
R_{rec,E} = μ_E^max(T) · K_E/(K_E + C) · E · h_E
```

**New:**
```
R_{rec,E} = μ_E^max(T) · K_E/(K_E + C_aq) · E · h_E
```

where `C_aq = C / (θ + ρ_b · k_d)`.

**Why:** The inhibition logic is "scavenge EPS when dissolved substrate is scarce." The sorbed fraction C_eq is unavailable for this purpose. Using C_aq is consistent with microbial uptake (R_B, R_F) which also uses C_aq.

**What to check in code:** In the EPS degradation function (biology/eps.jl or equivalent), replace the raw state variable C with C_aq. You should already be computing C_aq for microbial uptake at each node — reuse that value.

---

## Summary of C / C_aq / C_eq usage (complete)

All three are derived from the single state variable C:
```
C      = state variable (total mobile carbon)           [μg/mm³]
C_aq   = C / (θ + ρ_b·k_d)                              [μg/mm³ water]
C_eq   = k_d · C_aq = k_d·C / (θ + ρ_b·k_d)            [μg/mm³ solid]
```

| Where used | Which concentration |
|-----------|-------------------|
| Diffusion PDE (state variable) | C |
| Bacterial uptake R_B (Monod) | C_aq |
| Fungal uptake R_F (Monod) | C_aq |
| Bacterial maintenance R_{B,b} | C_B (constant parameter, not C_aq) |
| EPS degradation inhibition | C_aq ← **CHANGED** |
| MAOC equilibrium M_eq (Langmuir-Freundlich) | C_eq |
| MAOC rate law J_M | operates on M_eq which uses C_eq |
| S_C coupling to MAOC | −J_M·(θ + ρ_b·k_d)/k_d |
