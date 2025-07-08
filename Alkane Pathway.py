from optlang import glpk_interface

# --- 1. Create a solver object ---
solver = glpk_interface.Model()
print(f"Using solver: {glpk_interface.Model.__module__}")

# --- 2. Define Flux Variables for each reaction and add them to the solver ---

# Core Pathway Reactions
v_npado = glpk_interface.Variable("R_NpADO", lb=0.0, ub=1000.0)
solver.add(v_npado)

v_nadph_regen = glpk_interface.Variable("R_NADPH_regeneration", lb=0.0, ub=1000.0)
solver.add(v_nadph_regen)

# ATP Maintenance reaction (requires a minimum flux for cellular activity)
v_atp_maintenance = glpk_interface.Variable("R_ATP_Maintenance", lb=1.0, ub=1000.0) # Enforce a minimum ATP demand
solver.add(v_atp_maintenance)

# Simplified ATP Generation reaction
v_atp_generation = glpk_interface.Variable("R_ATP_Generation", lb=0.0, ub=1000.0)
solver.add(v_atp_generation)


# Exchange Reactions (inputs and outputs to/from the system)
wide_bound = 999999.0

# Precursor uptake: Specific fixed input for fa_ald_c16_c
v_ex_fa_ald = glpk_interface.Variable("EX_fa_ald_c16_c", lb=-10.0, ub=0.0) # Allows uptake up to 10 units.
solver.add(v_ex_fa_ald)

# Unlimited inputs (negative lower bound for uptake, positive upper bound for potential export)
v_ex_o2 = glpk_interface.Variable("EX_o2_c", lb=-wide_bound, ub=wide_bound)
solver.add(v_ex_o2)

# Bounds for NADPH and NADP exchange remain 0.0 to force internal regeneration
v_ex_nadph = glpk_interface.Variable("EX_nadph_c", lb=0.0, ub=0.0)
solver.add(v_ex_nadph)

v_ex_nadp = glpk_interface.Variable("EX_nadp_c", lb=0.0, ub=0.0)
solver.add(v_ex_nadp)

v_ex_h2o = glpk_interface.Variable("EX_h2o_c", lb=-wide_bound, ub=wide_bound)
solver.add(v_ex_h2o)

# **CRITICAL FIX HERE**: Removed the extra '.interface'
v_ex_formate = glpk_interface.Variable("EX_formate_c", lb=-wide_bound, ub=wide_bound)
solver.add(v_ex_formate)

# Unlimited input for H+
v_ex_h = glpk_interface.Variable("EX_h_c", lb=-wide_bound, ub=wide_bound)
solver.add(v_ex_h)

# Unlimited input for electrons
v_ex_e = glpk_interface.Variable("EX_e_c", lb=-wide_bound, ub=wide_bound)
solver.add(v_ex_e)

# Product export
v_ex_alkane_transport = glpk_interface.Variable("EX_alkane_c15_e_transport", lb=0.0, ub=wide_bound)
solver.add(v_ex_alkane_transport)


# --- 3. Define Metabolite Mass Balance Constraints ---

# Stoichiometric coefficients for existing reactions
# R_NpADO: fa_ald_c16_c + O2 + 2 NADPH + H+ + 4 e- -> alkane_c15_c + formate + H2O + 2 NADP+
npado_coeffs = {
    "fa_ald_c16_c": -1.0, "o2_c": -1.0, "nadph_c": -2.0, "h_c": -1.0, "e_c": -4.0,
    "alkane_c15_c": 1.0, "formate_c": 1.0, "h2o_c": 1.0, "nadp_c": 2.0
}

# R_NADPH_regeneration: 2 NADP+ + 2 H+ -> 2 NADPH
nadph_regen_coeffs = {
    "nadp_c": -2.0, "h_c": -2.0, "nadph_c": 2.0
}

# New: Stoichiometric coefficients for new reactions
# R_ATP_Maintenance: ATP_c + H2O_c -> ADP_c + PI_c + H_c
atp_maint_coeffs = {
    "atp_c": -1.0, "h2o_c": -1.0, "adp_c": 1.0, "pi_c": 1.0, "h_c": 1.0
}

# R_ATP_Generation: 4 e_c + 4 H_c + O2_c -> ATP_c + 2 H2O_c
atp_gen_coeffs = {
    "e_c": -4.0, "h_c": -4.0, "o2_c": -1.0, "atp_c": 1.0, "h2o_c": 2.0
}


constraints = []

# Metabolic Balances for all metabolites
# Note: For metabolites involved in new reactions, their balance equations are updated.
# New metabolites: ATP_c, ADP_c, PI_c

# balance_fa_ald_c16_c
constraints.append(glpk_interface.Constraint(
    v_npado * npado_coeffs["fa_ald_c16_c"] + v_ex_fa_ald * 1.0,
    lb=0, ub=0, name="balance_fa_ald_c16_c"
))

# balance_o2_c (now involved in NpADO, Exchange, and ATP_Generation)
constraints.append(glpk_interface.Constraint(
    v_npado * npado_coeffs["o2_c"] +
    v_atp_generation * atp_gen_coeffs["o2_c"] +
    v_ex_o2 * 1.0,
    lb=0, ub=0, name="balance_o2_c"
))

# balance_nadph_c (NpADO consumes, NADPH_regen produces, no exchange)
constraints.append(glpk_interface.Constraint(
    v_npado * npado_coeffs["nadph_c"] +
    v_nadph_regen * nadph_regen_coeffs["nadph_c"] +
    v_ex_nadph * 1.0,
    lb=0, ub=0, name="balance_nadph_c"
))

# balance_nadp_c (NpADO produces, NADPH_regen consumes, no exchange)
constraints.append(glpk_interface.Constraint(
    v_npado * npado_coeffs["nadp_c"] +
    v_nadph_regen * nadph_regen_coeffs["nadp_c"] +
    v_ex_nadp * 1.0,
    lb=0, ub=0, name="balance_nadp_c"
))

# balance_h2o_c (NpADO produces, ATP_Maintenance consumes, ATP_Generation produces, Exchange)
constraints.append(glpk_interface.Constraint(
    v_npado * npado_coeffs["h2o_c"] +
    v_atp_maintenance * atp_maint_coeffs["h2o_c"] +
    v_atp_generation * atp_gen_coeffs["h2o_c"] +
    v_ex_h2o * 1.0,
    lb=0, ub=0, name="balance_h2o_c"
))

# balance_formate_c
constraints.append(glpk_interface.Constraint(
    v_npado * npado_coeffs["formate_c"] + v_ex_formate * 1.0,
    lb=0, ub=0, name="balance_formate_c"
))

# balance_h_c (NpADO consumes, NADPH_regen consumes, ATP_Maintenance produces, ATP_Generation consumes, Exchange)
constraints.append(glpk_interface.Constraint(
    v_npado * npado_coeffs["h_c"] +
    v_nadph_regen * nadph_regen_coeffs["h_c"] +
    v_atp_maintenance * atp_maint_coeffs["h_c"] +
    v_atp_generation * atp_gen_coeffs["h_c"] +
    v_ex_h * 1.0,
    lb=0, ub=0, name="balance_h_c"
))

# balance_e_c (NpADO consumes, ATP_Generation consumes, Exchange)
constraints.append(glpk_interface.Constraint(
    v_npado * npado_coeffs["e_c"] +
    v_atp_generation * atp_gen_coeffs["e_c"] +
    v_ex_e * 1.0,
    lb=0, ub=0, name="balance_e_c"
))

# balance_alkane_c15_c
constraints.append(glpk_interface.Constraint(
    v_npado * npado_coeffs["alkane_c15_c"] + v_ex_alkane_transport * -1.0,
    lb=0, ub=0, name="balance_alkane_c15_c"
))

# NEW: balance_atp_c (ATP_Maintenance consumes, ATP_Generation produces)
constraints.append(glpk_interface.Constraint(
    v_atp_maintenance * atp_maint_coeffs["atp_c"] +
    v_atp_generation * atp_gen_coeffs["atp_c"],
    lb=0, ub=0, name="balance_atp_c"
))

# NEW: balance_adp_c (ATP_Maintenance produces)
constraints.append(glpk_interface.Constraint(
    v_atp_maintenance * atp_maint_coeffs["adp_c"],
    lb=0, ub=0, name="balance_adp_c"
))

# NEW: balance_pi_c (ATP_Maintenance produces)
constraints.append(glpk_interface.Constraint(
    v_atp_maintenance * atp_maint_coeffs["pi_c"],
    lb=0, ub=0, name="balance_pi_c"
))

solver.add(constraints)


# --- 4. Define the Objective Function ---
# Maximize the flux through EX_alkane_c15_e_transport
objective_expression = v_ex_alkane_transport
objective = glpk_interface.Objective(objective_expression, direction="max")
solver.objective = objective

# --- 5. Solve the problem ---
status = solver.optimize()

print(f"\nSolution status: {status}")
print(f"Objective Value (EX_alkane_c15_e_transport): {solver.objective.value}")

# --- 6. Print Key Reaction Fluxes ---
print("\nKey Reaction Fluxes:")
print(f"R_NpADO: {solver.primal_values.get('R_NpADO', 'N/A')}")
print(f"R_NADPH_regeneration: {solver.primal_values.get('R_NADPH_regeneration', 'N/A')}")
print(f"R_ATP_Maintenance: {solver.primal_values.get('R_ATP_Maintenance', 'N/A')}")
print(f"R_ATP_Generation: {solver.primal_values.get('R_ATP_Generation', 'N/A')}")
print(f"EX_alkane_c15_e_transport: {solver.primal_values.get('EX_alkane_c15_e_transport', 'N/A')}")


print("\nCofactor/Exchange Fluxes:")
print(f"EX_fa_ald_c16_c: {solver.primal_values.get('EX_fa_ald_c16_c', 'N/A')}")
print(f"EX_o2_c: {solver.primal_values.get('EX_o2_c', 'N/A')}")
print(f"EX_nadph_c: {solver.primal_values.get('EX_nadph_c', 'N/A')}")
print(f"EX_nadp_c: {solver.primal_values.get('EX_nadp_c', 'N/A')}")
print(f"EX_formate_c: {solver.primal_values.get('EX_formate_c', 'N/A')}")
print(f"EX_h2o_c: {solver.primal_values.get('EX_h2o_c', 'N/A')}")
print(f"EX_h_c: {solver.primal_values.get('EX_h_c', 'N/A')}")
print(f"EX_e_c: {solver.primal_values.get('EX_e_c', 'N/A')}")

# --- 7. Additional Debugging Check: Force NpADO flux and see if it's possible ---
# This part is just for diagnosis and won't affect the main optimization result.
# solver.problem.objective.expression = v_npado
# solver.problem.objective.direction = "max"
# solver.variables["R_NpADO"].lb = 1.0 # Try to force at least 1 unit of NpADO flux
# test_status = solver.optimize()
# print(f"\n--- Debugging: Attempting to force flux through R_NpADO directly ---")
# print(f"Test status (forcing R_NpADO=1.0): {test_status}")
# if test_status == "optimal":
#     print(f"R_NpADO flux when forced: {solver.primal_values.get('R_NpADO', 'N/A')}")
#     print(f"EX_h_c flux when NpADO forced: {solver.primal_values.get('EX_h_c', 'N/A')}")
#     print(f"EX_e_c flux when NpADO forced: {solver.primal_values.get('EX_e_c', 'N/A')}")
#     print(f"EX_o2_c flux when NpADO forced: {solver.primal_values.get('EX_o2_c', 'N/A')}")
# else:
#     print(f"Could NOT force flux through R_NpADO. This indicates a fundamental stoichiometric or cofactor balancing issue.")
# # Reset NpADO lower bound for the main objective if you were to run the original objective again
# solver.variables["R_NpADO"].lb = 0.0
# solver.problem.objective.expression = objective_expression
# solver.problem.objective.direction = "max"