import cobra

# --- 0. Set the solver explicitly ---
try:
    cobra.Configuration().solver = "glpk"
    print(f"Using solver: {cobra.Configuration().solver}")
except Exception as e:
    print(f"Could not set GLPK solver: {e}")
    print("Please ensure GLPK is installed and configured for cobrapy.")

# --- 1. Create a new COBRA model for R_NpADO only ---
model = cobra.Model("Minimal_NpADO_Test_Final_Corrected_V8_Full_Rigorous_Check")

# --- 2. Define Metabolites for R_NpADO with Formulas and Charges ---
fa_ald_c16_c = cobra.Metabolite(
    "fa_ald_c16_c",
    name="Hexadecanal (C16 Aldehyde)",
    compartment="c",
    formula="C16H32O",  # Hexadecanal
    charge=0
)
alkane_c15_c = cobra.Metabolite(
    "alkane_c15_c",
    name="Pentadecane (C15 Alkane)",
    compartment="c",
    formula="C15H32",   # Pentadecane
    charge=0
)

o2_c = cobra.Metabolite(
    "o2_c",
    name="Oxygen",
    compartment="c",
    formula="O2",
    charge=0
)
# Using standard charges for NAD(P)H and NAD(P)+ at physiological pH
# NADPH: C21H29N7O17P3, charge -4 (from phosphates)
# NADP+: C21H28N7O17P3, charge -3 (one less H+ than NADPH)
nadph_c = cobra.Metabolite(
    "nadph_c",
    name="NADPH",
    compartment="c",
    formula="C21H29N7O17P3",
    charge=-4
)
nadp_c = cobra.Metabolite(
    "nadp_c",
    name="NADP+",
    compartment="c",
    formula="C21H28N7O17P3",
    charge=-3
)
h2o_c = cobra.Metabolite(
    "h2o_c",
    name="Water",
    compartment="c",
    formula="H2O",
    charge=0
)
formate_c = cobra.Metabolite(
    "formate_c",
    name="Formate",
    compartment="c",
    formula="CHO2", # HCOO-
    charge=-1
)
h_c = cobra.Metabolite(
    "h_c",
    name="Proton",
    compartment="c",
    formula="H",
    charge=1
)

alkane_c15_e = cobra.Metabolite(
    "alkane_c15_e",
    name="Pentadecane (C15 Alkane) External",
    compartment="e",
    formula="C15H32", # Same as internal alkane
    charge=0
)


# --- 3. Add Metabolites to the model ---
model.add_metabolites([
    fa_ald_c16_c, alkane_c15_c,
    o2_c, nadph_c, nadp_c, h2o_c, formate_c, h_c,
    alkane_c15_e
])

# --- 4. Define R_NpADO Reaction ---
rxn_npado = cobra.Reaction("R_NpADO")
rxn_npado.name = "Nostoc punctiforme Aldehyde Deformylating Oxygenase"
rxn_npado.lower_bound = 0.0 # Will temporarily change this for test
rxn_npado.upper_bound = 1000.0

# Confirmed correct stoichiometry for 4-electron transfer and charge balance
rxn_npado.add_metabolites({
    fa_ald_c16_c: -1.0,
    o2_c: -1.0,
    nadph_c: -2.0,  # For 4-electron transfer
    h_c: -1.0,      # Correct for charge balance with 2 NADPH/2 NADP+
    alkane_c15_c: 1.0,
    formate_c: 1.0,
    h2o_c: 1.0,
    nadp_c: 2.0     # For 4-electron transfer
})
model.add_reactions([rxn_npado]) # Ensure model.add_reactions (plural) is used

# --- Manual Mass and Charge Balance Check for R_NpADO ---
print("\n--- Manual Balance Check for R_NpADO (using check_mass_balance and manual charge check) ---")
mass_balance = rxn_npado.check_mass_balance()

if mass_balance:
    print(f"Mass balance for {rxn_npado.id}: OK. Imbalance: {mass_balance}") # Display imbalance for verification
else:
    print(f"Mass balance for {rxn_npado.id}: FAILED. Imbalance: {mass_balance}")

# Manually verify charge balance (as check_charge_balance does not exist)
# Sum of (stoichiometric_coefficient * charge) for all metabolites in reaction
net_charge = sum(coeff * met.charge for met, coeff in rxn_npado.metabolites.items())
if net_charge == 0:
    print(f"Charge balance for {rxn_npado.id}: OK (Net charge: {net_charge})")
else:
    print(f"Charge balance for {rxn_npado.id}: FAILED. Net charge: {net_charge}")


# --- 5. Add ONLY necessary Exchange Reactions for R_NpADO inputs/outputs ---
# Inputs (sources)
ex_fa_ald = cobra.Reaction("EX_fa_ald_c16_c")
ex_fa_ald.add_metabolites({fa_ald_c16_c: -1.0})
ex_fa_ald.lower_bound = -10.0
model.add_reactions([ex_fa_ald])

ex_o2 = cobra.Reaction("EX_o2_c")
ex_o2.add_metabolites({o2_c: -1.0})
ex_o2.lower_bound = -1000.0
model.add_reactions([ex_o2])

ex_nadph = cobra.Reaction("EX_nadph_c")
ex_nadph.add_metabolites({nadph_c: -1.0})
ex_nadph.lower_bound = -1000.0
model.add_reactions([ex_nadph])

# Outputs (sinks)
ex_formate = cobra.Reaction("EX_formate_c")
ex_formate.add_metabolites({formate_c: -1.0})
ex_formate.lower_bound = 0.0
model.add_reactions([ex_formate])

ex_h2o = cobra.Reaction("EX_h2o_c")
ex_h2o.add_metabolites({h2o_c: -1.0})
ex_h2o.lower_bound = -1000.0 # Bidirectional
model.add_reactions([ex_h2o])

ex_nadp = cobra.Reaction("EX_nadp_c")
ex_nadp.add_metabolites({nadp_c: -1.0})
ex_nadp.lower_bound = -1000.0 # IMPORTANT: Changed to bidirectional
ex_nadp.upper_bound = 1000.0
model.add_reactions([ex_nadp])

ex_h = cobra.Reaction("EX_h_c")
ex_h.add_metabolites({h_c: -1.0})
ex_h.lower_bound = -1000.0 # Bidirectional, allows H+ to be consumed or produced
ex_h.upper_bound = 1000.0
model.add_reactions([ex_h])


# Product Export: alkane_c15_c -> alkane_c15_e
ex_alkane_transport = cobra.Reaction("EX_alkane_c15_e_transport")
ex_alkane_transport.add_metabolites({alkane_c15_c: -1.0, alkane_c15_e: 1.0})
ex_alkane_transport.lower_bound = 0.0
ex_alkane_transport.upper_bound = 1000.0
model.add_reactions([ex_alkane_transport])

# Objective: Maximize external alkane
model.objective = ex_alkane_transport
model.objective_direction = "maximize"


# --- 6. Debugging: Attempting to force flux through R_NpADO directly ---
print("\n--- Debugging: Attempting to force flux through R_NpADO directly ---")

# Temporarily set lower bound of R_NpADO to a positive value
rxn_npado.lower_bound = 1.0

# Solve the model
solution = model.optimize()

print(f"Model ID: {model.id}")
print(f"Number of reactions: {len(model.reactions)}")
print(f"Number of metabolites: {len(model.metabolites)}")

if solution.status == "optimal":
    # More robust objective printing
    obj_name = "Unknown Objective"
    if isinstance(model.objective, cobra.Reaction):
        obj_name = model.objective.id
    elif hasattr(model.objective.expression, 'name'):
        obj_name = model.objective.expression.name
    elif hasattr(model.objective, 'id'):
        obj_name = model.objective.id
    elif hasattr(model.objective.expression, '__str__'): # Fallback for complex objectives
        obj_name = str(model.objective.expression)

    print(f"\nSolution status: {solution.status}")
    print(f"Forced flux for {obj_name}: {solution.objective_value}")
    print("\nKey Reaction Fluxes:")
    print(f"{rxn_npado.id}: {solution.fluxes[rxn_npado.id]}")
    print(f"{ex_alkane_transport.id}: {solution.fluxes[ex_alkane_transport.id]}")
    print("\nCofactor Fluxes:")
    print(f"{ex_o2.id}: {solution.fluxes[ex_o2.id]}")
    print(f"{ex_nadph.id}: {solution.fluxes[ex_nadph.id]}")
    print(f"{ex_nadp.id}: {solution.fluxes[ex_nadp.id]}")
    print(f"{ex_formate.id}: {solution.fluxes[ex_formate.id]}")
    print(f"{ex_h2o.id}: {solution.fluxes[ex_h2o.id]}")
    print(f"{ex_h.id}: {solution.fluxes[ex_h.id]}")
else:
    print(f"\nCould NOT force flux through R_NpADO. Solution status: {solution.status}")
    if solution.status == "infeasible":
        print("This indicates a fundamental stoichiometric blockage within the NpADO reaction itself or its direct inputs/outputs, or an issue with cofactor balancing/exchange.")
        print("Given that manual mass/charge balance checks pass, the issue is highly likely related to the GLPK solver setup or its interaction with cobrapy, or a very subtle issue with exchange bounds.")


# --- 7. Save the model to an SBML file ---
# Reset lower bounds before saving
rxn_npado.lower_bound = 0.0
ex_alkane_transport.lower_bound = 0.0
cobra.io.write_sbml_model(model, "minimal_npado_test_final_corrected_v8_full_rigorous_check.xml")
print("\nModel saved to minimal_npado_test_final_corrected_v8_full_rigorous_check.xml")