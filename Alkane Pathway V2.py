import cobra
from cobra.core import Reaction, Metabolite
from cobra import Model
from libchebipy import ChebiEntity # Import libChEBIpy

# Custom Exception for clearer error handling
class MetaboliteNotFoundException(Exception):
    pass

# --- Helper function to fetch metabolite info from ChEBI (more robust) ---
def get_chebi_info_robust(name):
    """
    Fetches ChEBI ID, formula, and charge for a given metabolite name.
    Returns a dictionary {'id': chebi_id, 'formula': formula, 'charge': charge}
    Raises MetaboliteNotFoundException if not found or error occurs.
    """
    try:
        results = ChebiEntity.search_entity(name, "CHEBI_NAME", False)
        
        # Filter for exact name match (case-insensitive for robustness)
        exact_results = [r for r in results if ChebiEntity(r.chebi_id).get_name().lower() == name.lower()]

        if exact_results:
            chebi_id_found = exact_results[0].chebi_id
        elif results: # If no exact match, but some results exist, take the first
            print(f"Warning: No exact ChEBI match for '{name}'. Using first search result: {results[0].chebi_id}.")
            chebi_id_found = results[0].chebi_id
        else:
            raise MetaboliteNotFoundException(f"Could not find ChEBI ID for '{name}'.")

        entity = ChebiEntity(chebi_id_found)
        formula = entity.get_formula() if entity.get_formula() else ""
        charge = entity.get_charge() if entity.get_charge() is not None else 0 # Ensure charge is int or 0

        print(f"Found '{name}': CHEBI:{chebi_id_found}, Formula: {formula}, Charge: {charge}")
        return {'id': chebi_id_found, 'formula': formula, 'charge': charge}

    except Exception as e:
        print(f"Warning: Error fetching ChEBI info for '{name}': {e}. Falling back to hardcoded defaults for {name}.")
        # No re-raise here, we explicitly want to fall back to hardcoded if ChEBI lookup fails.


# --- Main COBRApy model construction ---
# --- 1. Create a new COBRA model ---
model = Model("Combined_Alkane_FAEE_Pathway_Model") # Changed model name to reflect combined pathways
model.solver = "glpk" # Explicitly set the solver

# Define compartment
model.compartments = {"c": "cytosol"}

# --- 2. Define Metabolites ---
# Define metabolites with their desired model IDs and ChEBI search names.
# Includes hardcoded fallbacks as before, and new ones for FAEE pathway.
metabolite_definitions = {
    # Alkane pathway metabolites
    "fa_ald_c16_c": {"name": "hexadecanal", "formula": "C16H32O", "charge": 0},
    "o2_c": {"name": "oxygen", "formula": "O2", "charge": 0},
    "nadph_c": {"name": "NADPH", "formula": "C21H29N7O17P3", "charge": -4},
    "h_c": {"name": "proton", "formula": "H", "charge": 1},
    "e_c": {"name": "electron", "formula": "", "charge": -1},
    "alkane_c15_c": {"name": "pentadecane", "formula": "C15H32", "charge": 0},
    "formate_c": {"name": "formate", "formula": "CHO2", "charge": -1},
    "h2o_c": {"name": "water", "formula": "H2O", "charge": 0},
    "nadp_c": {"name": "NADP(+)", "formula": "C21H26N7O17P3", "charge": -3},
    "atp_c": {"name": "ATP", "formula": "C10H16N5O13P3", "charge": -4},
    "adp_c": {"name": "ADP", "formula": "C10H15N5O10P2", "charge": -3},
    "pi_c": {"name": "phosphate", "formula": "H2O4P", "charge": -2},

    # New FAEE pathway metabolites
    "ethanol_c": {"name": "ethanol", "formula": "C2H6O", "charge": 0},
    "fatty_acyl_coa_c": {"name": "hexadecanoyl-CoA", "formula": "C27H44N7O17P3S", "charge": -5}, # Generic C16 acyl-CoA
    "faee_c": {"name": "ethyl hexadecanoate", "formula": "C18H36O2", "charge": 0}, # C16 FAEE (Hexadecanoic acid ethyl ester)
    "coa_c": {"name": "Coenzyme A", "formula": "C21H36N7O16P3S", "charge": -4} # ChEBI: 15346. Formula and charge from ChEBI are often without protonation state.
}

# Create Metabolite objects.
metabolites = []
for model_id, def_info in metabolite_definitions.items():
    met_name = def_info["name"]
    fetched_formula = def_info["formula"]
    fetched_charge = def_info["charge"]

    try:
        # Only try ChEBI for metabolites that are NOT explicitly hardcoded (like some aldehydes/alkanes)
        # or those where ChEBI lookups were consistently failing or hardcoded values are desired.
        if model_id not in ["fa_ald_c16_c", "alkane_c15_c", "e_c", "fatty_acyl_coa_c", "faee_c", "coa_c"]: # these were already hardcoded or now being explicitly hardcoded
            chebi_info = get_chebi_info_robust(met_name)
            fetched_formula = chebi_info["formula"]
            fetched_charge = chebi_info["charge"]
        else:
             print(f"Using hardcoded definition for {model_id} ({met_name}).")

    except MetaboliteNotFoundException as e:
        print(f"Warning: {e}. Falling back to hardcoded defaults for {model_id} ({met_name}).")
    except Exception as e:
        print(f"Warning: Error fetching ChEBI info for '{met_name}': {e}. Falling back to hardcoded defaults for {model_id} ({met_name}).")


    met = Metabolite(id=model_id,
                     name=met_name,
                     formula=fetched_formula,
                     charge=fetched_charge,
                     compartment="c")
    metabolites.append(met)

model.add_metabolites(metabolites)

# Convenient access to metabolites by their model ID
# Alkane pathway metabolites
fa_ald_c16_c = model.metabolites.get_by_id("fa_ald_c16_c")
o2_c = model.metabolites.get_by_id("o2_c")
nadph_c = model.metabolites.get_by_id("nadph_c")
h_c = model.metabolites.get_by_id("h_c")
e_c = model.metabolites.get_by_id("e_c")
alkane_c15_c = model.metabolites.get_by_id("alkane_c15_c")
formate_c = model.metabolites.get_by_id("formate_c")
h2o_c = model.metabolites.get_by_id("h2o_c")
nadp_c = model.metabolites.get_by_id("nadp_c")
atp_c = model.metabolites.get_by_id("atp_c")
adp_c = model.metabolites.get_by_id("adp_c")
pi_c = model.metabolites.get_by_id("pi_c")

# New FAEE pathway metabolites
ethanol_c = model.metabolites.get_by_id("ethanol_c")
fatty_acyl_coa_c = model.metabolites.get_by_id("fatty_acyl_coa_c")
faee_c = model.metabolites.get_by_id("faee_c")
coa_c = model.metabolites.get_by_id("coa_c")


# --- 3. Define Reactions ---

# Alkane Pathway Reactions (retained)
# R_NpADO: fa_ald_c16_c + O2 + 2 NADPH -> alkane_c15_c + formate_c + h2o_c + 2 nadp_c + 3 h_c + 4 e_c (Original balanced version)
rxn_npado = Reaction("R_NpADO")
rxn_npado.name = "NpADO Reaction"
rxn_npado.lower_bound = 0.0
rxn_npado.upper_bound = 1000.0
rxn_npado.add_metabolites({
    fa_ald_c16_c: -1, o2_c: -1, nadph_c: -2,
    alkane_c15_c: 1, formate_c: 1, h2o_c: 1, nadp_c: 2, h_c: 3, e_c: 4
})
model.add_reactions([rxn_npado])


# R_NADPH_regeneration: 2 NADP+ + 6 H+ + 8 e- -> 2 NADPH (CORRECTED stoichiometry and re-added)
rxn_nadph_regen = Reaction("R_NADPH_regeneration")
rxn_nadph_regen.name = "NADPH Regeneration"
rxn_nadph_regen.lower_bound = 0.0
rxn_nadph_regen.upper_bound = 1000.0
rxn_nadph_regen.add_metabolites({
    nadp_c: -2, h_c: -6, e_c: -8, 
    nadph_c: 2 
})
model.add_reactions([rxn_nadph_regen])


# R_ATP_Maintenance: ATP_c + H2O_c -> ADP_c + PI_c + H_c
rxn_atp_maintenance = Reaction("R_ATP_Maintenance")
rxn_atp_maintenance.name = "ATP Maintenance"
rxn_atp_maintenance.lower_bound = 0.0 # Temporarily set to 0.0
rxn_atp_maintenance.upper_bound = 1000.0
rxn_atp_maintenance.add_metabolites({
    atp_c: -1, h2o_c: -1,
    adp_c: 1, pi_c: 1, h_c: 1
})
model.add_reactions([rxn_atp_maintenance])


# New FAEE Pathway Reactions
# R_FAEE_synth: fatty_acyl_coa_c + ethanol_c -> faee_c + coa_c
rxn_faee_synth = Reaction("R_FAEE_synth")
rxn_faee_synth.name = "Fatty Acid Ethyl Ester Synthesis"
rxn_faee_synth.lower_bound = 0.0
rxn_faee_synth.upper_bound = 1000.0
rxn_faee_synth.add_metabolites({
    fatty_acyl_coa_c: -1, ethanol_c: -1,
    faee_c: 1, coa_c: 1
})
model.add_reactions([rxn_faee_synth])


# --- 4. Define Exchange Reactions ---
wide_bound = 999999.0

# Alkane pathway related exchange/demand (retained)
ex_fa_ald_c16_c = model.add_boundary(fa_ald_c16_c, type="demand", reaction_id="EX_fa_ald_c16_c", lb=-10.0, ub=0.0)
ex_alkane_c15_e_transport = model.add_boundary(alkane_c15_c, type="demand", reaction_id="EX_alkane_c15_e_transport", lb=0.0, ub=wide_bound) # Lower bound 0.0 for debugging

# General inputs/outputs (retained)
model.add_boundary(o2_c, type="exchange", reaction_id="EX_o2_c", lb=-wide_bound, ub=0.0)
model.add_boundary(h_c, type="exchange", reaction_id="EX_h_c", lb=-wide_bound, ub=0.0)
model.add_boundary(e_c, type="exchange", reaction_id="EX_e_c", lb=-wide_bound, ub=0.0)
model.add_boundary(h2o_c, type="exchange", reaction_id="EX_h2o_c", lb=-wide_bound, ub=wide_bound)
model.add_boundary(formate_c, type="exchange", reaction_id="EX_formate_c", lb=0.0, ub=wide_bound)

# Cofactor demand reactions (retained)
model.add_boundary(nadph_c, type="demand", reaction_id="DM_nadph_c", lb=0.0, ub=0.0)
model.add_boundary(nadp_c, type="demand", reaction_id="DM_nadp_c", lb=0.0, ub=0.0)
model.add_boundary(atp_c, type="demand", reaction_id="DM_atp_c", lb=0.0, ub=0.0)
model.add_boundary(adp_c, type="demand", reaction_id="DM_adp_c", lb=0.0, ub=0.0)
model.add_boundary(pi_c, type="demand", reaction_id="DM_pi_c", lb=0.0, ub=0.0)


# New FAEE pathway related exchange/demand reactions
ex_ethanol_c = model.add_boundary(ethanol_c, type="demand", reaction_id="EX_ethanol_c", lb=-10.0, ub=0.0) # Allow uptake of ethanol
ex_fatty_acyl_coa_c = model.add_boundary(fatty_acyl_coa_c, type="demand", reaction_id="EX_fatty_acyl_coa_c", lb=-10.0, ub=0.0) # Allow uptake of fatty acyl-CoA
ex_faee_c_export = model.add_boundary(faee_c, type="demand", reaction_id="EX_faee_c_export", lb=0.0, ub=wide_bound) # Export of FAEE, can be objective
model.add_boundary(coa_c, type="demand", reaction_id="DM_coa_c", lb=0.0, ub=0.0) # Demand for CoA for balance


# --- 5. Set the Objective Function ---
# Initially, we can set the objective to FAEE production to test it.
# We can switch this later to alkane production or a combined objective.
model.objective = ex_faee_c_export # Maximize FAEE export

# --- 6. Solve the problem ---
solution = model.optimize()

print(f"\nSolution status: {solution.status}")
print(f"Objective Value ({model.objective.id}): {solution.objective_value}") # Print the ID of the objective

# --- 7. Print Key Reaction Fluxes ---
print("\nKey Reaction Fluxes (Alkane Pathway):")
print(f"R_NpADO: {solution.fluxes.get('R_NpADO', 'N/A')}")
print(f"R_NADPH_regeneration: {solution.fluxes.get('R_NADPH_regeneration', 'N/A')}")
print(f"R_ATP_Maintenance: {solution.fluxes.get('R_ATP_Maintenance', 'N/A')}")
print(f"EX_alkane_c15_e_transport: {solution.fluxes.get('EX_alkane_c15_e_transport', 'N/A')}")

print("\nKey Reaction Fluxes (FAEE Pathway):")
print(f"R_FAEE_synth: {solution.fluxes.get('R_FAEE_synth', 'N/A')}")
print(f"EX_faee_c_export: {solution.fluxes.get('EX_faee_c_export', 'N/A')}")

print("\nCofactor/Exchange Fluxes:")
print(f"EX_fa_ald_c16_c: {solution.fluxes.get('EX_fa_ald_c16_c', 'N/A')}")
print(f"EX_o2_c: {solution.fluxes.get('EX_o2_c', 'N/A')}")
print(f"EX_ethanol_c: {solution.fluxes.get('EX_ethanol_c', 'N/A')}")
print(f"EX_fatty_acyl_coa_c: {solution.fluxes.get('EX_fatty_acyl_coa_c', 'N/A')}")
print(f"DM_nadph_c: {solution.fluxes.get('DM_nadph_c', 'N/A')}")
print(f"DM_nadp_c: {solution.fluxes.get('DM_nadp_c', 'N/A')}")
print(f"DM_atp_c: {solution.fluxes.get('DM_atp_c', 'N/A')}")
print(f"DM_adp_c: {solution.fluxes.get('DM_adp_c', 'N/A')}")
print(f"DM_pi_c: {solution.fluxes.get('DM_pi_c', 'N/A')}")
print(f"DM_coa_c: {solution.fluxes.get('DM_coa_c', 'N/A')}")
print(f"EX_formate_c: {solution.fluxes.get('EX_formate_c', 'N/A')}")
print(f"EX_h2o_c: {solution.fluxes.get('EX_h2o_c', 'N/A')}")
print(f"EX_h_c: {solution.fluxes.get('EX_h_c', 'N/A')}")
print(f"EX_e_c: {solution.fluxes.get('EX_e_c', 'N/A')}")


# --- 8. COBRApy Mass and Charge and Elemental Balance Checks ---
print("\n--- COBRApy Mass and Charge and Elemental Balance Checks ---")
for rxn in model.reactions:
    mass_balance = rxn.check_mass_balance()
    if mass_balance:
        print(f"Mass balance for {rxn.id}: FAILED. Imbalance: {mass_balance}")
    else:
        print(f"Mass balance for {rxn.id}: OK")