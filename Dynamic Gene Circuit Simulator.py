import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# --- Parameters and Initial Conditions ---
# All parameters are stored in a dictionary for easy access and modification.
# You can adjust these values to see how they affect the simulation.

params = {
    # FAEE Pathway Parameters
    'k_LacI_prod': 0.1,      # Production rate of LacI
    'd_LacI': 0.01,          # Degradation rate of LacI
    'k_Allo_prod': 0.5,      # Production rate of Allolactose from lactose
    'd_Allo': 0.05,          # Degradation rate of Allolactose
    'k_bind': 0.1,           # Binding rate of LacI to Allolactose
    'k_unbind': 0.01,        # Unbinding rate of LacI from Allolactose
    'k_mRNA_prod_AEAT': 1,   # Production rate of AEAT mRNA
    'd_mRNA_AEAT': 0.1,      # Degradation rate of AEAT mRNA
    'k_protein_prod_AEAT': 0.5, # Production rate of AEAT protein
    'd_protein_AEAT': 0.05,  # Degradation rate of AEAT protein
    'k_FAEE_prod': 0.5,      # Production rate of FAEE from AEAT
    'd_FAEE': 0.005,         # Degradation rate of FAEE
    'K_lac': 5,              # Dissociation constant of LacI from promoter
    'n_lac': 2,              # Hill coefficient for LacI repression
    'initial_lactose': 100,  # Initial concentration of lactose feedstock

    # Alkane Pathway Parameters
    'k_FadR_prod': 0.1,      # Production rate of FadR
    'd_FadR': 0.01,          # Degradation rate of FadR
    'k_FA_prod': 0.5,        # Production rate of Fatty Acids from triglycerides
    'd_FA': 0.05,            # Degradation rate of Fatty Acids
    'k_bind_FA': 0.1,        # Binding rate of FadR to Fatty Acids
    'k_unbind_FA': 0.01,     # Unbinding rate of FadR from Fatty Acids
    'k_mRNA_prod_NpADO': 1,  # Production rate of NpADO mRNA
    'd_mRNA_NpADO': 0.1,     # Degradation rate of NpADO mRNA
    'k_protein_prod_NpADO': 0.5, # Production rate of NpADO protein
    'd_protein_NpADO': 0.05, # Degradation rate of NpADO protein
    'k_Alkane_prod': 0.5,    # Production rate of Alkane from NpADO
    'd_Alkane': 0.005,       # Degradation rate of Alkane
    'K_FadR': 5,             # Dissociation constant of FadR from promoter
    'n_FadR': 2,             # Hill coefficient for FadR repression
    'initial_triglycerides': 100 # Initial concentration of triglycerides
}

# Time span for the simulation
t_span = (0, 100)
t_eval = np.linspace(t_span[0], t_span[1], 1000)

# --- ODE Functions for each pathway ---

def faee_odes(t, y, p):
    """
    Defines the system of ordinary differential equations for the FAEE pathway.
    y = [LacI, Allo, Complex, mRNA_AEAT, AEAT, FAEE]
    p = parameters dictionary
    """
    LacI, Allo, Complex, mRNA_AEAT, AEAT, FAEE = y

    # Free LacI concentration
    free_LacI = LacI - Complex

    # Hill function for gene repression
    promoter_term = (p['K_lac']**p['n_lac']) / (p['K_lac']**p['n_lac'] + free_LacI**p['n_lac'])

    # ODEs
    d_LacI = p['k_LacI_prod'] - p['d_LacI'] * LacI
    d_Allo = p['k_Allo_prod'] * p['initial_lactose'] - p['d_Allo'] * Allo
    d_Complex = p['k_bind'] * free_LacI * Allo - p['k_unbind'] * Complex
    d_mRNA_AEAT = p['k_mRNA_prod_AEAT'] * promoter_term - p['d_mRNA_AEAT'] * mRNA_AEAT
    d_AEAT = p['k_protein_prod_AEAT'] * mRNA_AEAT - p['d_protein_AEAT'] * AEAT
    d_FAEE = p['k_FAEE_prod'] * AEAT - p['d_FAEE'] * FAEE
    
    return [d_LacI, d_Allo, d_Complex, d_mRNA_AEAT, d_AEAT, d_FAEE]

def alkane_odes(t, y, p):
    """
    Defines the system of ordinary differential equations for the Alkane pathway.
    y = [FadR, FA, Complex_FA, mRNA_NpADO, NpADO, Alkane]
    p = parameters dictionary
    """
    FadR, FA, Complex_FA, mRNA_NpADO, NpADO, Alkane = y

    # Free FadR concentration
    free_FadR = FadR - Complex_FA

    # Hill function for gene repression
    promoter_term = (p['K_FadR']**p['n_FadR']) / (p['K_FadR']**p['n_FadR'] + free_FadR**p['n_FadR'])

    # ODEs
    d_FadR = p['k_FadR_prod'] - p['d_FadR'] * FadR
    d_FA = p['k_FA_prod'] * p['initial_triglycerides'] - p['d_FA'] * FA
    d_Complex_FA = p['k_bind_FA'] * free_FadR * FA - p['k_unbind_FA'] * Complex_FA
    d_mRNA_NpADO = p['k_mRNA_prod_NpADO'] * promoter_term - p['d_mRNA_NpADO'] * mRNA_NpADO
    d_NpADO = p['k_protein_prod_NpADO'] * mRNA_NpADO - p['d_protein_NpADO'] * NpADO
    d_Alkane = p['k_Alkane_prod'] * NpADO - p['d_Alkane'] * Alkane
    
    return [d_FadR, d_FA, d_Complex_FA, d_mRNA_NpADO, d_NpADO, d_Alkane]


# --- Main simulation and plotting function ---
def run_simulation_and_plot():
    """
    Runs the simulations for both pathways and plots the results in separate subplots.
    """
    # Define initial conditions
    initial_conditions_faee = [10, 0, 0, 0, 0, 0] # LacI, Allo, Complex, mRNA, AEAT, FAEE
    initial_conditions_alkane = [10, 0, 0, 0, 0, 0] # FadR, FA, Complex, mRNA, NpADO, Alkane

    # Solve the ODEs for the FAEE pathway
    sol_faee = solve_ivp(
        lambda t, y: faee_odes(t, y, params),
        t_span,
        initial_conditions_faee,
        t_eval=t_eval
    )

    # Solve the ODEs for the Alkane pathway
    sol_alkane = solve_ivp(
        lambda t, y: alkane_odes(t, y, params),
        t_span,
        initial_conditions_alkane,
        t_eval=t_eval
    )

    # Create a figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 6))
    fig.suptitle('Dynamic Gene Circuit Simulations', fontsize=16)

    # Plot results for the FAEE pathway
    ax1.plot(sol_faee.t, sol_faee.y[0], label='LacI Repressor', color='#3b82f6')
    ax1.plot(sol_faee.t, sol_faee.y[1], label='Allolactose Inducer', color='#f97316')
    ax1.plot(sol_faee.t, sol_faee.y[3], label='AEAT mRNA', color='#10b981')
    ax1.plot(sol_faee.t, sol_faee.y[4], label='AEAT Enzyme', color='#f59e0b')
    ax1.plot(sol_faee.t, sol_faee.y[5], label='FAEE Product', color='#ef4444')
    ax1.set_title('FAEE Production from Lactose')
    ax1.set_xlabel('Time (a.u.)')
    ax1.set_ylabel('Concentration (a.u.)')
    ax1.legend()
    ax1.grid(True, linestyle='--')

    # Plot results for the Alkane pathway
    ax2.plot(sol_alkane.t, sol_alkane.y[0], label='FadR Repressor', color='#3b82f6')
    ax2.plot(sol_alkane.t, sol_alkane.y[1], label='Fatty Acid Inducer', color='#f97316')
    ax2.plot(sol_alkane.t, sol_alkane.y[3], label='NpADO mRNA', color='#10b981')
    ax2.plot(sol_alkane.t, sol_alkane.y[4], label='NpADO Enzyme', color='#f59e0b')
    ax2.plot(sol_alkane.t, sol_alkane.y[5], label='Alkane Product', color='#ef4444')
    ax2.set_title('Alkane Production from Fatty Acids')
    ax2.set_xlabel('Time (a.u.)')
    ax2.set_ylabel('Concentration (a.u.)')
    ax2.legend()
    ax2.grid(True, linestyle='--')

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

# Run the simulation if the script is executed directly
if __name__ == '__main__':
    run_simulation_and_plot()
