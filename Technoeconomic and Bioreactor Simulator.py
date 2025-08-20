# ==============================================================================
# Biofuel Production Models with Data Visualization (Python Script)
#
# This script contains two computational models for a biofuel production project.
# It now includes data visualization using matplotlib to help interpret the
# model outputs and a new, interactive sensitivity analysis for key parameters.
#
# 1. Bioreactor Optimization Model: Calculates total biofuel yield and plots
#    production over time.
# 2. Technoeconomic Analysis (TEA) Model: Estimates financial viability and plots
#    Minimum Selling Price (MSP) against annual production.
# 3. Interactive Sensitivity Analysis: Allows the user to set parameters and
#    ranges to plot MSP vs. varying substrate cost and product yield.
#
# To run this script, you must have the 'matplotlib' library installed.
# You can install it by running: pip install matplotlib
# ==============================================================================

import matplotlib.pyplot as plt
import numpy as np

def get_valid_input(prompt, value_type=float, min_val=None, max_val=None):
    """
    Helper function to get and validate user input.
    Ensures the input is a valid number within an optional range.
    """
    while True:
        try:
            user_input = value_type(input(prompt))
            if min_val is not None and user_input < min_val:
                print(f"Error: Value must be at least {min_val}.")
            elif max_val is not None and user_input > max_val:
                print(f"Error: Value must not exceed {max_val}.")
            else:
                return user_input
        except ValueError:
            print("Invalid input. Please enter a valid number.")

def bioreactor_model():
    """
    Runs the Bioreactor Optimization Model with visualization.
    Prompts for bioreactor parameters, calculates total biofuel production,
    and plots production over batch time.
    """
    print("\n--- Bioreactor Optimization Model ---")
    print("Adjust the parameters to see the impact on total biofuel produced.")
    
    # Get user inputs for bioreactor parameters
    volume = get_valid_input("Bioreactor Volume (L, e.g., 100000): ", value_type=int, min_val=1)
    substrate_conc = get_valid_input("Substrate Concentration (g/L, e.g., 100): ", min_val=0)
    yield_factor = get_valid_input("Product Yield (g product / g substrate, e.g., 0.2): ", min_val=0.01)
    growth_rate = get_valid_input("Growth Rate (1/h, e.g., 0.5): ", min_val=0.1)
    batch_time = get_valid_input("Batch Time (h, e.g., 48): ", value_type=int, min_val=1)
    
    # Simplified calculation for substrate consumed (using exponential decay model)
    consumed_substrate = substrate_conc * (1 - np.exp(-growth_rate * batch_time))
    
    # Calculate the total biofuel produced
    total_biofuel = volume * consumed_substrate * yield_factor
    
    print("\n--- Results ---")
    print(f"Total Biofuel Produced: {total_biofuel:,.0f} g")
    print("(Note: This is a simplified model for demonstration purposes)")
    
    # --- Data Visualization for Bioreactor Model ---
    print("Generating plot for Bioreactor Production vs. Batch Time...")
    
    # Generate data for the plot
    time_points = np.linspace(0, batch_time * 1.5, 100) # Plot up to 150% of the entered time
    production_over_time = volume * (substrate_conc * (1 - np.exp(-growth_rate * time_points))) * yield_factor
    
    # Create the plot
    plt.figure(figsize=(10, 6))
    plt.plot(time_points, production_over_time, label='Biofuel Production')
    plt.plot(batch_time, total_biofuel, 'ro', label=f'Current Run: {batch_time}h') # Mark the user's specific point
    plt.title('Bioreactor Production Over Time', fontsize=16)
    plt.xlabel('Batch Time (h)', fontsize=12)
    plt.ylabel('Total Biofuel Produced (g)', fontsize=12)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

def tea_model():
    """
    Runs the Technoeconomic Analysis (TEA) Model with visualization.
    Prompts for CAPEX and OPEX values, calculates the Minimum Selling Price (MSP),
    and plots MSP against annual production.
    """
    print("\n--- Technoeconomic Analysis Model ---")
    print("Estimate the financial viability by calculating the Minimum Selling Price (MSP).")
    
    # --- Capital Expenditures (CAPEX) ---
    print("\n1. Capital Expenditures (CAPEX) - One-time costs")
    bioreactor_cost = get_valid_input("Bioreactor Cost ($, e.g., 2500000): ", value_type=int)
    dsp_cost = get_valid_input("Downstream Processing Cost ($, e.g., 1500000): ", value_type=int)
    plant_lifetime = get_valid_input("Plant Lifetime (years, e.g., 20): ", value_type=int, min_val=1)
    
    # Calculate annualized CAPEX
    annualized_capex = (bioreactor_cost + dsp_cost) / plant_lifetime
    
    # --- Operating Expenditures (OPEX) ---
    print("\n2. Operating Expenditures (OPEX) - Ongoing costs")
    annual_production = get_valid_input("Annual Biofuel Production (g, e.g., 1000000): ", value_type=int, min_val=1)
    annual_substrate_cost = get_valid_input("Annual Substrate Cost ($, e.g., 100000): ", value_type=int)
    annual_utility_cost = get_valid_input("Annual Utility & Chemical Cost ($, e.g., 50000): ", value_type=int)
    annual_labor_cost = get_valid_input("Annual Labor Cost ($, e.g., 200000): ", value_type=int)
    
    # Calculate total annual cost
    total_annual_cost = annualized_capex + annual_substrate_cost + annual_utility_cost + annual_labor_cost
    
    # Calculate Minimum Selling Price (MSP)
    msp = total_annual_cost / annual_production
    
    print("\n--- Results ---")
    print(f"Annualized CAPEX: ${annualized_capex:,.0f}")
    print(f"Total Annual Cost: ${total_annual_cost:,.0f}")
    print(f"Minimum Selling Price (MSP): ${msp:,.4f} / g")

    # --- Data Visualization for TEA Model ---
    print("Generating plot for MSP vs. Annual Production...")

    # Generate data for the plot
    production_points = np.linspace(annual_production * 0.5, annual_production * 2, 100)
    msp_over_production = total_annual_cost / production_points
    
    # Create the plot
    plt.figure(figsize=(10, 6))
    plt.plot(production_points, msp_over_production, label='MSP')
    plt.plot(annual_production, msp, 'ro', label=f'Current Run: {annual_production:,.0f} g/yr')
    plt.title('Minimum Selling Price vs. Annual Production', fontsize=16)
    plt.xlabel('Annual Production (g)', fontsize=12)
    plt.ylabel('Minimum Selling Price ($ / g)', fontsize=12)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

def sensitivity_analysis():
    """
    Performs and visualizes a sensitivity analysis of MSP vs. key parameters,
    with user-provided inputs.
    """
    print("\n--- Interactive Sensitivity Analysis ---")
    print("Please enter the base parameters for your model.")
    
    # Get base parameters from the user for the analysis
    bioreactor_cost = get_valid_input("Bioreactor Cost ($, e.g., 2500000): ", value_type=int)
    dsp_cost = get_valid_input("Downstream Processing Cost ($, e.g., 1500000): ", value_type=int)
    plant_lifetime = get_valid_input("Plant Lifetime (years, e.g., 20): ", value_type=int, min_val=1)
    base_production = get_valid_input("Base Annual Biofuel Production (g, e.g., 1000000): ", value_type=int, min_val=1)
    base_substrate_cost = get_valid_input("Base Annual Substrate Cost ($, e.g., 100000): ", value_type=int)
    base_yield = get_valid_input("Base Product Yield (g/g, e.g., 0.2): ", min_val=0.01)
    annual_utility_cost = get_valid_input("Annual Utility & Chemical Cost ($, e.g., 50000): ", value_type=int)
    annual_labor_cost = get_valid_input("Annual Labor Cost ($, e.g., 200000): ", value_type=int)
    
    # Get the range for the analysis
    print("\nNow, enter the range you want to analyze for Substrate Cost.")
    sub_cost_min = get_valid_input("Minimum Annual Substrate Cost ($, e.g., 50000): ", value_type=int, min_val=0)
    sub_cost_max = get_valid_input("Maximum Annual Substrate Cost ($, e.g., 200000): ", value_type=int, min_val=sub_cost_min)
    
    print("\nFinally, enter the range you want to analyze for Product Yield.")
    yield_min = get_valid_input("Minimum Product Yield (g/g, e.g., 0.1): ", min_val=0.01)
    yield_max = get_valid_input("Maximum Product Yield (g/g, e.g., 0.4): ", min_val=yield_min)
    
    annualized_capex = (bioreactor_cost + dsp_cost) / plant_lifetime

    # --- Varying Substrate Cost ---
    # Create a range of annual substrate costs based on user input
    substrate_costs = np.linspace(sub_cost_min, sub_cost_max, 50)
    msp_vs_substrate = (annualized_capex + substrate_costs + annual_utility_cost + annual_labor_cost) / base_production
    
    # --- Varying Product Yield ---
    # Vary the yield factor and calculate the corresponding annual production
    yield_factors = np.linspace(yield_min, yield_max, 50) # g product / g substrate
    production_from_yield = base_production * (yield_factors / base_yield)
    
    # Calculate the MSP for each production value, assuming a base substrate cost
    total_annual_cost_for_yield = annualized_capex + base_substrate_cost + annual_utility_cost + annual_labor_cost
    msp_vs_yield = total_annual_cost_for_yield / production_from_yield

    # --- Data Visualization for Sensitivity Analysis ---
    print("Generating sensitivity analysis plot...")

    plt.figure(figsize=(12, 8))
    
    # Plot MSP vs. Substrate Cost on the main axes and store the line object
    line1, = plt.plot(substrate_costs, msp_vs_substrate, 'b-', label='MSP vs. Substrate Cost')
    plt.title('Sensitivity Analysis: Minimum Selling Price (MSP)', fontsize=16)
    plt.xlabel('Annual Cost ($)', fontsize=12)
    plt.ylabel('MSP ($ / g)', fontsize=12)
    plt.grid(True)
    
    # Create a secondary x-axis and plot MSP vs. Product Yield, storing its line object
    ax2 = plt.twiny()
    ax2.set_xlabel('Product Yield (g/g)', color='r', fontsize=12)
    line2, = ax2.plot(yield_factors, msp_vs_yield, 'r--', label='MSP vs. Product Yield')
    
    # Create a single legend that includes both lines
    plt.legend(handles=[line1, line2], loc='upper right')

    plt.tight_layout()
    plt.show()

def main():
    """
    Main function to run the application loop.
    """
    while True:
        print("\n--- Biofuel Project Models ---")
        print("1. Run Bioreactor Optimization Model")
        print("2. Run Technoeconomic Analysis Model")
        print("3. Run Interactive Sensitivity Analysis")
        print("4. Exit")
        
        choice = input("Enter your choice (1, 2, 3, or 4): ")
        
        if choice == '1':
            bioreactor_model()
        elif choice == '2':
            tea_model()
        elif choice == '3':
            sensitivity_analysis()
        elif choice == '4':
            print("Exiting. Goodbye!")
            break
        else:
            print("Invalid choice. Please enter 1, 2, 3, or 4.")

# Entry point for the script
if __name__ == "__main__":
    main()
