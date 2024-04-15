import os
import sys

# Define the parameter names
parameters = ["NMCRY", "NMNAT", "HET", "MIG1", "MIG2"]

# Initialize a dictionary to store the values for each parameter
param_values = {param: [] for param in parameters}

# Iterate through each directory and collect parameter values
for i in range(1, 17):
    dir_name = f"bs{i}/bestrun/recent_geneflow.bestlhoods"
    if not os.path.exists(dir_name):
        print(f"Error: Directory or file '{dir_name}' does not exist. Exiting.")
        sys.exit(1)

    with open(dir_name, 'r') as file:
        lines = file.readlines()
        # Assuming the values are in the second line
        values = lines[1].split()
        for j, param in enumerate(parameters):
            param_values[param].append(float(values[j]))

# Calculate 95% confidence intervals
ci = {}
for param in parameters:
    sorted_values = sorted(param_values[param])
    lower = sorted_values[0]  # Approximation for the 2.5th percentile
    upper = sorted_values[-1] # Approximation for the 97.5th percentile
    ci[param] = (lower, upper)

# Print the confidence intervals
for param, interval in ci.items():
    print(f"{param} 95% CI: {interval}")

