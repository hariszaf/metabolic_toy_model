import os
import cobra
from pathlib import Path

# # Set the gurobi license
# os.environ['GRB_LICENSE_FILE'] = 'licenses/gurobi.lic'


root_dir = Path(__file__).resolve().parent.parent

modelfile = os.path.join(root_dir, "files/models/e_coli_core.xml")

model = cobra.io.read_sbml_model(modelfile)


solution = model.optimize()

print(f"flux balance analysis solution is {solution.objective_value}")

print("COBRApy is working", "\U0001F600")


# Set parameter WLSAccessID
# Set parameter WLSSecret
# Set parameter LicenseID to value 940603
# Academic license 940603 - for non-commercial use only - registered to da___@gmail.com
# flux balance analysis solution is 0.8739215069684301
# COBRApy is working 😀
# Warning: environment still referenced so free is deferred (Continue to use WLS)
