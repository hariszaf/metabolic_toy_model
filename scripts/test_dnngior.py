import os
import cobra
from pathlib import Path


# Get root directory
root_dir = Path(__file__).resolve().parent.parent
# Set the gurobi license
os.environ['GRB_LICENSE_FILE'] = os.path.join(root_dir, 'licenses/gurobi.lic')

import dnngior

path_to_blautia_model = os.path.join(root_dir, "files/models/bh_ungapfilled_model.sbml")

# Load draft reconstruction to be gap-filled
draft_reconstruction = cobra.io.read_sbml_model(path_to_blautia_model)

print(f"The solution of the draft reconstruction is:\n\n{draft_reconstruction.optimize()}",  "\U0001F63F", "\n")

# Gap fill with dnngior
gapfill = dnngior.Gapfill(
    draftModel = path_to_blautia_model,
    medium = None,
    objectiveName = 'bio1'
)

print("Number of reactions added:", len(gapfill.added_reactions), "\n~~")

# Check reactions added
for reaction in gapfill.added_reactions:
    print(f"Name: {gapfill.gapfilledModel.reactions.get_by_id(reaction).name}", "\n", f"Equation: {gapfill.gapfilledModel.reactions.get_by_id(reaction).build_reaction_string(use_metabolite_names=True)}")

# Get optimal value of the objective function for the gap-filled model
gf_model = gapfill.gapfilledModel.copy()

print(f"The solution of the gapfilled reconstruction is:\n\n{gf_model.optimize()}",  "\U0001F63A", "\n")

# The solution of the gapfilled reconstruction is:
# <Solution 154.228 at 0x189db678b10> 😺
