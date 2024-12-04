# %%
"""
Introduction class in metabolic modeling
Autor: Haris Zafeiropoulos
Date: Dec 2024
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import plotly.io as pio
import plotly.graph_objects as go
import plotly.io as pio


# %% 3D plot with uniform points and constraints

"""
In this plot, we visualize how the flux space of a system of 3 reactions looks like:

"""
def samples_in_constraint_3D(r1_bounds, r2_bounds, r3_bounds, r3_relationship):

    # Set the renderer to 'browser'
    pio.renderers.default = 'browser'

    # Define the range for R1 and R2
    R1_range = np.linspace(r1_bounds[0], r1_bounds[1], 100)
    R2_range = np.linspace(r2_bounds[0], r2_bounds[1], 100)

    # Create a meshgrid of R1 and R2
    R1, R2 = np.meshgrid(R1_range, R2_range)

    # Calculate R3_max (the plane) for each (R1, R2) pair
    R3_min = r3_bounds[0]
    R3_max = r3_bounds[1] + r3_relationship[0] * R1 + r3_relationship[1] * R2

    print(R3_max)

    # Number of random points to generate
    num_random_points = 5000

    # Randomly sample points for R1, R2, and R3
    R1_points = np.random.uniform(R1_range[0], R1_range[-1], num_random_points)
    R2_points = np.random.uniform(R2_range[0], R2_range[-1], num_random_points)

    # Calculate the optimal R3 (e.g., maximum R3)
    optimal_index = np.argmax(R3_max)
    optimal_R1 = R1.flatten()[optimal_index]
    optimal_R2 = R2.flatten()[optimal_index]
    optimal_R3 = R3_max.flatten()[optimal_index]

    R1_min, R1_max = R1_range[0], R1_range[-1]
    R2_min, R2_max = R2_range[0], R2_range[-1],

    # Calculate the corresponding R3 values, ensuring they are within the valid range (0 to R3_max)
    R3_points = np.random.uniform(R3_min, r3_bounds[1] + r3_relationship[0] * R1_points + r3_relationship[1] * R2_points)

    # Plot the random points in 3D space
    fig = go.Figure(data=[go.Scatter3d(
        x=R1_points,
        y=R2_points,
        z=R3_points,
        mode='markers',
        marker=dict(
            size=3,
            color=R3_points,  # Color by R3 values
            colorscale='Viridis',
            colorbar=dict(title='R3 Value', x=0),  # Add a colorbar with title
            opacity=0.8,
        )
    )])

    # Add the surface for R3_max (the constraint plane)
    fig.add_trace(go.Surface(
        x=R1,
        y=R2,
        z=R3_max,
        colorscale='Viridis',
        showscale=False,  # Disable the color bar
        opacity=0.8,
        cmin=0,
        cmax=20
    ))

    # R1 constraints (at R1_min and R1_max)
    fig.add_trace(go.Scatter3d(
        x=[R1_min, R1_max],
        y=[R2_min, R2_min],
        # z=[R1_max - 0.5*R1_min - 0.3*R2_min, 10 - 0.5*R1_max - 0.3*R2_min],
        z=[R1_max + r3_relationship[0]*R1_min + r3_relationship[1]*R2_min,
           R1_max + r3_relationship[0]*R1_max + r3_relationship[1]*R2_min],
        mode='lines',
        line=dict(color='red', width=8),
        name='R1 Constraint'
    ))

    # R2 constraints (at R2_min and R2_max)
    fig.add_trace(go.Scatter3d(
        x=[R1_min, R1_min],
        y=[R2_min, R2_max],
        z=[R1_max + r3_relationship[0]*R1_min + r3_relationship[1]*R2_min,
           R1_max + r3_relationship[0]*R1_min + r3_relationship[1]*R2_max],
        mode='lines',
        line=dict(color='green', width=8),
        name='R2 Constraint'
    ))

    # R3 constraints (at R3_min and R3_max)
    fig.add_trace(go.Scatter3d(
        x=[R1_min, R1_min],
        y=[R2_min, R2_min],
        z=[R3_min, np.max(R3_max)],
        mode='lines',
        line=dict(color='blue', width=8),
        name='R3 Constraint'
    ))

    # Add a marker for the optimal value of R3
    fig.add_trace(go.Scatter3d(
        x=[optimal_R1],
        y=[optimal_R2],
        z=[optimal_R3],
        mode='markers+text',
        marker=dict(size=8, color='orange'),
        text=['Optimal R3'],  # Text label for optimal point
        textposition='top center',
        name='Optimal R3'
    ))

    # Title and axis labels
    fig.update_layout(
        scene=dict(
            xaxis_title='R1',
            yaxis_title='R2',
            zaxis_title='R3'
        ),
        title="Solution Space Below the Plane"
    )

    return fig.show()


# %% FBA with 3 reactions and a plane

class FluxBalanceAnalysis3D:
    def __init__(self, stoich_matrix, bounds, objective):
        """
        Initialize the FBA class for 3D visualization.
        Args:
            stoich_matrix (np.ndarray): Stoichiometric matrix (m x n).
            bounds (list of tuples): Bounds for each flux (n).
            objective (np.ndarray): Coefficients of the objective function (n).
        """
        pio.renderers.default = 'notebook'
        self.S = stoich_matrix
        self.bounds = bounds
        self.objective = objective

    def optimize(self):
        """Simple optimization using numpy (toy solver)."""
        from scipy.optimize import linprog

        # Linear programming: maximize c^T v => minimize -c^T v
        result = linprog(
            c=-self.objective,
            A_eq=self.S,
            b_eq=np.zeros(self.S.shape[0]),
            bounds=self.bounds,
            method='highs'
        )
        self.solution = result
        return result

    def visualize_geometry_3d(self):
        """Visualize the flux space and feasible region in 3D."""
        if self.S.shape[1] != 3:
            raise ValueError("Visualization is limited to 3D problems.")

        # Generate a grid for v1 and v2 within bounds
        v1 = np.linspace(self.bounds[0][0], self.bounds[0][1], 50)
        v2 = np.linspace(self.bounds[1][0], self.bounds[1][1], 50)
        v1, v2 = np.meshgrid(v1, v2)

        # Calculate v3 based on S * v = 0
        v3 = -(self.S[0, 0] * v1 + self.S[0, 1] * v2) / self.S[0, 2]

        # Mask infeasible points (where v3 is outside its bounds)
        mask = (v3 >= self.bounds[2][0]) & (v3 <= self.bounds[2][1])
        v3 = np.where(mask, v3, np.nan)

        # Plot 3D feasible region
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(v1, v2, v3, color='lightblue', alpha=0.7, edgecolor='k', linewidth=0.2)

        # Add constraint lines
        for i in range(self.S.shape[0]):
            # Constraint plane: S[i,0]*v1 + S[i,1]*v2 + S[i,2]*v3 = 0
            v3_constraint = -(self.S[i, 0] * v1 + self.S[i, 1] * v2) / self.S[i, 2]
            ax.plot_wireframe(v1, v2, v3_constraint, color=f"C{i}", label=f"Constraint {i + 1}", linewidth=1)

        # Set axis labels
        ax.set_xlabel('Flux v1', fontsize=14, labelpad=20)
        ax.set_ylabel('Flux v2', fontsize=14, labelpad=20)
        ax.set_zlabel('Flux v3', fontsize=14, labelpad=20)

        # Adjust layout to avoid clipping
        fig.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.2)

        # Highlight the optimal solution if available
        if hasattr(self, 'solution') and self.solution.success:
            ax.scatter(
                self.solution.x[0], self.solution.x[1], self.solution.x[2],
                color="red", s=50, label="Optimal solution"
            )
            ax.text(
                self.solution.x[0], self.solution.x[1], self.solution.x[2],
                "Optimal", color="red"
            )

        plt.title("3D Flux Space with Constraints", pad=30)
        plt.legend()
        plt.show()


# %% 3 reactions with constraint plain

def constrain_plane_in_3D():
    pio.renderers.default = 'browser'

    # Define the constraints for the three reactions
    R1_min, R1_max = -10, 10
    R2_min, R2_max = -5, 15
    R3_min, R3_max = 0, 20

    # Generate grid points for the fluxes
    R1 = np.linspace(R1_min, R1_max, 50)
    R2 = np.linspace(R2_min, R2_max, 50)
    R1, R2 = np.meshgrid(R1, R2)

    # Example constraint: R3 = 10 - 0.5*R1 - 0.3*R2
    R3 = 10 - 0.5 * R1 - 0.3 * R2

    # Mask to ensure fluxes lie within bounds
    valid = (R3 >= R3_min) & (R3 <= R3_max)
    R1 = R1[valid]
    R2 = R2[valid]
    R3 = R3[valid]


    # Calculate the optimal R3 (e.g., maximum R3)
    R3_optimal_index = np.argmax(R3)
    R1_when_optimal_R3 = R1.flatten()[R3_optimal_index]
    R2_when_optimal_R3 = R2.flatten()[R3_optimal_index]
    optimal_R3 = R3.flatten()[R3_optimal_index]


    # Create a 3D scatter plot for valid points
    fig = go.Figure()

    # Add the feasible region as a scatter plot
    fig.add_trace(go.Scatter3d(
        x=R1.flatten(),
        y=R2.flatten(),
        z=R3.flatten(),
        mode='markers',
        marker=dict(
            size=4,
            color='lightblue',
            opacity=0.6
        ),
        name='Feasible Region'
    ))

    # Add lines for reaction constraints
    # R1 constraints (at R1_min and R1_max)
    fig.add_trace(go.Scatter3d(
        x=[R1_min, R1_max],
        y=[R2_min, R2_min],
        z=[10 - 0.5*R1_min - 0.3*R2_min, 10 - 0.5*R1_max - 0.3*R2_min],
        mode='lines',
        line=dict(color='red', width=8),
        name='R1 Constraint'
    ))

    # R2 constraints (at R2_min and R2_max)
    fig.add_trace(go.Scatter3d(
        x=[R1_min, R1_min],
        y=[R2_min, R2_max],
        z=[10 - 0.5*R1_min - 0.3*R2_min, 10 - 0.5*R1_min - 0.3*R2_max],
        mode='lines',
        line=dict(color='green', width=8),
        name='R2 Constraint'
    ))

    # R3 constraints (at R3_min and R3_max)
    fig.add_trace(go.Scatter3d(
        x=[R1_min, R1_min],
        y=[R2_min, R2_min],
        z=[R3_min, R3_max],
        mode='lines',
        line=dict(color='blue', width=8),
        name='R3 Constraint'
    ))

    # Add a marker for the optimal value of R3
    fig.add_trace(go.Scatter3d(
        x=[R1_when_optimal_R3],
        y=[R2_when_optimal_R3],
        z=[optimal_R3],
        mode='markers+text',
        marker=dict(size=8, color='orange'),
        text=['Optimal R3'],  # Text label for optimal point
        textposition='top center',
        name='Optimal R3'
    ))


    # Add axis labels and title
    fig.update_layout(
        title='Constraint-Based Model Visualization with Constraints',
        scene=dict(
            xaxis_title='Flux R1',
            yaxis_title='Flux R2',
            zaxis_title='Flux R3'
        ),
        width=1100,
        height=900
    )

    fig.show()


# %% Basic component

class FluxBalanceAnalysis:
    def __init__(self, stoich_matrix, bounds, objective):
        """
        Initialize the FBA class.
        Args:
            stoich_matrix (np.ndarray): Stoichiometric matrix (m x n).
            bounds (list of tuples): Bounds for each flux (n).
            objective (np.ndarray): Coefficients of the objective function (n).
        """
        self.S = stoich_matrix
        self.bounds = bounds
        self.objective = objective

    def optimize(self):
        """Simple optimization using numpy (toy solver)."""
        from scipy.optimize import linprog

        # Linear programming: maximize c^T v => minimize -c^T v
        result = linprog(
            c=-self.objective,
            A_eq=self.S,
            b_eq=np.zeros(self.S.shape[0]),
            bounds=self.bounds,
            method='highs'
        )
        self.solution = result
        return result
    def visualize_geometry(self):
        """Visualize the flux space (line) and feasible region (segment) in 2D."""
        if self.S.shape[1] != 2:
            raise ValueError("Visualization is limited to 2D problems.")

        v1 = np.linspace(self.bounds[0][0], self.bounds[0][1], 400)  # Flux v1 range based on its bounds
        v2 = -(self.S[0, 0] * v1) / self.S[0, 1]  # Flux space from S * v = 0

        # Filter feasible v2 based on bounds
        feasible = (v2 >= self.bounds[1][0]) & (v2 <= self.bounds[1][1])
        v1_feasible = v1[feasible]
        v2_feasible = v2[feasible]

        plt.figure(figsize=(8, 8))

        # Plot the full flux space (line)
        plt.plot(v1, v2, color="blue", linestyle="--", label="Flux space (S · v = 0)")

        # Highlight the feasible region (segment)
        plt.plot(v1_feasible, v2_feasible, color="blue", label="Feasible region")
        plt.scatter(v1_feasible[0], v2_feasible[0], color="green", label="Feasible region start")
        plt.scatter(v1_feasible[-1], v2_feasible[-1], color="purple", label="Feasible region end")

        # Add axes
        plt.axhline(0, color='black', linewidth=0.8)
        plt.axvline(0, color='black', linewidth=0.8)

        # Plot the optimal solution if available
        if hasattr(self, 'solution') and self.solution.success:
            plt.scatter(
                self.solution.x[0], self.solution.x[1], color="red", label="Optimal solution"
            )
            plt.annotate(
                "Optimal\n(%.2f, %.2f)" % tuple(self.solution.x),
                (self.solution.x[0], self.solution.x[1]),
                textcoords="offset points", xytext=(10, -15), ha="center"
            )

        plt.xlim(self.bounds[0][0], self.bounds[0][1])
        plt.ylim(self.bounds[1][0], self.bounds[1][1])
        plt.xlabel("Flux v1")
        plt.ylabel("Flux v2")
        plt.title("Flux Space and Feasible Region")
        plt.legend()
        plt.grid()
        plt.show()



if __name__ == "__main__":

    # Example Use of FluxBalanceAnalysis

    # Stoichiometric matrix for a single reaction: A -> B, constraints are -v1 + v2 = 0
    S = np.array([[1, -1]])

    # Bounds: flux 1 between 0 and 5, flux 2 between 0 and 5
    bounds = [(-2, 15), (0, 5)]

    # Objective: Maximize v2 (e.g., production of B)
    objective = np.array([0, 1])

    # Create the FBA instance
    fba = FluxBalanceAnalysis(S, bounds, objective)

    # Run optimization
    result = fba.optimize()
    print("Optimization Result:", result)

    # Visualize the geometry
    fba.visualize_geometry()





    # Example Usage FluxBalanceAnalysis3D

    # Stoichiometric matrix for a single reaction: A + B -> C
    S = np.array([[1, 1, -1]])  # A + B -> C

    # Bounds: fluxes v1, v2, and v3 between 0 and 10
    bounds = [(-2, 5), (-10, 0), (0, 12)]

    # Objective: Maximize v3 (e.g., production of C)
    objective = np.array([0, 0, 1])

    # Create the FBA instance
    fba = FluxBalanceAnalysis3D(S, bounds, objective)

    # Run optimization
    result = fba.optimize()
    print("Optimization Result:", result)

    # Visualize the geometry in 3D
    fba.visualize_geometry_3d()



