import numpy as np
from scipy.optimize import minimize
from itertools import combinations
import matplotlib.pyplot as plt

def objective(vars, d):
    # Split variables into points
    points = vars.reshape((d + 2, -1))  # d + 2 points in d-dimensional space

    common_face_points, apex1, apex2 = points[:d], points[d], points[d + 1]

    if d == 1:
        vec1, vec2 = apex1 - common_face_points[0], apex2 - common_face_points[0]
        vec1 /= np.linalg.norm(vec1)
        vec2 /= np.linalg.norm(vec2)
        angle = np.arccos(np.clip(np.dot(vec1, vec2), -1.0, 1.0))
    else:
        shared_face_vectors = [points[i] - points[0] for i in range(1, d)]
        P = np.linalg.qr(np.column_stack(shared_face_vectors))[0]
        apex_vector1 = apex1 - points[0]
        apex_vector2 = apex2 - points[0]
        proj_apex_vector1 = P @ (P.T @ apex_vector1)
        proj_apex_vector2 = P @ (P.T @ apex_vector2)
        res_apex_vector1 = apex_vector1 - proj_apex_vector1
        res_apex_vector2 = apex_vector2 - proj_apex_vector2
        res_apex_vector1 /= np.linalg.norm(res_apex_vector1)
        res_apex_vector2 /= np.linalg.norm(res_apex_vector2)
        angle = np.arccos(np.clip(np.dot(res_apex_vector1, res_apex_vector2), -1.0, 1.0))
    return angle


def constraint_edge_lengths(vars, d, q):
    """
    Inequality constraints ensuring the ratio of pairwise edge lengths is within [q, 1/q],
    excluding the apex-apex pair (d, d+1).
    """
    points = vars.reshape((d + 2, -1))
    constraints = []

    # Generate all pairs from {0, ..., d+1}, excluding (d, d+1)
    all_pairs = list(combinations(range(d + 2), 2))
    valid_pairs = [pair for pair in all_pairs if pair != (d, d + 1)]

    # Compute distances for valid pairs
    distances = {pair: np.linalg.norm(points[pair[0]] - points[pair[1]]) for pair in valid_pairs}

    # Add inequality constraints for each valid pair combination
    for pair_1, pair_2 in combinations(valid_pairs, 2):
        dist_ij, dist_kl = distances[pair_1], distances[pair_2]
        constraints.append(dist_ij - dist_kl * q)  # Ensure dist_ij >= dist_kl * q
        constraints.append(dist_kl - dist_ij * q)  # Ensure dist_kl >= dist_ij * q

    return constraints


#def constraint_bounds(vars, d, tau):
#    points = vars.reshape((d + 2, -1))
#    return [tau - abs(point[-1]) for point in points]

def constraint_bounds(vars, d, tau, D):
    """
    Updated constraint function to handle a general ambient dimension D >= d+1.
    - For dimensions from d+1 to D, constrain their values by tau / (D-d).
    """
    points = vars.reshape((d+2, -1))  # Reshape vars into points in D-dimensional space
    constraints = []

    # Apply constraints for dimensions d+1 to D
    #scaling_factor = tau / (D - d)

    #for point in points:
    #    for dim in range(d, D):
    #        constraints.append(scaling_factor - abs(point[dim]))
    for point in points:
        # Add the constraint that the sum of squares of dimensions d+1 to D is less than tau
        sum_of_squares = np.sum(point[d:D] ** 2)
        constraints.append(tau**2 - sum_of_squares)

    return constraints


def constraint_angle(vars, d):
    points = vars.reshape((d + 2, -1))
    common_face_points, apex1, apex2 = points[:d], points[d], points[d + 1]

    if d == 1:
        vec1, vec2 = apex1 - common_face_points[0], apex2 - common_face_points[0]
        vec1 /= np.linalg.norm(vec1)
        vec2 /= np.linalg.norm(vec2)
        angle = np.arccos(np.clip(np.dot(vec1, vec2), -1.0, 1.0))
    else:
        shared_face_vectors = [points[i] - points[0] for i in range(1, d)]
        P = np.linalg.qr(np.column_stack(shared_face_vectors))[0]
        apex_vector1 = apex1 - points[0]
        apex_vector2 = apex2 - points[0]
        proj_apex_vector1 = P @ (P.T @ apex_vector1)
        proj_apex_vector2 = P @ (P.T @ apex_vector2)
        res_apex_vector1 = apex_vector1 - proj_apex_vector1
        res_apex_vector2 = apex_vector2 - proj_apex_vector2
        res_apex_vector1 /= np.linalg.norm(res_apex_vector1)
        res_apex_vector2 /= np.linalg.norm(res_apex_vector2)
        angle = np.arccos(np.clip(np.dot(res_apex_vector1, res_apex_vector2), -1.0, 1.0))
    return angle - np.pi / 2

def generate_feasible_initial_guess(num_points, d, D, tau, scale=0.1):
    """
    Generate an initial guess for optimization that respects the sum-of-squares constraint
    for dimensions d+1 to D.

    Parameters:
    - num_points: Total number of points (d+2 for the simplex problem).
    - d: Dimension of the simplex.
    - D: Ambient dimension.
    - tau: Constraint parameter.
    - scale: Scaling factor for random initialization in dimensions 0 to d.

    Returns:
    - Initial guess as a flattened array.
    """
    points = np.zeros((num_points, D))

    # Random initialization for first d dimensions
    points[:, :d] = np.random.uniform(-scale, scale, size=(num_points, d))

    # For dimensions d+1 to D, generate vectors with norm <= tau
    for i in range(num_points):
        random_vector = np.random.normal(size=(D - d))
        random_vector /= np.linalg.norm(random_vector)  # Normalize
        random_vector *= np.sqrt(np.random.uniform(0, tau**2))  # Scale within [0, tau^2]
        points[i, d:D] = random_vector

    return points.flatten()



def optimize_d_simplex(d, D, e=1.0, tau=0.1, q=1, max_retries=100):
    """
    Retry optimization until it succeeds or reaches max_retries.
    """
    dim = D
    num_points = d + 2
    retries = 0

    while retries < max_retries:
        # Generate a new random initial guess

        print(f"Current attempt: {retries + 1}")

        #initial_guess = np.random.rand(num_points * dim).flatten() * 0.1
        initial_guess = generate_feasible_initial_guess(num_points, d, D, tau, scale=0.1)

        # Define constraints
        constraints = [
            {'type': 'eq', 'fun': lambda vars: np.linalg.norm(vars[dim:2*dim] - vars[:dim]) - e},
            {'type': 'ineq', 'fun': lambda vars: constraint_edge_lengths(vars, d, q)},
            {'type': 'ineq', 'fun': lambda vars: constraint_bounds(vars, d, tau, D)},
            {'type': 'ineq', 'fun': lambda vars: constraint_angle(vars, d)},
        ]

        # Attempt optimization
        result = minimize(
            lambda vars: objective(vars, d),
            initial_guess,
            constraints=constraints,
            method='trust-constr',
            tol=1e-6
        )

        if result.success:
            #print("Optimization succeeded!")

            #print(f"Initial guess:\n{initial_guess.reshape((d + 2, -1))}")

            return {
                "min_angle": result.fun,
                "wLAPD": np.pi - result.fun,
                "solution": result.x.reshape((d + 2, -1))
            }

        retries += 1
        #print("Optimization failed, retrying...")

    # If max retries reached, raise an error
    raise ValueError("Optimization failed after maximum retries!")

def optimize_multiple_runs(num_runs, d, D, e=1.0, tau=0.1, q=1):
    best_result = None
    largest_wLAPD = -np.inf  # Start with a very low value

    for i in range(num_runs):
        print(f"Run {i + 1}/{num_runs}")

        try:
            # Run the optimization
            current_result = optimize_d_simplex(d=d, D=D, e=e, tau=tau, q=q)

            # Check if the new complement angle is better
            if current_result["wLAPD"] > largest_wLAPD:
                result = current_result
                #largest_wLAPD = current_result["wLAPD"]
                best_result = current_result["wLAPD"]
                best_solution = current_result["solution"]

        except ValueError as e:
            print(f"Run {i + 1} failed: {e}")

    # After all runs, return the best result
    if best_result is not None:
        #print("Best result:")
        print(f"Largest wLAPD found: {best_result}")
        #print("Solution:")
        #print(best_solution)
    else:
        print("No successful optimizations after all runs.")

    return result


# Example usage
d, noise = 2, 0.1
final_result = optimize_multiple_runs(num_runs=1, d=d, D=d+1, e=1, tau=noise, q=1/1)
print(f"Largest wLAPD: {final_result['wLAPD']}")
print(f"Solution: {final_result['solution']}")


#Example usage
#result = optimize_d_simplex(d=3, D=4, e=1, tau=0.1, q=1/1)
#print(f"Minimum angle under constraints (in radians): {result['min_angle']}")
#print(f"Largest wLAPD: {result['wLAPD']}")
#print("Solution:")
#print(result['solution'])



################## plot ##############

# results = []
#
# # Loop through D from 3 to 30
# for D in range(3, 6):
#     print(f"Now running D={D}")
#     try:
#         result = optimize_multiple_runs(num_runs=100, d=2, D=D, e=1, tau=0.1, q=1)
#         complement_angle = result['wLAPD']
#         results.append((D, complement_angle))
#     except ValueError as e:
#         print(f"Optimization failed for D={D}: {e}")
#         results.append((D, None))
#
# # Separate valid results for plotting
# D_values = [item[0] for item in results if item[1] is not None]
# wLAPD = [item[1] for item in results if item[1] is not None]
#
# # Plot the results
# plt.figure(figsize=(10, 6))
# plt.plot(D_values, wLAPD, marker='o', label='WLAPD')
# plt.xlabel('Ambient Dimension (D)')
# plt.ylabel('WLAPD (radians)')
# plt.title('WLAPD vs Ambient Dimension (D)')
# plt.grid(True)
# plt.legend()
# plt.show()
