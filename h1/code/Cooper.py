from SingleWeber import LocPlanProb2d
import numpy as np
import matplotlib.pyplot as plt

# np.random.seed(3)

# Coopers algorithm:
def cooper(points, p, weights=None, init_X=None, print_plot=False, print_last=False):
    assert len(points) > p, "Number of facilities exceeds the number of points"
    colors = ["red", "blue", "yellow", "pink", "brown", "grey", "cyan", "purple", "orange", "black"]
    it = 0
    # Initial locations randomly in convex hull
    if init_X is None:
        X = []
        for i in range(p):
            scalars = np.exp(10 * np.random.random((len(points), 1)))
            scalars = scalars / np.sum(scalars)
            X.append((scalars.T @ points).reshape(2, ))
        X = np.array(X)
    else:
        X = init_X

    if print_plot:
        plt.plot(*points.T, "ro")
        plt.plot(*X.T, "go")
        plt.show()

    max_change = 1
    while it == 0 or max_change > 0.0001:
        # Assign to closest facility
        points_split = [[] for _ in range(p)]
        weights_split = [[] for _ in range(p)]
        for weight, point in zip(weights, points):
            dists = np.linalg.norm(point - X, axis=1)
            dmin = np.argmin(dists)
            points_split[dmin].append(point)
            weights_split[dmin].append(weight)

        # Plot assignment
        if print_plot:
            for i in range(p):
                plt.plot(*zip(*points_split[i]), colors[i])
                plt.plot(*X[i], colors[i], markersize=10)
            plt.plot(*X.T, "go")
            plt.show()

        # Solving single Weber
        new_X = [None for _ in range(p)]
        obj = 0
        max_d = 0
        for i in range(p):
            if len(points_split[i]) > 0:  # If cluster has no points move facility to a random new location
                sub_prob = LocPlanProb2d(points=points_split[i], weights=weights_split[i])
                # sub_prob.solve_ostresh(lam='drezner')
                sub_prob.solve_ostresh(lam=1)
                new_X[i] = sub_prob.sol
                obj += sub_prob.obj
                pot_max_d = np.max(np.linalg.norm(sub_prob.sol - points_split[i], axis=1))
                max_d = np.max([max_d, pot_max_d])
            else:
                scalars = np.exp(10 * np.random.random((len(points), 1)))
                scalars = scalars / np.sum(scalars)
                new_X[i] = (scalars.T @ points).reshape(2, )


        new_X = np.array(new_X)
        max_change = np.max(np.linalg.norm(new_X - X, axis=1))
        X = new_X
        it += 1



        # Plot facility locations
        if print_last:
            for i in range(p):
                plt.plot(*zip(*points_split[i]), "o", color=colors[i])
                plt.plot(*X[i], "o", color=colors[i], markersize=10)
            plt.plot(*X.T, "go")
            plt.show()

    return obj, X, max_d


# my_points = np.array([[1, 2], [2, 1], [2, 4], [3, 4], [3, 3], [5, 3], [5, 2], [4, 2]])
# my_points = np.random.randint(1, 20, (20, 2))
# my_weights = np.random.random((20, ))
# my_p = 3
#
# my_X = cooper(my_points, my_p, weights=my_weights)
# print(my_X)
