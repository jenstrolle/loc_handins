import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from math import sqrt
from sklearn.cluster import KMeans
import sympy

ZERO = 0.0001


def center_of_3(points, weights):
    Y1, Y2, Y3 = points
    w1, w2, w3 = weights
    if w1 == w2 == w3:  # Find intersection of 2 perpendicular bisectors
        v12 = Y2 - Y1
        v13 = Y3 - Y1
        v21 = Y1 - Y2
        v23 = Y3 - Y2
        if v12 @ v13 <= 0:  # has acute angle
            Y1 = Y3  # No need to update w's as w1=w2=w3
            X = (Y1 + Y2) / 2
            n = 2
        elif v21 @ v23 <= 0:
            Y2 = Y3
            X = (Y1 + Y2) / 2
            n = 2
        # elif v31 @ v32 <= 0:  # Do nothing in this case (wont happen as this implies Y3 already covered)
        else:  # no acute angles
            c1, c2, r = sympy.symbols("c1 c2 r", real=True)
            e1 = sympy.Eq(Y2[0]**2-2*Y2[0]*c1+Y2[1]**2-2*Y2[1]*c2-Y3[0]**2+2*Y3[0]*c1-Y3[1]**2+2*Y3[1]*c2, 0)
            e2 = sympy.Eq(Y2[0]**2-2*Y2[0]*c1+Y2[1]**2-2*Y2[1]*c2-Y1[0]**2+2*Y1[0]*c1-Y1[1]**2+2*Y1[1]*c2, 0)
            res = sympy.solve([e1, e2])
            x1, x2 = float(res[c1]), float(res[c2])

            X = np.array([x1, x2])
            n = 3

    elif len({w1, w2, w3}) == 2:  # Two weights are equal -> Find intersection of line and circle
        if w1 == w3:  # making sure w1 == w2 != w3 to make it compatible with the code
            Y2, Y3 = Y3, Y2
            w2, w3 = w3, w2
        elif w2 == w3:
            Y1, Y3 = Y3, Y1
            w1, w3 = w3, w1

        # Circle of equal weighted distance to Y2 and Y3
        if w3 > w2:
            r23 = w2 / w3
            C23 = (Y3 - r23 ** 2 * Y2) / (1 - r23 ** 2)
        else:
            r23 = w3 / w2
            C23 = (Y2 - r23 ** 2 * Y3) / (1 - r23 ** 2)
        rho23 = (np.linalg.norm(Y2 - Y3) * r23) / (1 - r23 ** 2)

        a = (Y2[1] - Y1[1]) ** 2 + (Y1[0] - Y2[0]) ** 2
        b = ((Y1[0] + Y2[0] - 2 * C23[0]) * (Y2[1] - Y1[1]) + (Y1[1] + Y2[1] - 2 * C23[1]) * (Y1[0] - Y2[0]))
        c = ((Y1[0] + Y2[0]) / 2 - C23[0]) ** 2 + ((Y1[1] + Y2[1]) / 2 - C23[1]) ** 2 - rho23 ** 2
        d = b**2 - 4 * a * c  # OBS sign!

        t1 = (-b + sqrt(d)) / (2 * a)
        t2 = (-b - sqrt(d)) / (2 * a)

        x1_t1 = (Y1[0] + Y2[0]) / 2 + t1 * (Y2[1] - Y1[1])
        x2_t1 = (Y1[1] + Y2[1]) / 2 + t1 * (Y1[0] - Y2[0])

        x1_t2 = (Y1[0] + Y2[0]) / 2 + t2 * (Y2[1] - Y1[1])
        x2_t2 = (Y1[1] + Y2[1]) / 2 + t2 * (Y1[0] - Y2[0])

        X1 = np.array([x1_t1, x2_t1])
        X2 = np.array([x1_t2, x2_t2])

        # Choose the solution in the convex hull of the 3 points
        hull_original = ConvexHull(np.array([Y1, Y2, Y3]))
        hull_X1 = ConvexHull(np.array([Y1, Y2, Y3, X1]))
        # hull_X2 = ConvexHull(np.array([Y1, Y2, Y3, X2]))

        if np.all(sorted(hull_original.vertices) == sorted(hull_X1.vertices)):
            X = X1
        else:
            X = X2

    else:  # Find intersection of two circles
        for_sort = [(w1, Y1), (w2, Y2), (w3, Y3)]
        (w1, Y1), (w2, Y2), (w3, Y3) = sorted(for_sort)
        # Circle of equal weighted distance to Y2 and Y3
        r23 = w2 / w3
        C23 = (Y3 - r23 ** 2 * Y2) / (1 - r23 ** 2)
        rho23 = (np.linalg.norm(Y2 - Y3) * r23) / (1 - r23 ** 2)

        # Circle of equal weighted distance to Y1 and Y2
        r12 = w1 / w2
        C12 = (Y2 - r12 ** 2 * Y1) / (1 - r12 ** 2)
        rho12 = (np.linalg.norm(Y1 - Y2) * r12) / (1 - r12 ** 2)
        d = np.linalg.norm(C12 - C23)
        temp = ((rho12 + rho23)**2 - d**2) * (d**2 - (rho12 - rho23)**2)
        # print(Y1, Y2, Y3, w1, w2, w3, rho12, rho23, C12, C23, d)
        K = 1/4 * sqrt(temp)

        x1_plus = 1/2 * (C12[0] + C23[0]) + 1/2 * (C23[0] - C12[0])/(d**2) * (rho12**2 - rho23**2) + 2*(C23[1]-C12[1]) * K/(d**2)
        x2_minus = 1/2 * (C12[1] + C23[1]) + 1/2 * (C23[1] - C12[1])/(d**2) * (rho12**2 - rho23**2) - 2*(C23[0]-C12[0]) * K/(d**2)
        x1_minus = 1/2 * (C12[0] + C23[0]) + 1/2 * (C23[0] - C12[0])/(d**2) * (rho12**2 - rho23**2) - 2*(C23[1]-C12[1]) * K/(d**2)
        x2_plus = 1/2 * (C12[1] + C23[1]) + 1/2 * (C23[1] - C12[1])/(d**2) * (rho12**2 - rho23**2) + 2*(C23[0]-C12[0]) * K/(d**2)

        X1 = np.array([x1_plus, x2_minus])
        X2 = np.array([x1_minus, x2_plus])

        # Choose the solution in the convex hull of the 3 points
        hull_original = ConvexHull(np.array([Y1, Y2, Y3]))
        hull_X1 = ConvexHull(np.array([Y1, Y2, Y3, X1]))
        # hull_X2 = ConvexHull(np.array([Y1, Y2, Y3, X2]))

        if np.all(sorted(hull_original.vertices) == sorted(hull_X1.vertices)):
            X = X1
        else:
            X = X2
    return X


def one_center_prob_w(points, weights):
    if len(points) == 1:
        return points[0]
    # Step 1
    n = 2
    if np.all(weights == 1):  # Unweighted case
        w1, w2 = 1, 1
        Y1, Y2 = None, None
        largest_dist = 0
        for i, point1 in enumerate(points[:-1]):
            for point2 in points[i+1:]:
                dist = np.linalg.norm(point1 - point2)
                if dist > largest_dist:
                    largest_dist = dist
                    Y1, Y2 = point1, point2

    else:
        i1, i2 = np.random.choice(np.arange(len(points)), 2, replace=False)
        Y1, Y2, w1, w2 = points[i1], points[i2], weights[i1], weights[i2]

    if w1 > w2:
        w1, w2 = w2, w1
        Y1, Y2 = Y2, Y1

    X = (w1 * Y1 + w2 * Y2) / (w1 + w2)
    RbX = np.linalg.norm(Y1 - X) * w1

    while True:
        if n == 2:
            # Step 2
            w_dists_to_center = np.linalg.norm(points - X, axis=1) * weights
            i_furthest = np.argmax(w_dists_to_center)
            Y3, w3 = points[i_furthest], weights[i_furthest]
            RX = w_dists_to_center[i_furthest]
            if RX - RbX < ZERO:
                return X

            # Step 3
            combinations = [(Y1, Y3), (Y2, Y3), (Y1, Y2, Y3)]
            not_in_combs = [Y2, Y1]
            weight_combs = [(w1, w3), (w2, w3), (w1, w2, w3)]
            not_in_w_combs = [w2, w1]
            for i in range(len(combinations)):
                if len(combinations[i]) == 2:
                    point1, point2 = combinations[i]
                    point3 = not_in_combs[i]
                    weight1, weight2 = weight_combs[i]
                    weight3 = not_in_w_combs[i]
                    X = (weight1 * point1 + weight2 * point2) / (weight1 + weight2)
                    wdist12 = np.linalg.norm(X - point1) * weight1
                    wdist3 = np.linalg.norm(X - point3) * weight3

                    if wdist3 - wdist12 < ZERO:
                        n = 2
                        Y1, Y2 = point1, point2
                        w1, w2 = weight1, weight2
                        break
                else:
                    X = center_of_3([Y1, Y2, Y3], [w1, w2, w3])
                    n = 3

            RbX = np.linalg.norm(Y1 - X) * w1

        if n == 3:
            # Step 4
            w_dists_to_center = np.linalg.norm(points - X, axis=1) * weights
            i_furthest = np.argmax(w_dists_to_center)
            Y4, w4 = points[i_furthest], weights[i_furthest]
            RX = w_dists_to_center[i_furthest]
            if RX - RbX < ZERO:
                return X

            # Step 5
            combinations = [(Y1, Y4), (Y2, Y4), (Y3, Y4), (Y1, Y2, Y4), (Y1, Y3, Y4), (Y2, Y3, Y4)]
            not_in_combs = [(Y2, Y3), (Y1, Y3), (Y1, Y2), Y3, Y2, Y1]
            weight_combs = [(w1, w4), (w2, w4), (w3, w4), (w1, w2, w4), (w1, w3, w4), (w2, w3, w4)]
            not_in_w_combs = [(w2, w3), (w1, w3), (w1, w2), w3, w2, w1]

            for i in range(len(combinations)):
                if len(combinations[i]) == 2:
                    point1, point2 = combinations[i]
                    point3, point4 = not_in_combs[i]
                    weight1, weight2 = weight_combs[i]
                    weight3, weight4 = not_in_w_combs[i]
                    X = (weight1 * point1 + weight2 * point2) / (weight1 + weight2)
                    wdist12 = np.linalg.norm(X - point1) * weight1
                    wdist3 = np.linalg.norm(X - point3) * weight3
                    wdist4 = np.linalg.norm(X - point4) * weight4
                    if wdist3 - wdist12 < ZERO and wdist4 - wdist12 < ZERO:  # Circle contains all 4 points
                        n = 2
                        Y1, Y2 = point1, point2
                        w1, w2 = weight1, weight2
                        break

                else:  # len(combinations[i] == 3:
                    X = center_of_3(combinations[i], weight_combs[i])
                    wdist_circ = np.linalg.norm(combinations[i][0] - X) * weight_combs[i][0]
                    wdist_last = np.linalg.norm(not_in_combs[i] - X) * not_in_w_combs[i]
                    if wdist_last - wdist_circ < ZERO:
                        Y1, Y2, Y3 = combinations[i]
                        w1, w2, w3 = weight_combs[i]
                        n = 3
                        break

            RbX = np.linalg.norm(Y1 - X) * w1


def p_center_prob_w(points, p, weights, init_X=None, print_plot=False, print_last=False):
    assert len(points) > p, "Number of facilities exceeds the number of points"
    colors = ["red", "blue", "yellow", "pink", "brown", "grey", "cyan", "purple", "orange", "black"]
    it = 0
    # Initial locations randomly in convex hull
    if init_X is None:
        # kmeans = KMeans(p).fit(points)
        # X = kmeans.cluster_centers_
        IDXS = np.random.choice(np.arange(len(points)), p, replace=False)
        X = points[IDXS]
    else:
        X = init_X

    if print_plot:
        plt.plot(*points.T, "ro")
        plt.plot(*X.T, "go")
        plt.show()

    max_change = 1
    while it == 0 or max_change > ZERO:
        # Assign to closest facility
        points_split = [[] for _ in range(p)]
        weights_split = [[] for _ in range(p)]
        for point, weight in zip(points, weights):
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

        # Solving one-center
        new_X = [None for _ in range(p)]
        max_d = 0
        max_wd = 0
        for i in range(p):
            if len(points_split[i]) > 0:  # If cluster has no points move facility to a random new location
                location = one_center_prob_w(points_split[i], weights_split[i])
                new_X[i] = location
                d = np.max(np.linalg.norm(location - points_split[i], axis=1))
                wd = np.max(np.linalg.norm(location - points_split[i], axis=1) * weights_split[i])
                max_d = np.max([max_d, d])
                max_wd = np.max([max_wd, wd])
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

    return X, max_d, max_wd


# my_points = np.array([[0, 0], [9, 0], [6, 3]])
# my_weights = np.array([2, 2, 3])
#
# one_center_prob(my_points, my_weights)
