import numpy as np
import matplotlib.pyplot as plt

ZERO = 0.0001


class LocPlanProb2d:
    def __init__(self, points=None, weights=None):
        if points is None:
            points = np.empty((0, 2))
        self.points = np.array(points)
        if weights is None:
            weights = [1] * len(points)
        self.weights = np.array(weights)
        self.sol = None
        self.sol_hist = []
        self.obj = None
        self.obj_hist = []
        self.iterations = 1

    def add_point(self, point, weight):
        self.points = np.vstack([self.points, point])
        self.weights = np.append(self.weights, weight)

    def set_sol(self, point):
        self.sol = point

    def init_sol(self):
        self.sol = 1 / np.sum(self.weights) * self.weights @ self.points
        self.sol_hist.append(self.sol)

    def update_sol_weiszfeld(self, hist=True):
        d = np.linalg.norm(self.sol - self.points, axis=1)
        if hist:
            self.obj = self.weights @ d
            self.obj_hist.append(self.obj)

        sX = np.sum(self.weights / d)
        alpha = (1 / sX) * (self.weights / d)

        self.sol = alpha @ self.points
        if hist:
            self.sol_hist.append(self.sol)

    def solve_weiszfeld(self):
        self.init_sol()
        self.update_sol_weiszfeld()
        while np.linalg.norm(self.sol - self.sol_hist[-2]) > ZERO:
            self.update_sol_weiszfeld()
            self.iterations += 1

    def update_sol_ostresh(self, lam, hist=True):
        d = np.linalg.norm(self.sol - self.points, axis=1)
        if hist:
            self.obj = self.weights @ d
            self.obj_hist.append(self.obj)
        if np.min(d) > ZERO ** 2:
            lam_val = lam
            if lam == 'drezner':  # Only used if sol is not on a customer point (else lambda = 1)
                lam_val = self.drezner_lam()
            sX = np.sum(self.weights / d)
            gradient = np.sum(self.weights / d * (self.sol - self.points).T, axis=1)
            self.sol = self.sol - lam_val / sX * gradient
        else:
            k = [i for i, point in enumerate(self.points) if np.linalg.norm(point - self.sol) < ZERO][0]
            w_not_k = np.concatenate((self.weights[:k], self.weights[k+1:]))
            sXk = np.sum(w_not_k / d[d >= ZERO])
            X_not_k = np.concatenate((self.sol[:k], self.sol[k+1:]))
            points_not_k = np.concatenate((self.points[:k], self.points[k+1:]))
            gradient_k = np.sum(w_not_k / d[d >= ZERO] * (X_not_k - points_not_k).T, axis=1)
            norm_grad = np.linalg.norm(gradient_k)
            if ZERO < norm_grad < self.weights[k]:
                self.sol = self.sol - 1 / sXk * (1 - self.weights[k] / norm_grad) * gradient_k
        if hist:
            self.sol_hist.append(self.sol)

    def solve_ostresh(self, lam='drezner'):
        self.init_sol()
        self.update_sol_ostresh(lam)
        while np.linalg.norm(self.sol - self.sol_hist[-2]) > ZERO:
            self.update_sol_ostresh(lam)
            self.iterations += 1

    def drezner_lam(self):
        d = np.linalg.norm(self.sol - self.points, axis=1)
        Xprime = self.sol.copy()
        self.update_sol_ostresh(lam=1, hist=False)
        X2prime = self.sol.copy()
        self.set_sol(Xprime)
        d2 = np.linalg.norm(self.points - X2prime, axis=1)

        fprime = self.weights @ d
        f2prime = self.weights @ d2

        temp = np.sum(self.weights / d) * np.linalg.norm(Xprime - X2prime) ** 2
        if f2prime - fprime + temp > ZERO:
            lam = 0.5 * temp / (f2prime - fprime + temp)
        else:
            lam = 1

        if np.min(d) < 0.1:  # To avoid overshooting when close to customer
            lam = 1

        print(lam)
        return lam

    def plot_hist(self):
        plt.plot(*zip(*self.points), "ro")
        for sol in self.sol_hist:
            plt.plot(*sol, "go")
        plt.show()

# ins = LocPlanProb2d([[0,0], [30,0], [0,40]], [5,3,4])
# ins = LocPlanProb2d([[0,0], [1,0], [0,1], [1, 1], [100, 100]], [1, 1, 1, 1, 4])
