import math
import numpy as np
from itertools import product, combinations, chain

class model_functions:


    # computes the denominator of the variance of the S1 estimate
    def compute_s1_var_coeff_c(self, breakpoint_, sampled_pts):
        sum_minus = 0
        num_minus = 0

        for x_i in sampled_pts[0]:
            if x_i > breakpoint_:
                sum_minus += x_i
                num_minus += 1

        if num_minus == 0:
            print(sampled_pts[0])
            print(breakpoint_)

        mean_minus = sum_minus / num_minus
        
        D_1 = 0
        for x_i in sampled_pts[0]:
            if x_i > breakpoint_:
                D_1 += (x_i - mean_minus) ** 2

        return(D_1)
        
    def compute_s2_var_coeff_c(self, breakpoint_, sampled_pts):
        
        alpha, beta_1, beta_2 = self.get_coeffs(breakpoint_, sampled_pts.T)\

        D_1 = self.compute_s1_var_coeff_c(breakpoint_, sampled_pts)

        sum_plus = 0
        num_plus = 0
        for x_i in sampled_pts[0]:
            if x_i <= breakpoint_:
                sum_plus += x_i
                num_plus += 1
        
        mean_plus = sum_plus / num_plus

        D_2 = 0
        for x_i in sampled_pts[0]:
            if x_i <= breakpoint_:
                D_2 += (x_i - mean_plus) ** 2
        if D_1 == 0:
            return D_2
        if D_2 == 0:
            return D_1
        
        L = (D_1 + D_2) / (D_1 * D_2)

        return(L)

    # given a set of points and an x-value of the breakpoint
    # estimate the alpha (y-int of clearance phase), 
    # beta_1 (slope in clearance phase), 
    # and beta_2 (difference between slope in prolif phase and beta_1) parameters of the underlying trajectory
    def get_coeffs(self, breakpoint_, points):

        # create array of all points above the breakpoint (in the clearance phase)
        above_row_mask = points[:, 0] > breakpoint_  
        above_break = points[above_row_mask]

        # get the mean x and y values of points above the breakpoint
        x_mean = sum(above_break[:,0]) / len(above_break[:,0])
        y_mean = sum(above_break[:,1]) / len(above_break[:,1])
        
        num = 0
        denom = 0
        for point in above_break:
            num += (point[0] - x_mean) * (point[1] - y_mean)
            denom += (point[0] - x_mean) ** 2
        
        beta_1 = (num / denom) # estimate the slope after the breakpoint

        alpha = y_mean - (beta_1 * x_mean)  # estimate the y-int of the line after the breakpoint

        # create array of all points below the breakpoint (in the prolif phase)
        below_row_mask = points[:, 0] < breakpoint_
        below_break = points[below_row_mask]
        # below_break = np.vstack((below_break, [breakpoint_, breakpoint_ * beta_1 + alpha]))
        
        # get the mean x and y values of points below the breakpoint
        x_mean = sum(below_break[:,0]) / len(below_break[:,0])
        y_mean = sum(below_break[:,1]) / len(below_break[:,1])
        
        num = 0
        denom = 0
        for point in below_break:
            num += (point[0] - x_mean) * (point[1] - y_mean)
            denom += (point[0] - x_mean) ** 2
        
        beta_2 = num / denom - beta_1 # estimate beta_2 (difference between slope in prolif phase and beta_1) 

        return(alpha, beta_1, beta_2)


    def compute_sse_full(self, sampled_pts, c, alpha, beta_1, beta_2):
        sse = 0
        for point in sampled_pts:
            sse += (point[1] - self.predictor(alpha, beta_1, beta_2, c, [point[0]])[0]) ** 2

        return(sse)
    
    # given the slope and intercept estimates and the breakpoint
    # and a set of timestamps x
    # predicts the values of y
    def predictor(self, alpha, beta_1, beta_2, breakpoint_, x):
        y = []
        for i in range(len(x)):
            if x[i] >= breakpoint_:
                y_i = alpha + (beta_1 * x[i])
            else:
                y_i = alpha + (beta_1 * x[i]) + (beta_2 * (x[i] - breakpoint_))
            
            y.append(y_i)
        
        return(y)
        
    def optimal_coeffs(self, breakpoint_, sampled_pts, num_process):
        # split the points into those below and above the breakpoint
        row_mask = sampled_pts.T[:, 0] < breakpoint_
        below_break = sampled_pts.T[row_mask]
        above_break = sampled_pts.T[~row_mask]

        # generate a list of all pairs of integers which sum to num_process
        # num process is the number of sampled points we will actually process (i.e. use for estimation) 
        pairs = []
        for sub in range(2, num_process - 1):
            pair = [sub, num_process - sub]
            pairs.append(pair)

        # generate a list of all subsets of num_process points
        # which contain at least 1 point in the proliferation phase
        # and 2 in the clearance phase
        all_combos = []
        for pair in pairs:
            combinations_above = list(combinations(below_break, pair[0]))
            combinations_below = list(combinations(above_break, pair[1]))
            combinationsab = list(product(combinations_above, combinations_below))
            all_combos = all_combos + [combo[0] + combo[1] for combo in combinationsab]
        all_combos = np.array(all_combos)

        # find the subset that minimizes the sum of the variances of the slope estimates
        overall_min = 100000000
        for combo in all_combos:
            L1 = self.compute_s1_var_coeff_c(breakpoint_, combo.T)
            L2 = self.compute_s2_var_coeff_c(breakpoint_, combo.T)
            if (1 / L1) + L2 < overall_min:
                overall_min = (1 / L1) + L2
                points3 = combo 

        # estimate the coefficients
        alpha, beta_1, beta_2 = self.get_coeffs(breakpoint_, points3)

        return([alpha, beta_1, beta_2, points3]) 


    def estimate(self, sampled_pts, c_val_step, num_process):
        c_values = np.arange(sorted(sampled_pts[0])[1] + c_val_step, sorted(sampled_pts[0])[-2], c_val_step)

        all_Fc = []
        min_Fc = 1000000
        for c in c_values:
            alpha, beta_1, beta_2, pts = self.optimal_coeffs(c, sampled_pts, num_process)
            Fc = self.compute_sse_full(pts, c, alpha, beta_1, beta_2)
            if Fc < min_Fc:
                min_Fc = Fc
                changepoint = c
                coeffs = [alpha, beta_1, beta_2]
                optimal_pts = pts
        return(changepoint, coeffs, optimal_pts)