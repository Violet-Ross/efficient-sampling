import math
import numpy as np
import itertools
import piecewise_regression

class testing_functions:

    # stochastically generates a ground truth trajectory 
    def generate_traj(self, num_traj):
        t = []
        vol = []
        trajec_num = []

        for i in range(num_traj):
            t_0 = np.random.uniform(2.5, 3.5)
            t.append(t_0)

            candidate_w_1 = 100
            while candidate_w_1 > 2.5:
                candidate_w_1 = np.random.gamma(1, 1.5)
            t.append(t_0 + candidate_w_1 + 0.5)

            w_2 = np.random.uniform(4, 9)
            t.append(t_0 + candidate_w_1 + 0.5 + w_2)

            v_peak = np.random.uniform(7, 11)
            vol.extend([3, v_peak, 6])

            trajec_num.extend([i, i, i])
        
        df = np.array([t, vol, trajec_num])
        return(df)
    
    # given a trajectory traj and set of timestamps x
    # generates the true value of y at each x
    def get_y(self, traj, x):
        m_1 = (traj[1,1] - traj[1,0]) / (traj[0,1] - traj[0,0]) # slope of first line seg
        b_1 = traj[1,1] - (m_1 * traj[0,1]) # y-int of first line seg

        m_2 = (traj[1,2] - traj[1,1]) / (traj[0,2] - traj[0,1]) # slope of second line seg
        b_2 = traj[1,2] - (m_2 * traj[0,2]) # y-int of second line seg

        y = []
        for i in range(len(x)):
            if x[i] <= traj[0, 1]:
                y_i = (m_1 * x[i]) + b_1
            else:
                y_i = (m_2 * x[i]) + b_2
            
            y.append(y_i)
        
        return(y)
    
    # samples num_samples points from the trajectory traj
    def sample_traj(self, traj, num_samples):
        x = np.linspace(traj[0, 0], traj[0, 2], num_samples)

        y = self.get_y(traj, x)

        sampled_pts = np.array([x, y])

        return(sampled_pts.T)
    
class fitting_functions:
    # given a set of points and an x-value of the breakpoint
    # estimate the alpha (y-int of clearance phase), 
    # beta_1 (slope in clearance phase), 
    # and beta_2 (difference between slope in prolif phase and beta_1) parameters of the underlying trajectory
    def get_coeffs(self, breakpoint_, points):

        # create array of all points below the breakpoint (in the prolif phase)
        below_row_mask = points[:, 0] <= breakpoint_  
        below_break = points[below_row_mask]

        # get the mean x and y values of points below the breakpoint
        x_mean = sum(below_break[:,0]) / len(below_break[:,0])
        y_mean = sum(below_break[:,1]) / len(below_break[:,1])
        
        num = 0
        denom = 0
        for point in below_break:
            num += (point[0] - x_mean) * (point[1] - y_mean)
            denom += (point[0] - x_mean) ** 2
        
        beta_1 = (num / denom) # estimate the slope below the breakpoint

        alpha = y_mean - (beta_1 * x_mean)  # estimate the y-int of the line before the breakpoint

        # create array of all points above the breakpoint (in the clearance phase)
        above_row_mask = points[:, 0] >= breakpoint_
        above_break = points[above_row_mask]
        
        # get the mean x and y values of points above the breakpoint
        x_mean = sum(above_break[:,0]) / len(above_break[:,0])
        y_mean = sum(above_break[:,1]) / len(above_break[:,1])
        
        num = 0
        denom = 0
        for point in above_break:
            num += (point[0] - x_mean) * (point[1] - y_mean)
            denom += (point[0] - x_mean) ** 2
        
        beta_2 = num / denom - beta_1 # estimate beta_2 (difference between slope in clearance phase and beta_1) 

        return(alpha, beta_1, beta_2)
    
    # given the slope and intercept estimates and the breakpoint
    # and a set of timestamps x
    # predicts the values of y
    def predictor(self, alpha, beta_1, beta_2, breakpoint_, x):
        y = []
        for i in range(len(x)):
            if x[i] <= breakpoint_:
                y_i = alpha + (beta_1 * x[i])
            else:
                y_i = alpha + (beta_1 * x[i]) + (beta_2 * (x[i] - breakpoint_))
            
            y.append(y_i)
        
        return(y)
    
    # computes the denominator of the variance of the S1 estimate
    def compute_s1_var_coeff_c(self, breakpoint_, sampled_t):
        sum_minus = 0
        num_minus = 0

        for x_i in sampled_t:
            if x_i <= breakpoint_:
                sum_minus += x_i
                num_minus += 1

        if num_minus == 0:
            print(sampled_t)
            print(breakpoint_)

        mean_minus = sum_minus / num_minus
        
        D_1 = 0
        for x_i in sampled_t:
            if x_i <= breakpoint_:
                D_1 += (x_i - mean_minus) ** 2

        return(D_1)
    
    def compute_s2_var_coeff_c(self, breakpoint_, sampled_t):
        D_1 = self.compute_s1_var_coeff_c(breakpoint_, sampled_t)

        sum_plus = 0
        num_plus = 0
        for x_i in sampled_t:
            if x_i >= breakpoint_:
                sum_plus += x_i
                num_plus += 1
        
        mean_plus = sum_plus / num_plus

        D_2 = 0
        for x_i in sampled_t:
            if x_i >= breakpoint_:
                D_2 += (x_i - mean_plus) ** 2
        if D_1 == 0:
            return D_2
        if D_2 == 0:
            return D_1
        
        #L = (D_1 + D_2) / (D_1 * D_2)

        return(D_2)
        
    def compute_sse(self, sampled_pts, c, alpha, beta_1, beta_2):
        sse = 0
        for point in sampled_pts:
            sse += (point[1] - self.predictor(alpha, beta_1, beta_2, c, [point[0]])[0]) ** 2

        return(sse)
    
    def estimate_bp(self, taken_points, c_val_step):
        c_values = np.arange(sorted(taken_points[:, 0])[1] + c_val_step, sorted(taken_points[:, 0])[-2], c_val_step)

        min_Fc = 1000000
        for c in c_values:
            alpha, beta_1, beta_2 = self.get_coeffs(c, taken_points)
            Fc = self.compute_sse(taken_points, c, alpha, beta_1, beta_2)
            if Fc < min_Fc:
                min_Fc = Fc
                changepoint = c
                coeffs = [alpha, beta_1, beta_2]
        return(coeffs[1], coeffs[2], changepoint)
    
    def ross_heuristic(self, sampled_pts, num_iter, c_val_step):
        all_coeff_ests = []

        sorted_pts = sorted(sampled_pts, key=lambda sampled_pt : sampled_pt[0])
        initial_pts = np.array([sorted_pts[0], sorted_pts[len(sorted_pts) // 3], sorted_pts[(len(sorted_pts) // 3) * 2], sorted_pts[-1]])
        beta_1, beta_2, cp = self.estimate_bp(initial_pts, c_val_step)
        all_coeff_ests.append([beta_1, beta_2, cp])

        taken_t = initial_pts[:,0]

        overall_min = 100000
        for i in range(num_iter):
            combos = [np.append(taken_t, point) for point in sampled_pts[:,0] if point not in taken_t]  
            overall_min = 100000
            for combo in combos:
                L1 = self.compute_s1_var_coeff_c(cp, combo)
                L2 = self.compute_s2_var_coeff_c(cp, combo)
                if L1 + L2 < overall_min:
                    overall_min = L1 + L2
                    best_combo = combo 
            taken_t = best_combo
            taken_points = np.array([point for point in sampled_pts if point[0] in taken_t])
            beta_1, beta_2, cp = self.estimate_bp(taken_points, c_val_step)
            all_coeff_ests.append([beta_1, beta_2, cp])
        
        return all_coeff_ests

    def optimal_pts_approach(self, sampled_pts, num_iter, c_val_step, true_coeffs):
        all_coeff_ests = []
        sorted_pts = sorted(sampled_pts, key=lambda sampled_pt : sampled_pt[0])
        initial_pts = np.array([sorted_pts[0], sorted_pts[len(sorted_pts) // 3], sorted_pts[(len(sorted_pts) // 3) * 2], sorted_pts[-1]])
        beta_1, beta_2, cp = self.estimate_bp(initial_pts, c_val_step)
        all_coeff_ests.append([beta_1, beta_2, cp])

        taken_t = initial_pts[:,0]

        for i in range(num_iter):
            combos = [np.append(taken_t, point) for point in sampled_pts[:,0] if point not in taken_t]  
            overall_min = 100000
            for combo in combos:
                combo_points = np.array([point for point in sampled_pts if point[0] in combo])
                beta_1, beta_2, cp = self.estimate_bp(combo_points, c_val_step)
                accuracy = math.dist(true_coeffs, [beta_1, beta_2, cp])
                if accuracy < overall_min:
                    best_combo = combo
                    overall_min = accuracy
            taken_t = best_combo
            taken_points = np.array([point for point in sampled_pts if point[0] in taken_t])
            beta_1, beta_2, cp = self.estimate_bp(taken_points, c_val_step)
            all_coeff_ests.append([beta_1, beta_2, cp])
        
        return all_coeff_ests
    
    def random_pts_approach(self, sampled_pts, num_iter, c_val_step):
        all_coeff_ests = []
        sorted_pts = sorted(sampled_pts, key=lambda sampled_pt : sampled_pt[0])
        initial_pts = np.array([sorted_pts[0], sorted_pts[len(sorted_pts) // 3], sorted_pts[(len(sorted_pts) // 3) * 2], sorted_pts[-1]])
        beta_1, beta_2, cp = self.estimate_bp(initial_pts, c_val_step)
        all_coeff_ests.append([beta_1, beta_2, cp])

        taken_t = initial_pts[:,0]
        for i in range(num_iter):
            combos = np.array([np.append(taken_t, point) for point in sampled_pts[:,0] if point not in taken_t])
            random_row_index = np.random.choice(combos.shape[0])
            taken_t = combos[random_row_index]
            taken_points = np.array([point for point in sampled_pts if point[0] in taken_t])
            beta_1, beta_2, cp = self.estimate_bp(taken_points, c_val_step)
            all_coeff_ests.append([beta_1, beta_2, cp])
        
        return all_coeff_ests
    
    def simple_varmin_approach(self, sampled_pts, num_iter, c_val_step):
        all_coeff_ests = []
        sorted_pts = sorted(sampled_pts, key=lambda sampled_pt : sampled_pt[0])
        initial_pts = np.array([sorted_pts[0], sorted_pts[len(sorted_pts) // 3], sorted_pts[(len(sorted_pts) // 3) * 2], sorted_pts[-1]])
        beta_1, beta_2, cp = self.estimate_bp(initial_pts, c_val_step)
        all_coeff_ests.append([beta_1, beta_2, cp])

        for i in range(num_iter):
            combos = list(itertools.combinations(range(len(sampled_pts)), 4 + i))

            maximum = 0
            for combo in combos:
                points = np.array([sampled_pts[j] for j in range(len(sampled_pts)) if j in combo])
                xbar = sum(points[:,0]) / len(points[:,0])
                varcoef = sum((points[:,0] - xbar) ** 2)
                if varcoef > maximum:
                    optimal_pts = points
                    maximum = varcoef
            pw_fit = piecewise_regression.Fit(optimal_pts[:,0], optimal_pts[:,1], n_breakpoints=1)
            if pw_fit.get_params() == {'converged' : False}:
                all_coeff_ests.append([np.nan, np.nan, np.nan])
            else:
                alpha, beta1, beta2 = self.get_coeffs(pw_fit.get_params()['breakpoint1'], sampled_pts)
                all_coeff_ests.append([beta1, beta2, pw_fit.get_params()['breakpoint1']])

        return all_coeff_ests