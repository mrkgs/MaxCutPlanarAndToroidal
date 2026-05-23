# Main script that verifies the solutions of GSet toroidal Max Cut problems, which are saved in ./sols
import os.path
import numpy as np
toroid_nums = np.array([11, 12, 13, 32, 33, 34, 48, 49, 50, 57, 62, 65, 66, 67, 72, 77, 81])
nums = toroid_nums
from file_lib import load_gset_problem

if __name__ == '__main__':
    for GNum in nums:
        H, edges = load_gset_problem(GNum)
        mat_sz = H.shape[0]
        cycle_size = []
        for cycle_start in range(0, mat_sz):
            ii_connected = np.argwhere(H[cycle_start, :] != 0)
            ii_connected = ii_connected[ii_connected > cycle_start]
            if len(ii_connected) > 1:
                if ii_connected[1] < mat_sz/2:
                    cycle_size.append(int(ii_connected[1]) - cycle_start)
        D = int(np.round(np.mean(cycle_size)))
        # print('GSet #'+str(GNum)+', #Nodes='+str(mat_sz)+', D='+str(D))

        fname = "sols/gset_{:d}_maxcut_groups_hex.txt".format(GNum)
        if os.path.exists(fname):
            with open(fname, "r") as f:
                solutions = f.read()
                solutions_split = solutions.split('\n')
                solutions = []
                for sol in solutions_split:
                    if len(sol)>0:
                        solutions.append(sol)
                grp_mat = []
                Cuts = []
                for sol1 in solutions:
                    if len(sol1)>0:
                        grp_vec = []
                        for x in sol1:
                            y = format(int(x, 16), '04b')
                            for ii in range(4):
                                grp_vec.append(int(y[ii]))
                        grp_vec_c = [1-x for x in grp_vec]
                        new_sol = True
                        for grp_prev in grp_mat:
                            if grp_prev==grp_vec or grp_prev==grp_vec_c:
                                new_sol = False
                        if new_sol:
                            grp_mat.append(grp_vec)
                            Cut = 0
                            for edge in edges:
                                if grp_vec[edge[0]] != grp_vec[edge[1]]:
                                    Cut += int(H[edge[0], edge[1]])
                            Cuts.append(Cut)
                print("GSet #{:d} #Nodes={:5d} D={:3d} [{:2d}/{:2d}]".format(GNum,mat_sz,D, len(Cuts), len(solutions)), Cuts)


