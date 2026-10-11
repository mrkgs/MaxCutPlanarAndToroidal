# Parameters

# Toroidal graphs : [11, 12, 13, 32, 33, 34, 48, 49, 50, 57, 62, 65, 66, 67, 72, 77, 81]
gsetnum = [11, 12, 13, 32, 33, 34, 48, 49, 50, 57, 62, 65, 66, 67, 72, 77, 81]
# gsetnum = [11, 12, 13, 32, 33, 34]
# gsetnum = [11]               # Comment this line to run all problems
use_multi_proc = 0
ncycles = 24
MinPathThreshhold = 1E6
MinPathThreshhold = 5        # Comment this line to work with full input matrix for MWM

alg_num = 5                  # Sum of the algorithms to run: 1 - Algorithm#1, 2 - Algorithm#2, 4 - Liers-Pardella (LP)
                             # 3 - Algorithms #1 and #2, 5 - Algorithm#1 and LP, 6 - Algorithm#2 and LP, 7 - all
lp_max_iters = 50            # LP: max number of improvement iterations
lp_max_stall = 10            # LP: stop after this number of iterations in a row without improvement
