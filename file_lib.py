# Files read/write functions
import numpy as np
from IPython.display import clear_output

def load_gset_problem(GNum):
    clear_output(wait=True)
    dirname = 'data/'
    filename = 'G' + str(GNum) + '.txt'
    f = open(dirname + filename, "r").readlines()
    max_i, max_j = 0, 0
    for ele in f[1:]:
        i, j, v = ele.split()
        i = int(i)
        j = int(j)
        if i > max_i:
            max_i = i
        if j > max_j:
            max_j = j
    mat_sz = max(max_i, max_j)
    H = np.zeros((mat_sz, mat_sz))
    edges = []
    for ele in f[1:]:
        i, j, v = ele.split()
        i = int(i)
        j = int(j)
        v = int(v)
        H[i - 1, j - 1] = v
        H[j - 1, i - 1] = v
        edges.append([i - 1, j - 1])
    return H, edges

def write_bin_to_hex_file(dd, fname, ext="a"):
    dd = dd.astype(int)
    assert dd.shape[1] / 4 - dd.shape[1] // 4 == 0
    with open(fname, ext) as flog:
        for rcnt in range(dd.shape[0]):
            for ccnt in range(dd.shape[1] // 4):
                newchar = hex(int(''.join(map(str, dd[rcnt, ccnt * 4 + np.arange(4)])), 2))[2:]
                flog.write(f"{newchar}")
            flog.write(f"\n")

def read_from_hex_file_to_bin(fname):
    with open(fname, "r") as f:
        data_all = f.read()
        data_all = data_all.split('\n')
        dd = []
        for data1 in data_all:
            if len(data1)>0:
                d = []
                for x in data1:
                    y = format(int(x, 16), '04b')
                    for ii in range(4):
                        d.append(int(y[ii]))
        dd.append(d)
        return np.array(dd)
