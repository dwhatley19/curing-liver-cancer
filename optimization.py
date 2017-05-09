import scipy.io as io
import numpy as np
import thread
import time
from gurobipy import *

# solving the problem for pre-set gantry and couch angles
# uncomment lines as needed
mat_contents = io.loadmat('TG119/Gantry0_Couch0_D.mat')
# mat_contents = io.loadmat('TG119/Gantry72_Couch0_D.mat')
# mat_contents = io.loadmat('TG119/Gantry144_Couch0_D.mat')
# mat_contents = io.loadmat('TG119/Gantry216_Couch0_D.mat')
# mat_contents = io.loadmat('TG119/Gantry288_Couch0_D.mat')

# modify to change start and end points of analysis
start = 0
end = mat_contents['D'].shape[0] - 1

oar_contents = io.loadmat('TG119/Core_VOILIST.mat')
oar_voxels = {}

tumor_contents = io.loadmat('TG119/OuterTarget_VOILIST.mat')
tumor_voxels = {}

frac_appl_oar = 0.0
for x in oar_contents['v']:
    if(x[0] > start and x[0] <= end):
        frac_appl_oar += 1
        oar_voxels[x[0]] = True

frac_appl_oar /= oar_contents['v'].size

for x in tumor_contents['v']:
    if(x[0] > start and x[0] <= end):
        tumor_voxels[x[0]] = True

# ACTUAL OPTIMIZATION CODE STARTS HERE
# do not modify anything unless otherwise specifically mentioned

#uncomment this line too
x = mat_contents['D'].shape

# UNCOMMENT THIS CODE WHEN OPTIMIZING
rows = end - start

m = Model("dose influence")
# variables, upper/lower bounds on x[i]
flu = [0] * x[1]
for i in range(x[1]):
    # modify 25 as needed
    flu[i] = m.addVar(lb=0, ub=25, name=('x'+str(i)))
fluence = np.array(flu)

arr = np.array(mat_contents['D'][start:end+1].todense())
nonzero = (np.sum(arr, axis=1) > 0).astype(int)
indices = np.cumsum(nonzero) - 1
arr = arr[np.sum(arr, axis=1) > 0]
d = [0] * (arr.shape[0])

# upper/lower bounds on d[i]
for i in range(arr.shape[0]):
    # modify 0.75 as needed
    d[i] = m.addVar(lb=0, ub=0.75, name=('d'+str(i)))

# dose-influence matrix stuff
t1 = time.time()
m.addConstrs((np.sum(arr[i] * fluence) == d[i] for i in range(arr.shape[0])))

# tumor region constraint
lb_tumor = 0.01 # modify as needed
for i in tumor_voxels:
    if(nonzero[i - start]):
        m.addConstr(d[indices[i - start]] >= lb_tumor)

# OAR region constraint on total sum
ub_OAR = 20 # modify as needed
oar_constr = 0
for i in oar_voxels:
    if(nonzero[i - start]):
        oar_constr += d[indices[i - start]]
m.addConstr(oar_constr <= ub_OAR)

# objective function
# sum of doses to all the voxels
obj = 0
for i in range(start, end):
    if(nonzero[i - start]):
        obj += d[indices[i - start]]
m.setObjective(obj, GRB.MAXIMIZE)

m.update()
# uncomment this line to print problem to "debug.lp"
# m.write("debug.lp")

# solve!
m.optimize()

total = 0

# uncomment this line to print total runtime
# print(time.time() - t1)

for v in m.getVars():
    if(v.varName[0] == 'x'): print(v.varName, v.x)
    if(v.varName[0] == 'd'): total += v.x

# compute average of all the OAR voxel doses
oar_total = 0
for i in range(start, end):
    if(i in oar_voxels and nonzero[i - start]):
        oar_total += d[indices[i - start]].x

# compute average of all the tumor voxel doses
tumor_total = 0
for i in range(start, end):
    if(i in tumor_voxels and nonzero[i - start]):
        tumor_total += d[indices[i - start]].x

print('Obj:', m.objVal)
print('Average tumor:', tumor_total / len(tumor_voxels))
print('Average OAR:', oar_total / len(oar_voxels))