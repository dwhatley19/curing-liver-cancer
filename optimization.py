import scipy.io as io
import numpy as np
from gurobipy import *

start = 1519000
end = 1520000

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

# solving the problem for pre-set gantry and couch angles
# uncomment lines as needed
mat_contents = io.loadmat('TG119/Gantry0_Couch0_D.mat')
# mat_contents = io.loadmat('TG119/Gantry72_Couch0_D.mat')
# mat_contents = io.loadmat('TG119/Gantry144_Couch0_D.mat')
# mat_contents = io.loadmat('TG119/Gantry216_Couch0_D.mat')
# mat_contents = io.loadmat('TG119/Gantry288_Couch0_D.mat')
# f = open('00.txt', 'w')

# sz = mat_contents['v'].size
# for i in range(sz):
#   if(mat_contents['v'][i].size > 1):
#       print(mat_contents['v'][i].size)

# this prints the max, average, and min doses when
#   fluencies are all 1

# ACTUAL OPTIMIZATION CODE STARTS ON LINE 30

#uncomment this line too
x = mat_contents['D'].shape

# doses = []
# for i in range(1600000, 1601000):
#   if(mat_contents['D'][i].size == 0): continue
#   r = np.dot(mat_contents['D'].getrow(i).todense(), np.ones(x[1]))
#   doses.append(r[0])

# print("Max dose:", max(doses))
# print("Avg dose:", sum(doses)/len(doses))
# print("Min dose:", min(doses))

#TODO:
# move around lower bounds
# introduce constraints such as "you can go lower than lower bounds
#   but that'll induce some suboptimality"
# 

# UNCOMMENT THIS CODE WHEN OPTIMIZING
rows = end - start

m = Model("dose influence")
# variables, constraint 5
fluence = [0] * x[1]
for i in range(x[1]):
    fluence[i] = m.addVar(lb=0, ub=2500, name=('x'+str(i)))

# constraint 2, 3, 4
d = [0] * (end - start)
for i in range(start, end):
    d[i - start] = m.addVar(lb=0, ub=1, name=('d'+str(i)))
    arr = mat_contents['D'].getrow(i).toarray()
    m.addConstr(quicksum(arr[0,j] * fluence[j] for j in range(x[1])) == d[i - start])

# constraint 6 not needed
# because we are pre-setting the beam angle

# tumor region constraint
lb_tumor = 0.65
for i in tumor_voxels:
    m.addConstr(d[i - start] >= lb_tumor)

# OAR region constraint on total sum
ub_OAR = 20 * frac_appl_oar
print(ub_OAR)
oar_constr = 0
for i in oar_voxels:
    oar_constr += d[i - start]
m.addConstr(oar_constr <= ub_OAR)

# constraint 8 not needed
# because not implementing constraint 1

# objective function
# leave alpha at 0, no need to consider y's
obj = 0
for i in range(start, end):
    obj += d[i - start]
m.setObjective(obj, GRB.MAXIMIZE)

m.update()
m.write("debug.lp")

# solve!
m.optimize()

total = 0

for v in m.getVars():
    print(v.varName, v.x)
    if(v.varName[0] == 'd'): total += v.x

print('Obj:', m.objVal)
print('Average tumor:', m.objVal / len(tumor_voxels))
print('Average OAR:', (total - m.objVal) / (end - start - len(oar_voxels)))