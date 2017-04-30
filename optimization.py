import scipy.io as io
import numpy as np
from gurobipy import *

start = 1600000
end = 1601000

voi_contents = io.loadmat('TG119/Body_VOILIST.mat')
voxels = {}

for x in voi_contents['v']:
	if(x[0] > start and x[0] <= end):
		voxels[x[0]-start] = True

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
    fluence[i] = m.addVar(lb=0, name=('x'+str(i)))

# constraint 2, 3, 4
d = [0] * (end - start)
for i in range(start, end):
    d[i - start] = m.addVar(lb=0, ub=100, name=('d'+str(i)))
    arr = mat_contents['D'].getrow(i).toarray()
    m.addConstr(quicksum(arr[0,j] * fluence[j] for j in range(x[1])) == d[i - start])

# constraint 6 not needed
# because we are pre-setting the beam angle

# constraint 7
ub_ntumor = 75
for i in range(start, end):
	if(i - start not in voxels):
		m.addConstr(ub_ntumor >= d[i - start])

# constraint 8 not needed
# because not implementing constraint 1

# objective function
# leave alpha at 0, no need to consider y's
obj = 0
for i in voxels:
	obj += d[i]
m.setObjective(obj, GRB.MAXIMIZE)

m.update()
m.write("debug.lp")

# solve!
m.optimize()

for v in m.getVars():
    print(v.varName, v.x)

print('Obj:', m.objVal)