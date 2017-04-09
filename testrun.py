import scipy.io as io
import numpy as np

# mat_contents = io.loadmat('TG119/Body_VOILIST.mat')
mat_contents = io.loadmat('TG119/Gantry0_Couch0_D.mat')
f = open('00.txt', 'w')

# sz = mat_contents['v'].size
# for i in range(sz):
# 	if(mat_contents['v'][i].size > 1):
# 		print(mat_contents['v'][i].size)

# this prints the max, average, and min doses when
#	fluencies are all 1
x = mat_contents['D'].shape
doses = []
print(mat_contents['D'].getrow(1602471).todense(), np.ones(x[1]))
print(np.dot(mat_contents['D'].getrow(1602471).todense(), np.ones(x[1])))
for i in range(1600000, 1601000):
	if(mat_contents['D'][i].size == 0): continue
	r = np.dot(mat_contents['D'].getrow(i).todense(), np.ones(x[1]))
	doses.append(r[0])

print("Max dose:", max(doses))
print("Avg dose:", sum(doses)/len(doses))
print("Min dose:", min(doses))