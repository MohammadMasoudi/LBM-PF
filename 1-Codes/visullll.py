import pickle
import numpy as np
import matplotlib.pyplot as plt
from skimage.measure import shannon_entropy

time = [i for i in range(1,8201)]
entropy= np.zeros([2,8200])
qq=-1
qq+=1

for i in time:
    with open('C:/Users/masoudi/Desktop/privat/Papers/Nature/Results/python/Sensivity-9 cases_run2/results_3/results_%s.pkl'%(i), 'rb') as f:
        struc = pickle.load(f)

    struc =struc .astype(int)    
    entropy[qq,i-1]=shannon_entropy(struc)
    
qq+=1
for i in time:
    with open('C:/Users/masoudi/Desktop/privat/Papers/Nature/Results/python/second duplicaterun/results_3/results_%s.pkl'%(i), 'rb') as f:
        struc = pickle.load(f)

    struc =struc .astype(int)    
    entropy[qq,i-1]=shannon_entropy(struc)

fig, ax = plt.subplots()

# for j in B:
#     # plt.figure()
    # plt.plot(time, entropy[0,:],label='case %s'%j)

ax.scatter(x=time, y=entropy[0,:])
ax.scatter(x=time, y=entropy[1,:])
# ax.scatter(x=time, y=entropy[2,:],label='case 7')
# ax.scatter(x=time, y=entropy[3,:],label='case 8')
# ax.scatter(x=time, y=entropy[4,:],label='case 11')
# ax.scatter(x=time, y=entropy[5,:],label='case 12')

ax.set_xlabel("time")
ax.set_ylabel("entropy")
ax.set_title("t-e")    
ax.legend()