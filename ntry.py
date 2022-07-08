from pymoo.factory import get_visualization, get_reference_directions
from pymoo.factory import get_performance_indicator
import numpy as np
from numpy import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
from pandas.plotting import scatter_matrix
import pandas as pd
# %matplotlib widget

def rn():
    return random.rand()
def random_cord(dim):
    rnums = np.random.rand(dim-1)
    rnums = np.append(rnums, np.array([0,1]))
    rnums.sort()
    diffs = []
    for i in range(len(rnums)-1):
        diffs.append(rnums[i+1] - rnums[i])
    return diffs

def random_sampling(M, N):
    nums = []
    for i in range(N):
        nums.append(random_cord(M))
    hv = get_performance_indicator("hv", ref_point=np.array([1.01] * M))
    # print("Effectiveness:", round(hv.do(np.array(nums)) * 100, 2), "%")
    return round(hv.do(np.array(nums)) * 100, 2)

def convergance():
    X = [] 
    Y = []
    for i in range(10, 1000, 10):
        metric = 0
        for j in range(10):
            metric += random_sampling(3, i)
        metric /= 10
        print("Round", i, "Results:", metric)
        X.append(i)
        Y.append(metric)
    plt.style.use('seaborn-whitegrid')
    fig = plt.figure()
    ax = fig.add_subplot(111)  
    ax.plot(X, Y)
    plt.show()
# convergance()

def dirichlet(a, m):
    return np.random.dirichlet([a] * m)
    
# dirichlet(1, 3)

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.plot([0,0,1, 0], [0,1,0, 0], [1,0,0,1], c='red')
# img = ax.scatter(rd[:,0], rd[:,1], rd[:,2], c=rd[:,3], cmap=plt.hot())
points = []
for i in range(100):    
    v = dirichlet(1, 3)
    points.append(np.array(v))
p = np.array(points)
cm = plt.cm.get_cmap('winter')
g = ax.scatter(p[:,0], p[:,1], p[:,2], c=p[:,2], cmap=cm)
fig.colorbar(g)
ax.view_init(30,30)
plt.show()

    
    
