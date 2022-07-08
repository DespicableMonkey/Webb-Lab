# %matplotlib widget

from pymoo.factory import get_visualization, get_reference_directions
from pymoo.factory import get_performance_indicator
import numpy as np
from numpy import random, prod
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pandas.plotting import scatter_matrix
import pandas as pd
from regex import F
from queue import Queue 
import matplotlib.cm as cm
from scipy.ndimage.filters import gaussian_filter
import math

''' Visualtion & Pymoo Methods '''
def getmap3D(viewInitX=30, viewInitY=45, title=""):
     fig,ax = plt.subplots(ncols=1,nrows=1,subplot_kw=dict(projection='3d'))
     ax.view_init(viewInitX, viewInitY)
     ax.set_title(title)
     return ax

def getmap2D():
    fig,ax=plt.subplots(ncols=1, nrows=1)
    # ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    return ax
def getmapdef2D():
    return plt.subplots(ncols=1, nrows=1)
def getmapaxdef2D():
    fig,ax= plt.subplots(ncols=1,nrows=1)
    return ax

def getLimitedMap3D(name=False,X=30,Y=45, title=""):
    ax = getmap3D(X, Y, title)
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ax.set_zlim(0,1)
    
    if(name):
       ax.set_xlabel('X axis')
       ax.set_ylabel('Y axis')
       ax.set_zlabel('Z axis') 
    return ax

def getheatmap(x, y, s, bins=1000):
    heatmap, xedges, yedges = np.histogram2d(x, y, bins=bins)
    heatmap = gaussian_filter(heatmap, sigma=s)

    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    return heatmap.T, extent
def drawheatmap(ax, X, Y, o):
    img, extent = getheatmap(X, Y, o)
    ax.imshow(img, extent=extent, origin='lower', cmap=cm.jet)
    

def hvobj(dim):
    return get_performance_indicator("hv", ref_point=np.array([1.001] * dim))

def hvobj_custom(dim, p):
    return get_performance_indicator("hv", ref_point=p)

def calchv(hvobj, points):
    return hvobj.do(points)

def printhv(hvobj, points):
    print("HV:", round(calchv(hvobj, points) * 100, 3),"%")
    
def plot_points3D(ax, points, color='black'):
    ax.plot(points[:,0], points[:,1], points[:,2],'o',color=color)
    
'''Sphere Methods'''
def randomSpherePoint():
    X = abs(np.random.normal(size=3))
    normalize = (np.sum(X**2))**(1/2)
    return X / normalize

def eluc_dis(a, b):
    return np.sum((a-b)**2)**(1/2)

def circle_dis(a,b,r=1):
    ratio = (eluc_dis(a,b)/2) / r
    angleratio = np.degrees(np.arcsin(ratio)) / 360
    length = (2 * np .pi * r) * angleratio
    return length

def maximinSphere(pts, its, distance, randomizer):
    return maximin(pts,its, distance,randomizer)**2

def maximinAvgSphere(pts, its, distance, randomizer):
    return maximin_avg(pts,its, distance,randomizer)**2

def SeqmaximinSphere(pts, its, distance, randomizer, past=[]):
    return SeqMaximin(pts,its, distance,randomizer, past)**2

'''Space-Filling Methods'''
#maxiavg
def all_dis(pt, points, index, dfunc):
    distances = []
    for i in range(len(points)):
        if(i != index):
            distances.append(dfunc(pt, points[i]))
    return np.array(distances)

def getAllDis(pts):
    results = np.zeros(shape=(len(pts)))
    for i in range(len(pts)):
        results[i] = min_distance(pts[i], pts, i, eluc_dis)
    return results

def min_distance(pt, points, index, dfunc):
    mind = 10000
    for i in range(len(points)):
        if(i != index):
            dis = dfunc(pt, points[i])
            mind = min(mind, dis)
    return mind

def avg_distance(pt, points, index, dfunc):
    mind = 0
    if(len(points) == 0): return 0
    for i in range(len(points)):
        if(i != index):
            dis = dfunc(pt, points[i])
            mind += dis
    return (mind / len(points))

def maximin(pts, its, distance, randomizer):
    D = 3
    N = pts
    points = np.array(np.arange(D * N).reshape(N, D), dtype='float')
    for i in range(len(points)):
        points[i] = randomizer()
    current_mins = [-1.0] * N
    for _ in range(its):
        index = int(np.random.rand() * N)
        if(current_mins[index] < 0):
            current_mins[index] = min_distance(points[index], points, index, distance)

        canidate = randomizer()
        canidate_min = min_distance(canidate, points, index, distance)

        if(canidate_min > current_mins[index]):
            points[index] = canidate
            current_mins[index] = canidate_min
    return points;

def maximin_avg(pts, its, distance, randomizer):
    D = 3
    N = pts
    points = np.array(np.arange(D * N).reshape(N, D), dtype='float')
    for i in range(len(points)):
        points[i] = randomizer()
    current_mins = [-1.0] * N
    for _ in range(its):
        index = int(np.random.rand() * N)
        if(current_mins[index] < 0):
            current_mins[index] = avg_distance(points[index], points, index, distance)

        canidate = randomizer()
        canidate_min = avg_distance(canidate, points, index, distance)

        if(canidate_min > current_mins[index]):
            points[index] = canidate
            current_mins[index] = canidate_min
    return points;

def SeqMaximin(pts, its, distance, randomizer, past=[]):
    D = 3
    N = pts
    points = np.array(np.arange(D * N).reshape(N, D), dtype='float')
    for i in range(len(points)):
        points[i] = randomizer()
    current_mins = [-1.0] * N
    canidate_min = -1
    for _ in range(its):
        index = int(np.random.rand() * N)
        if(current_mins[index] < 0):
            current_mins[index] = min(min_distance(points[index], points, index, distance),
                                      min_distance(points[index], past, -1, distance))

        canidate = randomizer()
        canidate_min = min(min_distance(canidate, points, index, distance),
                           min_distance(canidate, past, -1, distance))

        if(canidate_min > current_mins[index]):
            points[index] = canidate
            current_mins[index] = canidate_min
    return points;

'''Simplex Methods'''
def random2SimplexPoint():
    rnums = np.random.rand(2)
    rnums = np.append(rnums, np.array([0,1]))
    rnums.sort()
    diffs = np.array([0.0] * 3)
    for i in range(len(rnums)-1):
        diffs[i] = (rnums[i+1] - rnums[i])
    return diffs

def randomSimpexPoint(dimensions): #dimensions=2, is a 2-simplex=triangle
    rnums = np.random.rand(dimensions)
    rnums = np.append(rnums, np.array([0,1]))
    rnums.sort()
    diffs = np.array([0.0] * (dimensions + 1))
    for i in range(len(rnums)-1):
        diffs[i] = (rnums[i+1] - rnums[i])
    return diffs

def randomSimplexPrismPoint(dimensions): #dimensions-1=simplex dimnsions
    return np.append(randomSimpexPoint(dimensions-1), np.random.rand())
        
def map_hypercuble_to_simplex(points):
    rpoints = np.array(np.arange(len(points) * (len(points[0]) + 1)).reshape(len(points), len(points[0]) + 1),dtype=float)
    for i in range(len(points)):
        pts = np.copy(points[i])
        pts = np.append(pts, np.array([0,1]))
        pts = np.sort(pts)
        for j in range(1, len(pts)):
            rpoints[i][j-1] = pts[j] - pts[j-1]
        # print(pts, rpts)
    return rpoints
                
def dich(num):
    points = np.array([[0.0,0.0,0.0]] * (num))
    for i in range(num):
        p = np.random.dirichlet((1,1,1))
        points[i] = np.array(p)
    return points;

def getRD(M, N, seed=128):
    return get_reference_directions("energy", M, N, seed=seed)

'''Hypercube Methods'''
def unif(num):
    points = np.array([[0.0,0.0,0.0]] * (num))
    for i in range(num):
        p = np.random.uniform(size=num)
        points[i] = np.array(p)
    return points;


'''performace metric methods'''

class performance_metrics:
    type = ""
    radius = 0
    def __init__(self):
        pass
    def get_indicator(name, radius=0.0):
        obj = performance_metrics()
        obj.type = name
        obj.radius = radius
        return obj
        
    def do(self, points):
        if(self.type == "rad-simplex"):
            if(len(points) == 0):  return 0  
            n = len(points[0])-1
            area = ((n+1)**(1/2))/(math.factorial(n))
            return self.do_radial_performance_monto_carlo(points, self.radius, randomSimpexPoint, area)
        if(self.type == "rad-simplexprism"):
            if(len(points) == 0):  return 0  
            n = len(points[0])-2
            area = ((n+1)**(1/2))/(math.factorial(n)) * 1.0
            return self.do_radial_performance_monto_carlo(points, self.radius, randomSimplexPrismPoint, area)
        if(self.type == "dmin"):
            return self.do_dmin_performance(points)
        if(self.type == "vgm"):
            return self.do_vgm_performance(points)
    
    def do_radial_performance_monto_carlo(self, points, radius, point_generator, shapearea):
        constants = [0,0,1,4/3,1/2,8/15,1/6]
        accuracy_iterations = 10000
        total_points_in_circles, points_in_at_least_one_circle, point_overlapping_circles = 0,0,0
        
        dimensions = len(points[0])
        for _ in range(accuracy_iterations):
            point = point_generator(dimensions-1)
            count = 0
            for p in points:
                hypot = eluc_dis(p, point)
                if(hypot <= radius):
                    count += 1
            if(count > 0):
                point_overlapping_circles += count - 1
                total_points_in_circles += count
                points_in_at_least_one_circle += 1
        
        volume_ratio = points_in_at_least_one_circle / accuracy_iterations # ratio of covered volume / total volume 
        volume_overlap_ratio = point_overlapping_circles / total_points_in_circles # ratio: # points overlap / # points

        total_volume = len(points) * (radius)**(dimensions-1) * np.pi * constants[dimensions-1]
        total_volume_in_simplex = (total_points_in_circles / (accuracy_iterations)) * shapearea

        volume_overflow = total_volume - total_volume_in_simplex
        volume_overflow_ratio = volume_overflow / total_volume # ratio of overflow volume / total volume by spheres

        return [volume_ratio,volume_overlap_ratio,volume_overflow_ratio]
        # extra = (p[-1] - point[-1])
        # dis = (hypot**2 - extra**2)**(1/2)
        # total_volume_covered = (3**(1/2))/2 * volume_ratio
        # total_volume_overlap = (3 ** (1/2))/2 * (point_overlapping_circles / (accuracy_iterations * len(points)))
        
        
    def do_dmin_performance(self, points):
        mdis = 1000
        for i in range(len(points)):
            mdis = min(mdis, min_distance(points[i], points, i, eluc_dis))
        return mdis
   
    def getVar(self, points):
        mean = np.sum(points) / len(points)
        top = np.sum((points-mean)**2)
        bottom = len(points)-1
        SD = (top/bottom)
        return SD

    def getGeoMean(self, points):
        geomean = prod(points)**(1/3)

        return geomean
  
    def getGeoVar(self, points):
        geoVar = np.exp((len(points)-1) / len(points) * self.getVar(np.log(points)))
        # geoVar2 = np.exp(np.log((self.getGeoSD(points) ** 2))) #untested
        return geoVar
    
    def getGeoSD(self, points):
        geoSD = np.exp(np.sqrt((len(points)-1) / len(points) * self.getVar(np.log(points))))
        # geoSD2 = np.exp((sum(np.log(points / self.getGeoMean(points))**2) / 3)**(1/2))
        return geoSD
    
    def do_vgm_performance(self, points):
        means = []
        for i in range(len(points)):
            k = len(points[0])-1
            distances = all_dis(points[i], points, i, eluc_dis)
            indx = np.argpartition(distances, k-1)
            part_pts = points[indx[:k]]
            means.append(self.getGeoMean(part_pts))
        return self.getVar(np.array(means))
            
    def do_radial_performance_partial(circles, radius):
        pass
   
    def do_radial_performance_systematic(circles, radius):
        pass
   
    def do_radial_performance_integration(circles, radius):
        pass




''' DP Methods '''


def layer_wise_maximin(points, low, high, otherpoints=np.empty((0,3))):
    points = np.append(points, np.random.uniform(low=low, high=high, size = (len(points), 1)), axis=1)
    iterations = 10000
    
    
    for i in range(len(points)):
        points[i][-1] = np.random.uniform(low=low, high=high)
    
    for _ in range(iterations):
        i = int(np.random.rand() * len(points))
        dis = min(min_distance(points[i], points, i, eluc_dis),
                  min_distance(points[i], otherpoints, -1, eluc_dis))
        save = points[i][-1];
        points[i][-1] =  np.random.uniform(low=low, high=high)
        dis2 = min(min_distance(points[i], points, i, eluc_dis), 
                   min_distance(points[i], otherpoints, -1, eluc_dis))
        
        if(dis > dis2):
            points[i][-1] = save
            
    return points

def multilayer_projection_DP_Maximin_Energy(n_per_layer, n_bins):
    bin_starts = np.linspace(0, 1, n_bins+1)
    bin_sz = bin_starts[1] - bin_starts[0]
    
    retpoints = np.empty((0, 3))
    seqpoints = np.empty((0,3))
    
    for bin_st in bin_starts[:-1]:
        points = getRD(2, n_per_layer, seed=np.random.randint(low=10, high=1000))
        points_full = layer_wise_maximin(points, bin_st, bin_st + bin_sz, otherpoints=seqpoints)
        retpoints = np.append(retpoints, points_full, axis=0)
        # seqpoints = np.append(seqpoints, (bin_st + bin_sz) - (points_full - bin_st), axis=0)
        seqpoints = np.append(seqpoints, points_full, axis=0)
        seqpoints[:,2] += bin_sz
        #invert prev row 
        #SAW ROW, DIFF POINTS
        
    return retpoints
        
def segmented_simplex_projection_DP_Maximin_Energy(n_total, n_bins):
    bin_starts_x = np.linspace(0, 1, n_bins+1)
    bin_starts_y = np.linspace(1, 0, n_bins+1)
    diff = bin_starts_x[1] - bin_starts_x[0]
    points = getRD(2, n_total, seed=np.random.randint(low=10, high=1000))
    
    # sorted_points = points[np.argsort(points[:, 0])]
    final_points = np.empty((0,3))
    for bin_st in bin_starts_x[:-1]:
        points_in = points[(points[:,0] >= bin_st) & (points[:,0] <= bin_st+diff)]
        # points[(points[:,0] > bin_starts_x) & (points[:,0] < (bin_starts_x+diff)) & (points[:,2] > 10)]
        new_points = layer_wise_maximin(points_in, 0, 1)
        final_points = np.append(final_points, new_points, axis=0)
    return final_points
            
def duplicated_projected_DP_Energy(n_per_layer, n_bins):
    bin_starts = np.linspace(0, 1, n_bins+1)
    bin_sz = bin_starts[1] - bin_starts[0]
    retpoints = np.empty((0, 3))
     
    base_points = getRD(2, n_per_layer, seed=np.random.randint(low=10, high=1000))
    
    basepts = np.append(base_points, np.zeros((len(base_points), 1), dtype=float), axis=1)
    
    for bin_st in bin_starts[:-1]:
            retpoints = np.append(retpoints, basepts, axis=0)
            retpoints[:,2] += bin_sz
            
    return retpoints   


#make squaree go all thee way arround
def no_segmentation_projection_DP_Maximin_Energy(n, sq_pts=10,dis=0.03):
    total = np.linspace(0, 2**(1/2) * 2 + 2, sq_pts+1)[:-1]
    
    sq = 2**(1/2)
    sin2 = sq / 2
    
    
    X1 = total[total <= sq]
    Y1 = total[(total > sq) & (total <= sq+1)] - sq;
    X2 = sq - ((total[(total > sq+1) & (total <= 2*sq+1)]) - (sq+1));
    Y2 = total[total > 2*sq+1] - (2*sq+1)
    
    
    pts = np.array([])
    
    a1 = np.reshape(X1,( len(X1),1)) * sin2
    a2 = np.reshape(X2, (len(X2), 1)) * sin2
    
    
    r1 = np.append(np.append(a1, 1-a1, axis=1), np.full((len(a1),1),-dis), axis=1)
    r2 = np.append(np.append(a2, 1-a2, axis=1),np.full((len(a2),1),1+dis), axis=1)
    
    r3 = np.concatenate((np.full((len(Y1),1),1+dis), 
                         np.full((len(Y1),1), -dis),
                         np.reshape(Y1, (len(Y1),1))
                         ),axis=1)
    
    
    r4 = np.concatenate((np.full((len(Y2),1),-dis),
                         np.full((len(Y2),1),1+dis),
                         1 - np.reshape(Y2, (len(Y2),1))
                         ),axis=1)


    pts = np.concatenate((r1,r3,r2,r4));
    
    points = getRD(2, n, seed=np.random.randint(low=10, high=1000))
    new_points = layer_wise_maximin(points,0,1,otherpoints=pts)
    return new_points;
    
def getMetricString(pts, vgm=None,rad=None, dmin=None):
    res = ''
    if vgm is not None:
        res += 'VGM(< Better) -> ' + str(round(vgm.do(pts) * 10 ,2)) + '\n';
    if rad is not None:
        rs = rad.do(pts)
        res += 'Vol(> Better) -> ' + str(round(rs[0] ,2)) + '\n';
        res += 'Ovlap(< Better) -> ' + str(round(rs[1] ,2)) + '\n';
        res += 'Ovflw(< Better) -> ' + str(round(rs[2] ,2)) + '\n';
    if dmin is not None:
        res += 'Dmin(> Better) -> ' + str(round(dmin.do(pts),2)) + '\n';
    
    return res;


# dis = sqrt((1-x)^2+(y)^2)


# dis = m.getAllDis(pts)
# ax3.tricontourf(V1, pts[:,2], dis , levels=np.linspace(0, .3, 100))
# ax3.set(xlim=(0, 2**(1/2)), ylim=(0, 1))

#marginals
#fix holE(√)
#the paper about optimizations
#remove bin(√))