from random import uniform
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')
import numpy as np
from mpl_toolkits import mplot


def random(D, N):
    points = []
    for i in range(D):
        points.append([])
    for i in range(N):
        for j in range(D):
            points[i][j].append(rnum())
    
    if(D == 2):
        plt.plot(points[0], points[1], 'o', color='black')
        
        
    

def main(): 
    D = input("Dimensions: ")
    N = input("Points:")
    random(D, N);
    
    
    
            
def rnum():
    return uniform(0, 1);
            
if __name__ == "__main__":
    main()