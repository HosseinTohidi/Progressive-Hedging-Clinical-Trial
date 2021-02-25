### grouping

from copy import deepcopy
import numpy as np

def myClusters(X, k, oper = 'min'):
    """ 
    Returns a dictionary containing k (or k+1) clusters
    in case oper = 'max' it will return k+1 clusters and the k+1th is the dissimilar group
    and makes the clusters a cover not a partition. The rest oper s will return a partition of clusters.
    
    Parameters
    ----------
    X : 2D list/numpy 2D array
        a list or array to be clustered!
    k : integer num
        the number of clusters
        
    oper: str
        can be either 'min', 'max', 'random'
    """
    def dist(a, b, ax=1):
        return np.linalg.norm(a - b, axis=ax)
    # X coordinates of random centroids
    def find_cores(X, k):
        new_inds = []
        new_X = deepcopy(X)
        cores = []
        idx = np.random.randint(0,len(X))
        cores.append(new_X[idx])
        new_X = np.delete(new_X, idx, axis = 0)  #remove the rows if already selected as a core       
        W_distance = dist(new_X, cores[0])
        W_distance = W_distance / W_distance.sum() 
        chance = W_distance.cumsum()
        for kk in range(k-1):
            rnd = np.random.random()
            idx = len(chance[chance<rnd])-1
            cores.append(new_X[idx])
            new_X = np.delete(new_X, idx, axis = 0)  #remove the rows if already selected as a core
        return np.array(cores)       
    #C = np.random.randint(0, np.max(X)//2, size=[k,len(X[0])])
    if oper != 'random':   # min or max
        C = find_cores(X, k)
        # To store the value of centroids when it updates
        C_old = np.zeros(C.shape)
        # Cluster Lables(0, 1, 2)
        clusters = np.zeros(len(X))
        # Error func. - Distance between new centroids and old centroids
        error = dist(C, C_old, None)
        counter = 0
        # Loop will run till the error becomes zero
        while error >= 0.000001:
            counter +=1
            # Assigning each value to its closest cluster
            for i in range(len(X)):
                distances = dist(X[i], C)
                cluster = np.argmin(distances)                
                clusters[i] = cluster
            # Storing the old centroid values
            C_old = deepcopy(C)
            # Finding the new centroids by taking the average value
            for i in range(k):
                points = [X[j] for j in range(len(X)) if clusters[j] == i]
                C[i] = np.mean(points, axis=0)
            error = sum(dist(C, C_old))
        cDict = {}
        for i in range(len(X)):
            if cDict.get(clusters[i],[]) == []:
                cDict[clusters[i]] = [i]
            else:
                cDict[clusters[i]].append(i)               
        if oper == 'max': #create dissimilarity group
            l = len(cDict.keys())
            cDict[l]= []
            for i in range(k):
                points = np.array([X[j] for j in range(len(X)) if clusters[j] == i])
                mean_point = points.mean(axis=0)
                idx = np.argmin(dist(mean_point,points))
                b = cDict[i][idx]
                cDict[l].append(b)   
                cDict[i].pop(idx) #comment this line if you want a cover! not a partition
    else:
        cDict = {}
        rnd = np.argsort(np.random.random(len(X)))
        length = len(X)// k
        for i in range(k-1):
            cDict[i] = [rnd[j] for j in range(i*length,(i+1)*length)]
        cDict[k-1] = [rnd[j] for j in range((k-1)*length,len(X))]    
    return cDict

        




    
    
    
    