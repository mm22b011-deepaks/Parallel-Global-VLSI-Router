import numpy as np
import sys

def read_out_file(filename):
    wiring = []
    with open(filename, 'r') as file:
        # Skip lines starting with ';' and empty lines until the first non-empty, non-comment line
        line = file.readline().strip()
        while line.startswith(';') or not line:
            line = file.readline().strip()
        
        # Extract M and N
        M, N = map(int, line.split())
        
        # Read data for array Gridx
        Gridx = np.zeros((M, N+1), dtype=int)
        for i in range(M):
            line = file.readline().strip()
            Gridx[i] = list(map(int, line.split()))
        
        # Read data for Gridy
        Gridy = np.zeros((M+1,N), dtype=int)
        for i in range(N):
            line = file.readline().strip()
            Gridy[i] = list(map(int, line.split()))
        
        # Read data for the list of arrays
        while True:
            line = file.readline().strip()
            if not line:
                break
            wiring.append(np.array(list(map(int, line.split()))))
    
    return Gridx, Gridy, wiring


def read_in_file(filename):
    nets = []
    with open(filename, 'r') as file:
        # Skip lines starting with ';' and blank lines
        lines = [line.strip() for line in file if line.strip() and not line.strip().startswith(';')]
        
        # Extract M, N, C, and v
        M, N, C, v = map(int, lines[0].split())
        
        # Read data into a NumPy array
        for line in lines[1:]:
            nets.append(np.array(list(map(int, line.split()))))
    
    return M, N, C, v, np.array(nets)

import numpy as np

def overflow(Gridx, Gridy, C):
    # Calculate the overflow based on the difference between the value and the threshold C
    overflow = np.sum(np.maximum(Gridx - C, 0)) + np.sum(np.maximum(Gridy - C, 0))
    return overflow

def cost(paths):
    cost = 0
    for i in paths:
        n = len(i)//2
        # wire length cost
        for j in range(n-1):
            if (i[2*j]==i[2*j+2]):
                cost += np.abs(i[2*j+1]-i[2*j+3])
            elif (i[2*j+1]==i[2*j+3]):
                cost += np.abs(i[2*j]-i[2*j+2])
            else:
                print("\nBad wiring; ",i)
                return -1
        cost += v*(n-2) # via cost
    return cost

def Chkgridgraph(Gridx,Gridy,paths):
    X = np.zeros_like(Gridx)
    Y = np.zeros_like(Gridy)
    for i in paths:
        n = len(i)//2
        for j in range(n-1):
            if (i[2*j]==i[2*j+2]):
                cost += np.abs(i[2*j+1]-i[2*j+3])
            elif (i[2*j+1]==i[2*j+3]):
                cost += np.abs(i[2*j]-i[2*j+2])
            else:
                print("\nBad wiring; ",i)
                return -1
        cost += v*(n-2) # via cost


# Usage
filename = sys.argv[1]
M, N, C, v, data_array = read_in_file(filename)

print("M:", M)
print("N:", N)
print("C:", C)
print("v:", v)
print("Data array:")
print(data_array)

filename = sys.argv[2]  
Gridx, Gridy, wiring = read_out_file(filename)

print("Gridx:")
print(Gridx)
print("\nGridy:")
print(Gridy)
print("\nList of wiring paths:")
for arr in wiring:
    print(arr)

print("\nShorts = ",overflow(Gridx, Gridy, C))
print("\ncost = ",cost(wiring))