#! /usr/bin/python3

import argparse
from mpl_toolkits import mplot3d 
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

parser = argparse.ArgumentParser(description='Plots the spherical plot of orientaitons.')
parser.add_argument('filename', metavar='f', type=str, help='input file name')
args = parser.parse_args()
  
# Creating dataset 
#z = 4 * np.tan(np.random.randint(10, size =(500))) + np.random.randint(100, size =(500)) 
#x = 4 * np.cos(z) + np.random.normal(size = 500) 
#y = 4 * np.sin(z) + 4 * np.random.normal(size = 500) 
#c = (x + y + z)

x = []
y = []
z = []
c = []

fi = open(args.filename)
for line in fi:
  cacca = line.split()
  if cacca == []:
    continue
  x.append( float ( cacca[0] ) )
  y.append( float ( cacca[1] ) )
  z.append( float ( cacca[2] ) )
  c.append( float ( cacca[4] ) )
  
x = np.array(x)
y = np.array(y)
z = np.array(z)
c = np.array(c)
  
# Creating figure 
fig = plt.figure(figsize = (16, 9)) 
ax = plt.axes(projection ="3d") 
    
# Add x, y gridlines  
ax.grid(b = True, color ='grey',  
        linestyle ='-.', linewidth = 0.3,  
        alpha = 0.2)  
  
  
# Creating color map
#my_cmap = ListedColormap(["white", "white", "white", "white", "red", "black"])
#my_cmap = plt.get_cmap('BuPu') 
colors=["white", "white", "white", "white", "white", "white", "orange", "red", "black"]
my_cmap = LinearSegmentedColormap.from_list("cacca", colors)
  
# Creating plot 
sctt = ax.scatter3D(x, y, z, 
                    alpha = 0.8, 
                    c = c,  
                    cmap = my_cmap,  
                    marker ='o') # . o 8 
  
plt.title("simple 3D scatter plot") 
ax.set_xlabel('X-axis', fontweight ='bold')  
ax.set_ylabel('Y-axis', fontweight ='bold')  
ax.set_zlabel('Z-axis', fontweight ='bold') 
fig.colorbar(sctt, ax = ax, shrink = 0.5, aspect = 5) 


plt.axis('off')
plt.grid(b=None)
  
# show plot 
plt.show() 

