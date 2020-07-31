#! /usr/bin/python3

import argparse
from mpl_toolkits import mplot3d 
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.colors import LinearSegmentedColormap

parser = argparse.ArgumentParser(description='Plots the spherical plot of orientaitons.')
parser.add_argument('filename', metavar='f', type=str, help='input file name')
args = parser.parse_args()

azimu = [] # 0 -> 2pi longitude
polar = [] # 0 -> pi  latitude
color = []

fi = open(args.filename)
#isFirst = True
for line in fi:
  cacca = line.split()
  azimu.append( float ( cacca[0] ) )
  polar.append( float ( cacca[1] ) )
  color.append( float ( cacca[2] ) )

x = np.zeros( (80, 40) )
y = np.zeros( (80, 40) )
z = np.zeros( (80, 40) )
c = np.zeros( (80, 40) )

minC = 1
maxC = 0
for i in range(80):
  for j in range(40):
    x[i,j] = np.cos(azimu[i*40])*np.sin(polar[j])
    y[i,j] = np.sin(azimu[i*40])*np.sin(polar[j])
    z[i,j] = np.cos(polar[j])
    aa     = color[i*40+j]
    c[i,j] = aa
    if aa < minC:
      minC = aa
    if aa > maxC:
      maxC = aa
cRange = maxC - minC

# Creating color map
#colors=["white", "white", "white", "white", "white", "white", "orange", "red", "black"]
colors=["black", "blue", "cyan", "white", "yellow", "orange", "red"]
my_cmap = LinearSegmentedColormap.from_list("cacca", colors)

print(minC, maxC, cRange)

my_facecolors = my_cmap(  (c - minC)/cRange  )

# Creating figure 
fig = plt.figure()
ax = plt.axes(projection ="3d") 
    
# Add x, y gridlines  
ax.grid(b = True, color ='grey',  
        linestyle ='-.', linewidth = 0.3,  
        alpha = 0.2)


#plot = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=my_cmap, linewidth=0, antialiased=False, alpha=0.5)
plot = ax.plot_surface(x, y, z, rstride=1, cstride=1, facecolors=my_facecolors, linewidth=0, antialiased=False, alpha=0.5)

plt.axis('off')
plt.grid(b=None)

plt.show()
