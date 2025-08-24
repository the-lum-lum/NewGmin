#This script will read in the coords files to visuliase the elastic filaments
import numpy as np
import matplotlib.pyplot as plt


#1) Choose the file we want to visualise
path = r"C:\Users\callu\Documents\NewGmin\user\lamina_phase_diagram"
filename=r"\lamina_NROD10_NSEG_200_RAT5\Iteration_4_plot"
mycoords=(path + filename + "\\lowest")
nheader=4
myinits=(path + filename + "\\initpos")
nheaderinit=0

#2) User inputs
N_ROD=10
N_SEG=200
L=5

#3) Read in the coordinates
fid=open(mycoords,'r')
angles=np.loadtxt(fid,skiprows=nheader)
fid.close
print(len(angles))

fid=open(myinits,'r')
inits=np.loadtxt(fid,skiprows=nheaderinit)
fid.close

#4) Restructure the data and convert to xy
LS=L/N_SEG
coords=np.zeros([N_ROD,N_SEG,2])
coords[:,0,0]=inits[:,0]
coords[:,0,1]=inits[:,1]
s=0
for i in range(0,N_ROD):
	for j in range(1,N_SEG):
		coords[i,j,0]=coords[i,j-1,0]+LS*np.cos(angles[s])
		coords[i,j,1]=coords[i,j-1,1]+LS*np.sin(angles[s])
		s=s+1

#5) View the results
fig1=plt.figure()
ax1 = plt.axes()

for i in range(0,N_ROD):
	plt.plot(coords[i,:,0],coords[i,:,1])

plt.show()

