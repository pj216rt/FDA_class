# https://fdasrsf-python.readthedocs.io/en/latest/time_warping_example.html
	# Documentation / User Guide
# https://github.com/jdtuck/fdasrsf_python/tree/master
	# Github front page
# https://github.com/jdtuck/fdasrsf_python/blob/master/notebooks/example_warping.ipynb
	# Examples
# https://github.com/jdtuck/fdasrsf_python/blob/master/fdasrsf/time_warping.py
	# Used circa lines 75-160 to make pairwise warp example

# Install was simple for me with
	# pip install fdasrsf

# I couldn't find the simu_data.npz anywhere along with the install, 
# but was able to use it after downloading it separately at
# https://github.com/jdtuck/fdasrsf_python/blob/master/bin/simu_data.npz
# and load it manually from a directory.

# I wasn't able to get a non-srvf version out of Dr. Tucker's code, but I translated a (poor-quality) Matlab version that I wrote myself for a class. 
# It runs for me as written here, but I haven't tested it beyond this. I do know that changing the maxStep value can break it, but the default of 4 worked for me.

# To run the parts of this file:
# First run the function at the bottom of the file. (Only needed for the Function F warping).
# Put a '#' in front of the ''' at the beginning of a section/subsection to uncomment a section/subsection block.
# The ExampleFromGithub section does common registration of multiple functions to a mean, as well as some other things we haven't covered in class yet.
# The Function F warping subsection runs a (bad) DP registration algorithm to directly register one function to another, producing a pinching effect.
# The SRVF warping subsection uses a piece of Dr. Tucker's code to properly register one function directly to another using their srvfs.

##############################################################################
# Contents 
##############################################################################

## << ExampleFromGithub>
## << PairwiseWarping>
##		< Function F warping
##		< SRVF Q warping

##############################################################################
##############################################################################

##############################################################################
''' %% ExampleFromGithub> 
##############################################################################
import os
import fdasrsf as fs
import numpy as np
import matplotlib.pyplot as plt

# Load
os.chdir(yourDataFolder)
data = np.load('./simu_data.npz')
#
time = data['arr_1']
f = data['arr_0']

# Alignment
obj = fs.fdawarp(f,time) # Create fdawarp object
obj.srsf_align(parallel=True) # Align functions
obj.plot()

# Vertical PCA - analyze amplitude
vpca = fs.fdavpca(obj)
vpca.calc_fpca(no=10)
vpca.plot()

# Horizontal PCA - analyze phase
hpca = fs.fdahpca(obj)
hpca.calc_fpca(no=3)
hpca.plot()

# Joint fPCA - analyze amplitude and phase jointly
jpca = fs.fdajpca(obj)
jpca.calc_fpca(no=10)
jpca.plot()

##############################################################################
### >>'''
##############################################################################

##############################################################################
''' #%% PairwiseWarping> 
##############################################################################
import os
import numpy as np
import fdasrsf as fs
import fdasrsf.utility_functions as uf
import matplotlib.pyplot as plt

# Load
os.chdir(yourDataFolder)
data = np.load('./simu_data.npz')
#
fTime = data['arr_1']
f = data['arr_0']

fSet = [3,18]
# # fSet = np.arange(0,20)
# fig,ax = plt.subplots()
# for fun in fSet:
# 	ax.plot(fTime,f[:,fun], label=str(fun))
# plt.legend()
# # plt.show()
# savePlotToFile(fileName1,svDir,True)

fSmall = f[:,fSet]
M = fSmall.shape[0]
N = fSmall.shape[1]


#""" # # # # Function F warping # # # # # # # # # # # # # # # # # # # # # # # # $ 

fWarpedF = np.zeros([M,N])
gamF = np.zeros([M,N])
fWarpedF[:,0],gamF[:,0] = badDP_F(fSmall[:,1],fSmall[:,0],fTime)
fWarpedF[:,1],gamF[:,1] = badDP_F(fSmall[:,0],fSmall[:,1],fTime)

# The legends here aren't working for me. Let me know if you know why.
fig,(ax1,ax2,ax3) = plt.subplots(3,1, sharey=True)
for k in range(0,N):
	ax1.plot(fTime,fSmall[:,k], label=str(k+1))
plt.legend()
ax2.plot(fTime,fWarpedF[:,0], label='1 warped')
ax2.plot(fTime,fSmall[:,1], label='2')
plt.legend()
ax3.plot(fTime,fSmall[:,0], label='1')
ax3.plot(fTime,fWarpedF[:,1], label='2 warped')
plt.legend()
plt.show()
#
fig,ax4 = plt.subplots(1,1, sharey=True)
ax4.plot(fTime,gamF[:,0], label='gamma 1')
ax4.plot(fTime,gamF[:,1], label='gamma 2')
plt.legend()
plt.show()


### """ # # # # # # # # # # # # # # # # # # # # # # # # 

#"""  # # # # SRVF Q warping # # # # # # # # # # # # # # # # # # # # # # # # # # #

omethod = "DP2"
grid_dim = 7
smoothdata = False
lam = 0.0 
eps = np.finfo(np.double).eps

fSmallQ,g,_ = uf.gradient_spline(fTime,fSmall,smoothdata)
q = g/np.sqrt(abs(g) + eps)
gamQ = np.zeros([M,N])
gamQ[:,0] = uf.optimum_reparam(q[:,1],fTime,q[:,0], omethod,lam,grid_dim) # Warp first function to match second.
gamQ[:,1] = uf.optimum_reparam(q[:,0],fTime,q[:,1], omethod,lam,grid_dim) # Warp second function to match first.

fSmallWarpedQ = np.zeros([M,N])
for k in range(0,N):
	fSmallWarpedQ[:,k] = np.interp( (fTime[-1]-fTime[0]) * gamQ[:,k] + fTime[0], fTime,fSmallQ[:,k])

# The legends here aren't working for me. Let me know if you know why.
fig,(ax1,ax2,ax3) = plt.subplots(3,1, sharey=True)
for k in range(0,N):
	ax1.plot(fTime,fSmallQ[:,k], label=str(k+1))
plt.legend()
ax2.plot(fTime,fSmallWarpedQ[:,0], label='1 warped')
ax2.plot(fTime,fSmallQ[:,1], label='2')
plt.legend()
ax3.plot(fTime,fSmallQ[:,0], label='1')
ax3.plot(fTime,fSmallWarpedQ[:,1], label='2 warped')
plt.legend()
plt.show()

### """ # # # # # # # # # # # # # # # # # # # # # # # # 

################################################################################
### >>'''
################################################################################

################################################################################
#'''<< Makeshift time-register DP
################################################################################

def badDP_F(f1,f2,t, maxStep=4):
	'''Creates a warping function gamma to f2 to register to f1.
		Returns f2Warped,gamma.
		Hastily translated from amateur Matlab code. NOT an example of good python code.'''
	import math
	import numpy as np

	grain = t.shape[0]
	S = np.inf*np.ones([grain,grain])
	S[-1,-1] = 0
	# Allowed node-to-node steps in (x,y) position changes.
	# Contains redundancies that will slow things down.
	stepRn = np.vstack(np.arange(maxStep)+1)
	stepOnes = np.ones([maxStep,1])
	steps = np.hstack([np.kron(stepOnes,stepRn), np.kron(stepRn,stepOnes)])

	stepGrain = 1 # Will be lcm of step size integers.
	for i in stepRn:
		stepGrain = stepGrain*i/math.gcd(int(stepGrain),int(i))
	stepGrain=stepGrain[0]; stepGrain=int(stepGrain)

	# Determine the optimal node values.
	crumbs = np.inf*np.ones([grain-1,grain-1]) # Leave some crumbs behind to mark the trail.
	stepsHeight = steps.shape[0]
	h = np.zeros(stepsHeight)
	# for I in reversed(np.arange(gridSize)+1):
	for I in reversed(np.arange(grain-1)):
		J = I
		j = J; j=int(j)
		# Traverse the highest unevaluated row horizontally.
		for i in reversed(np.arange(I+1)):
			i=int(i)
			# Find potential values.
			for s in np.arange(stepsHeight):
				h[s] = np.inf # Deals with non-accessible nodes.
				k = i+steps[s,0]; k=int(k)
				l = j+steps[s,1]; l=int(l)
				if all([k<=grain-1, l<=grain-1]):
					# Use interpolation to match the intervals for f1(t) and f2(t).
					numSteps = stepGrain*(k-i)
					x1 = t[i:k]
					v1 = f1[i:k]
					xf1 = np.concatenate([np.arange(t[i],t[k], (t[k]-t[i])/numSteps), np.array([t[k]]) ]) # np.arange stops short of the final value.
					vf1 = np.interp(xf1,x1,v1)
					x2 = t[j:l]
					v2 = f2[j:l]
					xf2 = np.concatenate([np.arange(t[j],t[l], (t[l]-t[j])/numSteps), np.array([t[l]]) ]) # np.arange stops short of the final value.
					vf2 = np.interp(xf2,x2,v2)
					# The potential node value is the new link cost plus the optimum of the destination node.
					linkCost = np.trapz((vf1-vf2)**2, xf1)
							# For proper SRVF version, use 
							# linkCost = np.trapz((vf1-vf2*np.sqrt((l-j)/(k-i)))**2, xf1)
							# (I haven't tested this in the python version)
					h[s] = linkCost+S[k,l]
			# Select the optimal value and record its trail.
			S[i,j] = np.min(h)
			crumbs[i,j] = np.argmin(h)
		#
		i = I
		# Traverse the rightmost unevaluated column vertically.
		# This is exactly the same as the above loop.
		for j in reversed(np.arange(J)):
			j=int(j)
			for s in np.arange(stepsHeight):
				h[s] = np.inf # Deals with non-accessible nodes.
				k = i+steps[s,0]; k=int(k)
				l = j+steps[s,1]; l=int(l)
				if all([k<=grain-1, l<=grain-1]):
					# Use interpolation to match the intervals for f1(t) and f2(t).
					numSteps = stepGrain*(k-i)
					x1 = t[i:k]
					v1 = f1[i:k]
					xf1 = np.concatenate([np.arange(t[i],t[k], (t[k]-t[i])/numSteps), np.array([t[k]])]) # np.arange stops short of the final value.
					vf1 = np.interp(xf1,x1,v1)
					x2 = t[j:l]
					v2 = f2[j:l]
					xf2 = np.concatenate([np.arange(t[j],t[l], (t[l]-t[j])/numSteps), np.array([t[l]])]) # np.arange stops short of the final value.
					vf2 = np.interp(xf2,x2,v2)
					# The potential node value is the new link cost plus the optimum of the destination node.
					linkCost = np.trapz((vf1-vf2)**2, xf1)
							# For proper SRVF version, use 
							# linkCost = np.trapz((vf1-vf2*np.sqrt((l-j)/(k-i)))**2, xf1)
							# (I haven't tested this in the python version)
					h[s] = linkCost+S[k,l]
			# Select the optimal value and record its trail.
			S[i,j] = np.min(h)
			crumbs[i,j] = np.argmin(h)

	# Construct gamma.
	gammaSparse = np.zeros([2,1])
	counter = 0
	while gammaSparse[0,-1]!=grain-1:
		tPoint = int(gammaSparse[0,-1])
		gPoint = int(gammaSparse[1,-1])
		gammaSparse = np.hstack([gammaSparse,np.array([
			[tPoint+steps[int(crumbs[tPoint,gPoint]),0]],
			[gPoint+steps[int(crumbs[tPoint,gPoint]),1]]
			])])
	gammaSparse = list(gammaSparse)
	gammaSparse[0] = [int(a) for a in gammaSparse[0]]
	gammaSparse[1] = [int(a) for a in gammaSparse[1]]
	gamma = np.interp(t, t[gammaSparse[0]],t[gammaSparse[1]])
	gamma = np.round(gamma,8)
	# Warp f2.
	f2Warped = np.interp(gamma, t,f2)
	#
	return f2Warped,gamma

##############################################################################
### >>'''
##############################################################################