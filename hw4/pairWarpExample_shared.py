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
#%% Contents 
##############################################################################

## << ExampleFromGithub>
## << PairwiseWarping>
##		< Function F warping
##		< SRVF Q warping

##############################################################################
##############################################################################