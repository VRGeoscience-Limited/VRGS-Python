import vrgs
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import pandas as pd

# TODO : Add some Python code

def classify(dip, coplanarity):
    if dip > 80 and coplanarity > 6:
        return 1
    else:
        return 0 
    
def generate_crossplot(azimuth, dip, coplanarity):    
	fig, ax = plt.subplots()
	size = np.asarray(coplanarity)
	cmap = mpl.cm.inferno
	ax.scatter( azimuth,dip, coplanarity, s=size, cmap=cmap, alpha=0.05)
	ax.set_xlabel('Dip', fontsize=15)
	ax.set_ylabel('Azimuth', fontsize=15)
	ax.set_title('Classified')
	ax.grid(True)
	fig.tight_layout()
	plt.show()
    
def generate_intensity_plot(azimuth, dip, coplanarity):  
	fig, ax = plt.subplots()
	cmap = mpl.cm.inferno
	plt.hist2d(azimuth,dip,   bins=36, cmap=cmap)
	ax.set_xlabel('Dip', fontsize=15)
	ax.set_ylabel('Azimuth', fontsize=15)
	ax.set_title('Classified')
	ax.grid(True)
	cb = plt.colorbar()
	cb.set_label('counts in bin')
	fig.tight_layout()
	plt.show()
		
meshes =  vrgs.mesh_list()
if meshes:
	print ("List of meshes",meshes)
	attributes = vrgs.mesh_attribute_list(meshes[0][0])
	if attributes:
		print(meshes[0][0], attributes)
		dip = vrgs.mesh_attribute(meshes[0][0], "Dip")
		azimuth = vrgs.mesh_attribute(meshes[0][0], "Azimuth")
		coplanarity = vrgs.mesh_attribute("Color_prerift_02 - Scan01_TIN", "CoPlanarity")
		normalx = vrgs.mesh_attribute("Color_prerift_02 - Scan01_TIN", "NormalX")
		normaly = vrgs.mesh_attribute("Color_prerift_02 - Scan01_TIN", "NormalY")
		normalz = vrgs.mesh_attribute("Color_prerift_02 - Scan01_TIN", "NormalZ")
		deltacurve = vrgs.mesh_attribute("Color_prerift_02 - Scan01_TIN", "ChangeOfCurvature")
		result = list(map(classify, dip, coplanarity))
		if result:
			vrgs.mesh_add_attribute("Color_prerift_02 - Scan01_TIN", "classified", result)
##		generate_crossplot(azimuth, dip, coplanarity)
		generate_intensity_plot(azimuth, dip, coplanarity)
		