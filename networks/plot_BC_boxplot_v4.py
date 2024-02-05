"""
Written by Riya Samanta

Boxplots for the S100B/Ca and S100B/TRTK/Ca for the n replicas. Currently n = 2*4 replicas (2 subunits in each replica)
And plot them for different concentrations.
Confidence cutoff added.
Only those residues are shown that showed statistical difference. Applies Welch's t-test which is to 
compare two samples with unequal variance and small sample sizes.

"""
#Define modules

import numpy as np
import os
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-bc','--bc', default="gn")
parser.add_argument('-c1', '--cutoff1', default=2) #nearest neighbor
parser.add_argument('-c2', '--cutoff2', default=5) #
parser.add_argument('-m1', '--mincutoff', default=0.2)
parser.add_argument('-m2', '--maxcutoff', default=1.0) #currently testing with 0.9, 1.0
parser.add_argument('-nmax','--nmax', default='yes') #yes,no
parser.add_argument('-conc','--conc1',default='17') #42,71,180,180,240
parser.add_argument('-alpha','--alpha', default=0.1) #0.05,0.1,0.15

#Initialize
args = parser.parse_args()
cutoff2 = float(args.cutoff2)
cutoff1 = int(args.cutoff1)
bc_type = str(args.bc)
mincutoff = float(args.mincutoff)
maxcutoff = float(args.maxcutoff)
nmax = str(args.nmax)
conc1 = str(args.conc1)
alpha = str(args.alpha)
#------------------------------------------------------------
#Define function

def split_list(x): #x is a list
	half = len(x)//2
	return [x[:half],x[half:]]

def reformat(x1,x2,nn=9,prot="S100B_Ca"):
	if prot=="S100B_Ca":
		return list(x1)+[np.nan for x in range(nn)]+list(x2)
	elif prot=="S100B_TRTK_Ca":
		return list(x1)+list(x2)

#Output formatted
def reformat2(x1,x2,nn=9,prot="S100B_Ca"):
	if prot=="S100B_Ca":
		nan_array = np.empty((2*len(conc_array),nn))
		nan_array[:] = np.nan
		return np.column_stack((x1,nan_array,x2))
	elif prot=="S100B_Ca":
		return np.column_stack((x1,x2))

#Output concatenated, placement in the figures formatted
def reformat3(x1,x2,nn=9,prot="S100B_Ca"): #x1, x2 lists, nn size of TRTK, prot = S100B_Ca, S100B_TRTK_Ca
	n1,n2,n3=88,9,2
	if prot=="S100B_Ca":
		res_list=list(range(n1))+list(range(n1+n2,n1+n2+n3))
	elif prot=="S100B_TRTK_Ca":
		res_list=list(range(n1+n2+n3))
	return np.column_stack((x1,x2)),res_list

#------------------------------------------------------------
BC_values = [0 for i in range(99)] # S100B + TRTK + 2 Ca
pt_type = ["dimer"]
conc=conc1
conc_dict = {'17':['v16','v17','v30','v31'],\
		'240':['v18','v19','v33','v39'],\
		'42':['v20','v21','v34','v40'],\
		'120':['v22','v23','v35','v41'],\
		'71':['v24','v25','v36','v42'],\
		'180':['v26','v27','v37','v43']}
if bc_type=="gn":
	if nmax=='yes':
		cutoff=0.25
	else:
		cutoff=0.032
elif bc_type=="cf":
	cutoff=0.06
elif bc_type=="cl":
	cutoff=0.6
elif bc_type=="dg":
	cutoff=0.065
#------------------------------------------------------------
bc_values={}
protein_array=['S100B_Ca','S100B_TRTK_Ca']
for protein in protein_array:
	bc1_list,bc2_list=[],[]
	for version in conc_dict[conc]:
		for protein_type in pt_type:
		
		#Energy-based, FOR PAPER FIGURES
			bc=np.loadtxt("/home/riya/scratch/S100B/energy_analysis/"+protein+"/with_Ca_analysis/"+protein_type+"/"+version+"/pairwise_mean_analysis/"+version+'_'+protein+"_"+protein_type+"_node_"+bc_type+"_betweenness_normalized_nncutoff_"+str(cutoff1)+"_heavy_atoms_cutoff_"+str(cutoff2)+"_minmax"+str(mincutoff)+"_"+str(maxcutoff)+"_v3_different_contact_array.txt")
			bc1,bc2 = bc[:-4],bc[-4:]
			if nmax=='yes':
				bc1,bc2 = bc1/max(bc1),bc2/max(bc2) #protein values, Ca 
			bc1_split,bc2_split = split_list(bc1),split_list(bc2) 
			for x in range(len(bc1_split)):
				bc1_list.append(bc1_split[x])
				bc2_list.append(bc2_split[x])
	bc_values[protein+",protein"]=bc1_list
	bc_values[protein+",Ca"]=bc2_list
#print(bc_values)
#--------------------------------------------------------
#Create Image directory

#-------------------------------------------------------
#Create figure
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import scipy as sp

fig, ax = plt.subplots(figsize=(28,10), dpi= 120, tight_layout=True)
plt.rc('axes', titlesize=40)     # fontsize of the axes title
plt.rc('legend', fontsize=40)    # legend fontsize
xres = range(0,99)
ax.grid(True)
n1,n2,n3=88,9,2
#colors_array = ["pink","tab:blue","green","tab:gray","black"]
colors_array=["tab:cyan","darkgray"]
markers_array=['tab:blue','black']
values_list = []
values2_list= []
elements = []
labels_name_list = []
values2_list = [] 
for i in range(len(protein_array)):
	protein=protein_array[i]
	prot_val,Ca_val = np.median(bc_values[protein+',protein'],axis=0),np.median(bc_values[protein+',Ca'],axis=0) #plot the median values
	prot_std,Ca_std = np.std(bc_values[protein+',protein'],axis=0),np.std(bc_values[protein+',Ca'],axis=0)
	values,res_list2 = reformat3(bc_values[protein+',protein'],bc_values[protein+',Ca'],nn=n2,prot=protein)
	label_name=protein.replace('_','/')+'$^{2+}$'
	values2 = reformat(prot_val,Ca_val,nn=n2,prot=protein)

	print(res_list2)
	elements.append(ax.boxplot(values,positions=res_list2,notch=True,widths=0.5,patch_artist=True,whis=1.5,\
		whiskerprops={'linewidth':3.5,'color':colors_array[i],'linestyle':'--'},\
		boxprops={'facecolor':colors_array[i],'edgecolor':markers_array[i],'alpha':0.3},\
		showmeans=False,medianprops={'marker':'o','markerfacecolor':markers_array[i],'markeredgecolor':markers_array[i],'markersize':12}))

	if protein=='S100B_Ca':
		values_ref = values
			
	else:
		values_mod = np.column_stack((values[:,:n1],values[:,-n3:]))
		p_values = sp.stats.ttest_ind(values_mod,values_ref,axis=0,equal_var=False).pvalue

	labels_name_list.append(label_name)
	values2_list.append(values2)

values_max=np.nanmax(np.asarray(values2_list),axis=0)
res_list3=list(range(n1))+list(range(n1+n2,n1+n2+n3)) #for S100B/Ca2+ listing
text_list=list(range(1,89))+list(range(3,12))+['EFII','EFI']
output = [res_list3[idx] for idx, element in enumerate(p_values) if p_values[idx]<float(alpha)] #find residues with values greater than cutoff

for idx in output:
	ax.text(idx, values_max[idx]+.2, s=text_list[idx], horizontalalignment= 'center', verticalalignment='bottom', fontsize=35,c='black')

#ax.set_title('conc='+conc+'mM,'+chr(945)+'='+alpha+',nn>'+str(cutoff1)+',min='+str(mincutoff)+',max='+str(maxcutoff))
ax.set_title(conc+'mM')

if bc_type=="cf":
	ax.set_ylabel('Current Flow Centrality Score',fontsize=30)
elif bc_type=="gn":
	ax.set_ylabel("$C_{B}$", fontsize=40)

ax.set_xlabel('Residue Number',fontsize=40)
ax.set_xticks([0] + [10*i-1 for i in range(1,9)] + [89,98])
ax.set_xticklabels(['1']+[str(10*i) for i in range(1,9)] + ['TRTK','Ca$^{2+}$'],Fontsize=40)
plt.yticks(fontsize=40)
plt.xticks(fontsize=40)
plt.tick_params(axis='both', which='minor', labelsize=30)
plt.axvline(x = 88,c='k')
plt.axvline(x = 96,c='k')
plt.axhline(y = 0,c='k')
if bc_type=="cl":
	ax.set_ylim(0,1)
else:
	ax.set_ylim(0,0.2)

ax.set_xlim([-1,101])
ax.legend([element["boxes"][0] for element in elements],labels_name_list)
#ax.legend()
if nmax=='yes':
	ax.set_ylim(-0.1,1.4)
	fig.savefig(alpha+'_'+conc1+'_BC_'+bc_type+'_normalized_to_max_cutoff_'+str(cutoff1)+'_heavy_atoms_cutoff2_subunitI'+str(cutoff2)+'_min_'+str(mincutoff)+'_max_'+str(maxcutoff)+'_boxplots.tif')
else:
	fig.savefig(alpha+'_'+conc1+'_BC_'+bc_type+'_normalized_cutoff_'+str(cutoff1)+'_heavy_atoms_cutoff2_subunitI'+str(cutoff2)+'_min_'+str(mincutoff)+'_max_'+str(maxcutoff)+'_boxplots.tif')


plt.close()
