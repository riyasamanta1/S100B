"""
The goal of this code is to reformulate how energy networks are created. 
input energy file: pairwise_interaction_energy.txt
imput contact list: contact_prob.txt
(different cutoffs tried, with probabilities 0.2<=p<=0.9 for residue-residue interactions, and 0.2<=p<=1
for Ca2+ interactions with rest of the system.
Not different from v3, just changed how to write the output name

"""

#Load modules
import sys,os,argparse
import numpy as np
import networkx as nx
from networkx.algorithms import community
from networkx import edge_betweenness_centrality as eb
from networkx import betweenness_centrality as gn_bc #betweenness centrality
from networkx import current_flow_betweenness_centrality as cf_bc #current flow betweenness centrality
from networkx import closeness_centrality as cl_bc #closeness centrality
from networkx import degree_centrality as dg_bc #degree centrality
import pandas as pd
#--------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('-p', '--protein', default="S100B_Ca") # S100B_Ca, S100B_TRTK_Ca
parser.add_argument('-t', '--protein_type', default="dimer") # monomer, dimer
parser.add_argument('-c', '--cutoff', default=2) # 1,2,3 how many nearest neigbors ignored 
parser.add_argument('-c2', '--cutoff2', default=5) # 7, 7.5, 8, 9
parser.add_argument('-s', '--structure', default="topol.tpr")
parser.add_argument('-f', '--trajectory', default="md.xtc")
parser.add_argument('-cc','--conc', default='17') #combine,17,42,120,240
parser.add_argument('-m1','--mincutoff', default=0.2) #testing with min 20%
parser.add_argument('-m2','--maxcutoff', default=0.9) #cutoff for interactions within the protein
parser.add_argument('-m3','--maxcutoff2', default=1) #cutoff for any interactions that involve Ca2+.
#Initialize
args = parser.parse_args()
protein = str(args.protein) #S100B/Ca or S100B/TRTK/Ca
protein_type = str(args.protein_type) # Dimer/monomer
cutoff = int(args.cutoff) #Number of residues NN >  cutoff
cutoff2 = float(args.cutoff2) # cutoff for the contact list createed using C-alpha
structure = str(args.structure)
trajectory = str(args.trajectory)
conc = str(args.conc)
mincutoff = float(args.mincutoff)
maxcutoff = float(args.maxcutoff)
maxcutoff2 = float(args.maxcutoff2)

if protein=="S100B_Ca":
	nres=180
elif protein=="S100B_TRTK_Ca":
	nres=198


#--------------------------------------------------------------------------------------
#Define functions
tol=1e-15
def w1(A):
	B=abs(A)
	(x,y) = np.shape(B)
	for i in range(x):
		for j in range(y):
			if abs(i-j)>cutoff:
				if B[i,j]>tol:
					B[i,j] = (B[i,j])**(-1.0/3)
				else:
					B[i,j] = 0
			else:
				B[i,j] = 0
	return B

def w2(A):
	B = np.where(abs(A)>tol,abs(A),0)
	x,y = np.shape(B)
	for i in range(x):
		for j in range(y):
			if abs(i-j)<=cutoff:
				B[i,j]=0

	return B

def w3(A):
	B=abs(A)
	(x,y) = np.shape(B)
	for i in range(x):
		for j in range(y):
			if abs(i-j)>cutoff:
				if B[i,j]>tol:
					B[i,j] = -np.log(B[i,j])
				else:
					B[i,j] = 0
			else:
				B[i,j] = 0
	return B


def special_contact_array(A,c1=0.2,c2=0.9,c3=1,Ca=4): #contact probabilities array, mincutoff1,maxcutoff1,maxcutoff2, number of Ca
	sp_contact_array=np.zeros(np.shape(A),dtype=int)
	p_a=A[:-Ca,:-Ca] #between protein residues
	p_c_a=A[:-Ca,-Ca:] # between protein and Ca2+
	c_a=A[-Ca:,-Ca:] # between Ca2+
	protein_array = np.where((p_a>=c1) & (p_a<=c2),1,0)
	protein_Ca_array = np.where((p_c_a>=c1) & (p_c_a<=c3),1,0)
	Ca_array = np.where((c_a>=c1) & (c_a<=c3),1,0)
	sp_contact_array[:-Ca,:-Ca] = protein_array
	sp_contact_array[:-Ca,-Ca:] = protein_Ca_array
	sp_contact_array[-Ca:,:-Ca] = protein_Ca_array.transpose()
	sp_contact_array[-Ca:,-Ca:] = Ca_array
	return sp_contact_array
#----------------------------------------------------------------------------------------



energy=pd.read_csv('/home/riya/scratch/S100B/energy_analysis/'+protein+'/with_Ca_analysis/'+protein_type+'/'+conc+'/pairwise_interaction_energy.txt',sep=" ",header=None,index_col=False)
contact_file='/home/riya/scratch/S100B/energy_analysis/contact_files_for_plot/'+protein+'_'+protein_type+'_contact_prob_heavy_atoms_cutoff_'+str(cutoff2)+'_conc_'+conc+'.txt'

contact_lines=np.loadtxt(contact_file)
#original contact array where all edges follow the same rule
#contact_array=np.where((contact_lines>=mincutoff) & (contact_lines<=maxcutoff),1,0)

#new contact array where 0.2<=p<=0.9 between protein residues. And 0.2<=p<=1 between Ca2+ and rest.
contact_array=special_contact_array(contact_lines,c1=mincutoff,c2=maxcutoff,c3=maxcutoff2,Ca=4)

#create an interaction energy array
energy_array=energy.to_numpy()

#Filter those edges that don't have contacts
energy_contact_array=np.multiply(energy_array,contact_array)
"""
energy_contact_array=np.zeros((nres,nres))
for i in range(nres):
	for j in range(nres):
		energy_contact_array[i,j]  = contact_array[i,j]*energy_array[i,j] 
G=nx.from_numpy_array(w2(energy_array))
output_file=open('energy_final_w2.txt','w+')
for x in G.edges(data=True):
	output_file.write(str(x)+'\n')
"""
#-----------------------------------------------------------------------------
#Extract graph G 
#G1
"""
gn=girvan-newman betweenness centrality, function of weight w1,w3
cl=closeness centrality, w1,w3
dg=degree centrality, function of weight can be anything
cf=current flow centrality, function of weight w2
"""
bw_array=["gn"] 
w_array={"cf":w2,"gn":w1,"gn1":w3,"cl":w1,"dg":w1} #weight functions
bc_array={"cf":cf_bc,"gn":gn_bc,"gn1":gn_bc,"cl":cl_bc,"dg":dg_bc} # betweenness centrality functions
for i in bw_array:
	weight_array = w_array[i](energy_contact_array)
	G=nx.from_numpy_array(weight_array)

	if i=="cf": #Current-Flow Centrality
		Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
		G0 = G.subgraph(Gcc[0])
		G_nodes = bc_array[i](G0,weight='weight',normalized=True)
		initial_node_values=[0 for i in range(nres)]
		for kk in G_nodes.keys():
			initial_node_values[kk] = G_nodes[kk]
		initial_node_values=np.asarray(initial_node_values)

	else: #Girvan-Newman, closeness centrality, degree centrality
		if i=="cl":
			G_nodes = bc_array[i](G,distance='weight')
		elif i=="dg":
			G_nodes = bc_array[i](G)
		else:
			G_nodes = bc_array[i](G,weight='weight',normalized=True)
		initial_node_values = np.asarray(list(G_nodes.values()))


	final_node_values = initial_node_values
	residue_prob_file_str= "/home/riya/scratch/S100B/energy_analysis/"+protein+"/with_Ca_analysis/"+protein_type+"/"+conc+"/pairwise_mean_analysis/"+conc+"_"+protein+"_"+protein_type+"_node_"+i+"_betweenness_normalized_nncutoff_"+str(cutoff)+"_heavy_atoms_cutoff_"+str(cutoff2)+"_minmax"+str(mincutoff)+"_"+str(maxcutoff)+"_v3_different_contact_array.txt"
	#residue_prob_file_str= "/home/riya/scratch/S100B/energy_analysis/"+protein+"/with_Ca_analysis/"+protein_type+"/"+conc+"/pairwise_mean_analysis/"+protein+"_"+protein_type+"_node_"+i+"_betweenness_normalized_nncutoff_"+str(cutoff)+"_heavy_atoms_cutoff_"+str(cutoff2)+"_minmax"+str(mincutoff)+"_"+str(maxcutoff)+"_v3.txt"

	np.savetxt(residue_prob_file_str,final_node_values,fmt="%0.4f")
