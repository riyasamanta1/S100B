"""
Written by Riya Samanta

python segregate_trajectory_surface_area.py -s topology -f trajectory -i what to plot -p protein -t mutation -d version
-i:
"H1": #1-19, 28-40, 69-88 (H1,H2,H4)
        hel_atoms=u.select_atoms('protein and backbone and (resid 1-19 or resid 28-40 or resid 69-88)') 
"H2": #28-40, 49-61, 69-88 (H2,H3,H4)
        hel_atoms=u.select_atoms('protein and backbone and (resid 28-40 or resid 49-61 or resid 69-88)') 
"EF2": # (H3+Loop+H4)
        hel_atoms=u.select_atoms('protein and backbone and (resid 49-88)')
"H3_H4":
        hel_atoms=u.select_atoms('protein and backbone and (resid 49-61 or resid 69-88)')
"EF2_l":
        hel_atoms=u.select_atoms('protein and backbone and (resid 61-72)')
"H1_if": #Helical interface
        if protein=="S100B_Ca":
                hel_atoms=u.select_atoms('protein and backbone and (resid 1-19 or resid 89-107)')
        else:
             	hel_atoms=u.select_atoms('protein and backbone and (resid 1-19 or resid 98-116)')

-p:
S100B_Ca, S100B_TRTK_Ca
-t
dimer, D61A, monomer 
-d:
combine, v9,v10..
-cutoff
any float number, 1200 is used for the current case of surface area. We can also segregate the trajetcory according to the radius of gyration.

Goal of the code:
1. Calculate volume/surface area of the three helices (28-40,49-61,69-88)
2. Ensure that the structure file is in tpr format: otherwise, indices might be wrong
3. Segregate the trajetcory according to the surface area criterion into two states, based on 1200 A2 cutoff.
We want to do this because we want to compare the contact maps of distinct states at different salt concentrations, and calculate Ca2+ coordinaion
for the those states. 

"""
import numpy as np
import MDAnalysis as MDA
import scipy as sp
import argparse
from scipy.spatial import ConvexHull, convex_hull_plot_2d
#########################################################################################
#########################################################################################
#Default settings

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--trajectory', default='dynamic.xtc')
parser.add_argument('-s', '--topology', default='dynamic.tpr')
parser.add_argument('-ti', '--initial_frame', default=0000)
parser.add_argument('-tf', '--final_frame', default=500)
parser.add_argument('-np', '--numoftasks', default=10)
parser.add_argument('-dt', '--skip', default=1)
parser.add_argument('-i', '--hel1', default="H1")
parser.add_argument('-p','--protein',default="S100B_Ca")
parser.add_argument('-t','--protein_type',default="dimer")
parser.add_argument('-d','--d_type',default="combine")
parser.add_argument('-cutoff','--cutoff', default=1200)
##########################################################################################
##########################################################################################
#Initialize
args = parser.parse_args()
top = str(args.topology)
traj = str(args.trajectory)
ti = int(args.initial_frame)
tf = int(args.final_frame)
nt = int(args.numoftasks)
skip = int(args.skip)
hel1 = str(args.hel1)
protein = str(args.protein)
protein_type = str(args.protein_type)
d_type = str(args.d_type)
cutoff = float(args.cutoff)
##########################################################################################
#Define functions to calculate volumes
def convex_hull_volume(pts):
	ch= ConvexHull(pts)
	return ch.volume,ch.simplices,ch.area #Volume, simplices , area

##########################################################################################
#Load the system and trajectory
u = MDA.Universe(top,traj,tpr_resid_from_one=True)
#/home/riya/scratch/S100B/simulations/S100B_Ca/dimer/AMBER/v16/
n1,n2=88,9

def get_hel_atoms(hel=hel1,ut='1'):
	if hel=="H1": #1-19, 28-40, 69-88 (H1,H2,H4)
		hel_atoms=u.select_atoms('protein and backbone and (resid 1-19 or resid 28-40 or resid 69-88)') 
	elif hel=="H2": #28-40, 49-61, 69-88 (H2,H3,H4)
		hel_atoms=u.select_atoms('protein and backbone and (resid 28-40 or resid 49-61 or resid 69-88)') 
	elif hel=="EF2": # (H3+Loop+H4)
		hel_atoms=u.select_atoms('protein and backbone and (resid 49-88)')
	elif hel=="H3_H4":
		numbers=[49,61,69,88]
		if ut=='1':
			r_list=numbers
		elif ut=='2':
			if protein=='S100B_Ca':
				r_list=[int(x+n1) for x in numbers]
			elif protein=='S100B_TRTK_Ca':
				r_list=[int(x+n1+n2) for x in numbers]
		r1,r2,r3,r4=r_list[0],r_list[1],r_list[2],r_list[3]
		hel_atoms=u.select_atoms('protein and backbone and (resid '+str(r1)+'-'+str(r2)+' or resid '+str(r3)+'-'+str(r4)+')')


	elif hel=="EF2_l": #Second Ca2+
		helres=[61, 63, 65, 72, 67]
		if ut=='1':
			helresstr = [str(x) for x in helres]
		elif ut=='2':
			if protein=="S100B_Ca":
				helresstr=[str(x+n1) for x in helres]
			elif protein=="S100B_TRTK_Ca":
				helresstr=[str(x+n1+n2) for x in helres]
		helresstr1,helresstr2=' '.join(helresstr[:-1]),helresstr[-1]
		hel_atoms = u.select_atoms('protein and ((resid '+helresstr1+' and name OE* OD*) or (resid '+helresstr2+' and name O))') 

	elif hel=="EF1_l": #First Ca2+
		hel_atoms=u.select_atoms('protein and ((resid 18 21 23 26 and name O) or (resid 31 and name OE* OD*))')
	elif hel=="H1_if": #Helical interface
		if protein=="S100B_Ca":
			hel_atoms=u.select_atoms('protein and backbone and (resid 1-19 or resid 89-107)')
		else:
			hel_atoms=u.select_atoms('protein and backbone and (resid 1-19 or resid 98-116)')
	return hel_atoms

states=[]
for t in u.trajectory:
	hel_atoms1 = get_hel_atoms(hel=hel1,ut='1') #subunit 1
	hel_atoms2 = get_hel_atoms(hel=hel1,ut='2') #subunit 2
	V1,simp1,A1 = convex_hull_volume(hel_atoms1.positions)
	V2,simp2,A2 = convex_hull_volume(hel_atoms2.positions)
	if (A1>=cutoff) and (A2>=cutoff): #both open
		states.append(4)
	elif (A1<cutoff) and (A2<cutoff): #both closed
		states.append(1)
	elif (A1<cutoff) and (A2>=cutoff): #first closed, second open
		states.append(2)
	elif (A1>=cutoff) and (A2<cutoff): #first open, second closed
		states.append(3)
 
direc='/home/riya/scratch/S100B/Codes/dis_files/'
np.savetxt(direc+protein_type+'_'+protein+'_Helical_area_backbone_'+hel1+'_'+d_type+'_states.txt',states)

