"""
Written by Riya Samanta
I am using backbone atoms instead of Calpha atoms
python calculate_volume.py -s topology -f trajectory -i what to plot -p protein -t mutation -d version
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

Goal of the code:
1. Calculate volume of the three helices (28-40,49-61,69-88)
2. Ensure that the structure file is in tpr format: otherwise, indices might be wrong
3. Segregate trajectory according to the surface area criterion.
 
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
parser.add_argument('-i', '--hel', default="H1")
parser.add_argument('-p','--protein',default="S100B_Ca")
parser.add_argument('-t','--protein_type',default="dimer")
parser.add_argument('-d','--d_type',default="combine")
parser.add_argument('-ut','--ut', default="1")
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
hel = str(args.hel)
protein = str(args.protein)
protein_type = str(args.protein_type)
d_type = str(args.d_type)
ut = str(args.ut)
##########################################################################################
#Define functions to calculate volumes
def convex_hull_volume(pts):
	ch= ConvexHull(pts)
	return ch.volume,ch.simplices,ch.area #Volume, simplices , area

##########################################################################################
u = MDA.Universe(top,traj,tpr_resid_from_one=True)
n1,n2=88,9
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
	print(hel_atoms.resids)

elif hel=="EF1_l": #First Ca2+
	hel_atoms=u.select_atoms('protein and ((resid 18 21 23 26 and name O) or (resid 31 and name OE* OD*))')
elif hel=="H1_if": #Helical interface
	if protein=="S100B_Ca":
		hel_atoms=u.select_atoms('protein and backbone and (resid 1-19 or resid 89-107)')
	else:
		hel_atoms=u.select_atoms('protein and backbone and (resid 1-19 or resid 98-116)')
	
hel_volume=[]
hel_area=[]
for t in u.trajectory:
	if hel=="EF2_l_2": # Different definitions of EFII loop coordinating Oxygens around the Ca2+
		hel_atoms=u.select_atoms('protein and name OE* OD* O*  and around 3 (resname CA and resid 177)')
	hel_positions = hel_atoms.positions
	V,simp,A = convex_hull_volume(hel_positions)
	hel_volume.append(V)
	hel_area.append(A)


direc='/home/riya/scratch/S100B/Codes/dis_files/'
np.savetxt(direc+protein_type+'_'+protein+'_Helical_area_backbone_'+hel+'_'+d_type+'_ut_'+ut+'_v4.txt',hel_area)
#np.savetxt(direc+protein_type+'_'+protein+'_Helical_volume_backbone_'+hel+'_'+d_type+'_ut_'+ut+'_v4.txt',hel_volume)
