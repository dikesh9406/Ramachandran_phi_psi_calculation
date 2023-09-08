#python program to calculate phi and psi dihedral angles

f=open('1asy.pdb')
x=f.readlines()
f.close()
y=[]
for i in x:
	if i.startswith('ATOM'):
		y.append(i)
		
class atom:
	def __init__(self, line):
		#atom name , residue name, residue number, x, y, z we need
		#atom name always consider 3 letter length (nmaximum)
		self.atomnum=int (line[6:11].strip())
		self.atomname=(line[11:17].strip())
		self.resname=(line[17:21].strip())
		self.chain=(line[21].strip())
		self.resnum=int(line[22: 28].strip())
		self.x=float(line[30:38].strip())
		self.y=float(line[38:46].strip())
		self.z=float(line[46:54].strip())
				
				
atomlist=[]
for i in y:
	atomlist.append(atom(i))				

		
		
proobjs=[]
rnaobjs=[]
# only protein
for i in atomlist:
	if len(i.resname)==3:
		proobjs.append(i)
	else:
		rnaobjs.append(i)	




class residue:
	def __init__(self, listofatomobjs):
		self.resname=listofatomobjs[0].resname
		self.resnum=listofatomobjs[0].resnum
		self.chain=listofatomobjs[0].chain
		self.phy=None
		self.psi=None
		self.n=None
		self.ca=None
		self.c=None
		for i in  listofatomobjs:
			if i.atomname=='N':
				self.n=i
				continue
			if i.atomname=='CA':
				self.ca=i
				continue
			if i.atomname=='C':
				self.c=i
					
		
		
		
residues=[]
temp=[]
for i in range(len(proobjs)):
	temp.append(proobjs[i])
	if i==len(proobjs)-1 or proobjs[i].resnum!=proobjs[i+1].resnum:
		residues.append(residue(temp))
		temp=[]



import math
def dihedral(atom1, atom2, atom3, atom4):
	v1=[atom2.x-atom1.x, atom2.y-atom1.y, atom2.z-atom1.z]
	v2=[atom3.x-atom2.x, atom3.y-atom2.y, atom3.z-atom2.z]
	v3=[atom4.x-atom3.x, atom4.y-atom3.y, atom4.z-atom3.z]
	n1=[v1[1]*v2[2]-v2[1]*v1[2], -1*(v1[0]*v2[2]-v2[0]*v1[2]), v1[0]*v2[1]-v2[0]*v1[1]]
	
	n2=[v2[1]*v3[2]-v3[1]*v2[2], -1*(v2[0]*v3[2]-v3[0]*v2[2]), v2[0]*v3[1]-v3[0]*v2[1]]
	
	num=n1[0]*n2[0]+n1[1]*n2[1]+n1[2]*n2[2]
	den=math.sqrt((n1[0]**2+n1[1]**2+n1[2]**2)*(n2[0]**2+n2[1]**2+n2[2]**2))
	angle=math.acos(num/den)*(180/math.pi)
	n3=[n1[1]*n2[2]-n2[1]*n1[2], -1*(n1[0]*n2[2]-n2[0]*n1[2]), n1[0]*n2[1]-n2[0]*n1[1]]
	costheta=(v2[0]*n3[0]+v2[1]*n3[1]+v2[2]*n3[2])/math.sqrt((v2[0]**2+v2[1]**2+v2[2]**2)*(n3[0]**2+n3[1]**2+n3[2]**2))
	if costheta<0:
		angle*=-1
	
	
	return round(angle, 3)




for i in range(len(residues)):
	if i==0 or residues[i-1].chain!=residues[i].chain:
		residues[i].phi=180
		residues[i].psi=dihedral( residues[i].n, residues[i].ca, residues[i].c, residues[i+1].n,)
	elif i==len(residues)-1 or residues[i+1].chain!=residues[i].chain:
		residues[i].psi=180
		residues[i].phi=dihedral(residues[i-1].c, residues[i].n, residues[i].ca, residues[i].c)			
	else:
		residues[i].phi=dihedral(residues[i-1].c, residues[i].n, residues[i].ca, residues[i].c)
		residues[i].psi=dihedral( residues[i].n, residues[i].ca, residues[i].c, residues[i+1].n,)
		

phis=[]
psis=[]
for i in residues:
		phis.append(i.phi)
		psis.append(i.psi)
   

#print('number phi psi')
#for i in range(len(phis)):
#	print(str(i)+'\t' str(phis[i]))		


from matplotlib import pyplot as plt

plt.title("Ramachandran plot", fontsize=15)
plt.xlabel('phi', fontsize=13)
plt.ylabel('psi', fontsize=13)

plt.axis([-180, 180, -180, 180])
plt.axhline(0, color='black', lw=1)
plt.axvline(0, color='black', lw=1)
plt.grid(True)
plt.plot(phis, psis, 'co')
plt.savefig('1asy_rmc.png')
plt.show()



