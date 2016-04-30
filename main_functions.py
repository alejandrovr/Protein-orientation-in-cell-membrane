from htmd import *
import math

###---FUNCTIONS---###
#Function to define membrane object
def dummy_leaflet(pmin,pmax,z,spacing=3.0):
    '''Function used to generate planes formed of dummy atoms in a space.
    It requires 3 arguments and a last optional argument can be given.
    These are:
        - pmin: list of two coordinates stating the x,y origin of the plane.
        - pmax: list of two coordinates stating the x,y end of the plane.
        - z: z coordinate of the plane.
    Lastly, the optional argument is:
        -spacing: space between atoms.
    Which is defaulted to 3.0.
    By playing with pmin and pmax, we can adjust the plane's size.'''
    mem=Molecule()
    xatoms=int(((pmax[0]-pmin[0])/spacing))+1
    yatoms=int(((pmax[1]-pmin[1])/spacing))+1
    mylist=list()
    mem.empty(xatoms*yatoms)
    for row in range(yatoms):
        for col in range(xatoms):
            mylist.append([pmin[0]+(col*spacing),pmin[1]+(spacing*row),z]) 
    mem.set('coords', mylist)
    mem.set('beta',0)
    mem.set('resname','DUM')
    mem.set('record','ATOM')
    mem.set('name','DUAT')
    residlist=[i for i in range(xatoms*yatoms)]
    mem.set('resid',residlist)
    mem.set('chain','?')
    return mem

#Function to select aromatic rings within membrane
def select_rings(molecule_o, number_atoms_structure, z_1=None, z_2=None):
    '''Function to select a given structure between two z coordinates.
    It takes two required argument and two optional.
    The first required argument is a molecule object.
    The second one is an integer with the number of atoms the structure has.
    Both of the optional arguments are z coordinates, the first being
    larger than the second.
    It returns a list with every structure that fulfills the atom
    number requirements.'''
    #OBTAIN RINGS FROM AROMATIC RESIDUES
    if z_1 != None:
        z_1 = math.floor(np.asscalar(z_1))
        z_2 = math.floor(np.asscalar(z_2))
        myat = molecule_o.get('coords', "z > " + str(z_1) + " and z < " + str(z_2))
    else:
        myat = molecule_o.get('coords')
    for atom in myat:
        #Checking that any given residue has all atoms of the filtered structure.
        if len(list(atom)) == number_atoms_structure:
            try:
                all_rings.append(list(atom))
            except NameError:
                all_rings = list()
                all_rings.append(list(atom))
    return all_rings

#Function to rotate molecule object around X axis
def rotate(molecule_o, struc_array):
    '''Function to rotate a molecule around a given director vector.
    First, it takes an array with structures' coordinates and finds the
    director vector of the mean plane of them.
    Then, it uses this director vector to rotate the molecule object around the
    X axis.
    It returns the rotated molecule object.'''
    #Placing valid rings into array and creating zero-like matrix.
    coord = np.array(struc_array)
    directors = np.empty([3,int(len(coord)/3)])
    #Solve determinant to get director vector of ring's plane
    for i in range(0,len(coord),3):
        directors[0][int(i/3)] = ((coord[i+1][1]-coord[i][1])*(coord[i+2][2]-coord[i][2]))-((coord[i+1][2]-coord[i][2])*(coord[i+2][1]-coord[i][1]))
        directors[1][int(i/3)] = ((coord[i+1][2]-coord[i][2])*(coord[i+2][0]-coord[i][0]))-((coord[i+1][0]-coord[i][0])*(coord[i+2][2]-coord[i][2]))
        directors[2][int(i/3)] = ((coord[i+1][0]-coord[i][0])*(coord[i+2][1]-coord[i][1]))-((coord[i+1][1]-coord[i][1])*(coord[i+2][0]-coord[i][0]))
    vdir = [np.mean(directors[0]),np.mean(directors[1]),np.mean(directors[2])]

    #ROTATION OF PROTEIN AROUND X AXIS
    angle = np.arccos(vdir[2]/math.sqrt(vdir[0]**2+vdir[1]**2+vdir[2]**2))
    molecule_o.rotateBy(rotationmatrix.rotationMatrix([1, 0, 0],angle+math.pi/2),center= [0,0,0])
    return (molecule_o, angle)

#Function to rotate molecule object with a given angle
def angle_rotate(molecule_o, angle, it, max_it):
    '''Function to rotate a molecule around a given angle.'''
    #Placing valid rings into array and creating zero-like matrix.
    #ROTATION OF PROTEIN AROUND X AXIS
    angle = math.degrees(angle)*(((max_it + 1)/360)*it)
    molecule_o.rotateBy(rotationmatrix.rotationMatrix([1, 0, 0],math.radians(angle)+math.pi/2),center= [0,0,0])
    return (molecule_o, angle)

# Function to compute hydrophobicity of molecule (with respect to current membranes)
def get_hydro(mol1, hydrophobies_dict, scale_dict, solvation):
    '''Function to get residue names and their corresponding z.
    It takes as an argument a molecule object and a dictionary of residues
    with their corresponding hydrophibicities.
    It returns a list of tuples with the structure [(residue name, z coordinate)].'''
    if solvation == True:
        molecule_o= solvate(mol1, pad=10)
        molecule_o.filter('protein within 4 of water')
    else:
        molecule_o=mol1.copy()
    sets = [set(molecule_o.get('resid', 'not water'))]
    chains = set(molecule_o.get('chain'))
    for se in sets:
        my_list = []
        for chain in chains:
            for atom in se:
                try:
                    ret = "resid " + str(int(atom)) + " and name CA and chain " + str(chain)
                    z = molecule_o.get('coords', ret)[2]
                    res_type = molecule_o.get('resname', ret)[0]
                    if res_type in hydrophobies_dict.keys():
                        if scale_dict != None:
                            my_list.append((res_type, z, scale_dict[chain][atom]))
                        else: 
                            my_list.append((res_type, z, 1))
                except IndexError:
                    pass
    return my_list

# Compute best Z for membrane based on hydrophobicity 
def membrane_z(molecule_o, hydrophobies_dict, window, scale_dict, solvation):
    '''This function calculates z coordinates in basis of the minimization of a function.
    First, it gets a list of attributes with coordinates, finds the relative min and max
    taking into account a window variable and scans through the coordinates.
    Lastly, it returns the coordinates that minimize the objective values.
    It takes a molecule object, a dictionary with values and a window size as arguments
    and returns two coordinates.'''
    hydro = get_hydro(molecule_o, hydrophobies_dict, scale_dict, solvation)
    hydro = sorted(hydro, key=lambda x: x[1])
    min_z = (hydro[0][1] - window)
    max_z = (hydro[-1][1] + window)
    global_hydro = []
    for i in np.arange(min_z, max_z, 0.01):
        hydrophobicy = 0
        for res in hydro:
            if res[1] > i and res[1] < int(i + window): #inside membrane
                hydrophobicy = hydrophobicy - hydrophobies_dict[res[0]]*res[2]
            else: #outside membrane
                hydrophobicy = hydrophobicy + hydrophobies_dict[res[0]]*res[2]
        global_hydro.append((i, hydrophobicy))
    global_hydro = sorted(global_hydro, key=lambda x: x[1])
    mem_z = global_hydro[0][0]
    mem_z2 = mem_z + window
    return (mem_z, mem_z2, global_hydro[0][1])

def apply_main(molecule_o, filtering, number_atoms, hydrophobies_dict, window, scale_dict, solvation, max_it = None, style=0):
    '''Function that performs an affine transformation of a molecule object.
    This affine transformation can be looped through the use of the max_it.
    It has three mandatory arguments and four optional ones:
        - It needs a complete molecule object.
        - It needs a molecule object which has been filtered somehow.
        - It needs a number of atoms of the structure we have used to filter.
        - We can pass to it an integer number of iterations to loop through
        the process.
        - We can pass two initial z coordinates. With only 1, it will crash.
        - It takes the current iteration number.'''
    #Generating filtered molecule object
    a = molecule_o.copy()
    a.filter(filtering)
    sel_rings = select_rings(a, number_atoms)
    scans = []
    angle = None
    c_mol = None
    if max_it == None:
       max_it = 0
    it = 0
    while it <= max_it:
       if style == 1 and c_mol != None:
           a = c_mol.copy()
           a.filter(filtering)
           sel_rings = select_rings(a, number_atoms)
       mol = molecule_o.copy()
       it = it + 1
       if style == 1:
           rot_result = rotate(mol, sel_rings)
       if style == 0:
           if angle != None:
               rot_result = angle_rotate(mol, angle, it, max_it)
           else:
               rot_result = rotate(mol, sel_rings)
       c_mol = rot_result[0]
       angle = rot_result[1]
       (c_z_1, c_z_2, c_hydrophobicy) = membrane_z(c_mol, hydrophobies_dict, window, scale_dict, solvation)
       scan = (c_hydrophobicy, c_mol, c_z_1, c_z_2)
       try:
           scans.append(scan)
       except NameError:
           scans = [scan]
    scans = sorted(scans, key=lambda x: x[0])
    return (scans[0][1], scans[0][2], scans[0][3])
