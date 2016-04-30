from asa_config import *
from collections import defaultdict
import numpy as np
import math

#--Adaptation of code by Bosco K. Ho's found in 'http://boscoh.com/protein/asapy.html'  

###---CLASSES---###--CUSTOM
class vector(object):
    ##--ATTRIBUTE INITIALIZATION--##
    def __init__(self, coords):
        self.__coordinates = coords
        for coord in self.__coordinates:
            if not isinstance(coord, int) and not isinstance(coord, np.float32) and not isinstance(coord, float):
                raise ValueError("Impossible to create object: %s not possible." % type(coord))
        magnitude2 = 0
        for k, e in enumerate(self.__coordinates):
            magnitude2 = magnitude2 + e*e
        self.__magnitude2 = magnitude2
        self.__magnitude = math.sqrt(magnitude2)
    ##--METHODS--##
    def __iter__(self):
        return iter(self.__coordinates)
    def get_coords(self):
        return self.__coordinates
    def get_magnitude(self):
        return self.__magnitude
    def get_magnitude2(self):
        return self.__magnitude2
    def get_n(self, n):
        e = "Input value for method 'get_n' must be an integer within list scope."
        if type(n) is int:
            try:
                return self.__coordinates[n]
            except IndexError:
                raise ValueError(e)
        else:
            raise ValueError(e)
    def scale(self, value):
        new_coords = []
        e = "Invalid input type for method: vector scaling."
        if type(value) is int or type(value) is float or (type(value) is list and len(value) == 1):
            for k, i in enumerate(self.__coordinates):
                new_coords.append(i * value)
            return new_coords
        else:
            raise ValueError(e)
    def __add__(self, value):
        new_coords = []
        e = "Invalid input type for method"
        if type(value) is int or type(value) is float:
            for k, i in enumerate(self.__coordinates):
                new_coords.append(i + value)
            return new_coords
        elif type(value) is list or type(value) is vector:
            size1 = len(self.__coordinates)
            if type(value) is vector:
                value = value.get_coords()
            size2 = len(value)
            if size1 >= size2:
                new_coords = self.__coordinates
                short_coords = value
            else:
                new_coords = value
                short_coords = self.__coordinates
            for k in range(len(short_coords)):
                new_coords[k] = new_coords[k] + short_coords[k]
            return new_coords
        else:
            raise ValueError(e)
    def __sub__(self, value):
        new_coords = []
        e = "Invalid input type for method"
        if type(value) is int or type(value) is float:
            for k, i in enumerate(self.__coordinates):
                new_coords.append(i - value)
            return new_coords
        elif type(value) is list or type(value) is vector:
            size1 = len(self.__coordinates)
            if type(value) is vector:
                value = value.get_coords()
            size2 = len(value)
            if size1 >= size2:
                new_coords = self.__coordinates
                short_coords = value
            else:
                new_coords = value
                short_coords = self.__coordinates
            for k in range(len(short_coords)):
                new_coords[k] = new_coords[k] - short_coords[k]
            return new_coords
        else:
            raise ValueError(e)
    def __len__(self):
        return len(self.__coordinates)
        
            
class atom(object):
    ##--ATTRIBUTE INITIALIZATION--##
    def __init__(self, identifier, name, vector, residue_id=None, chain=None):
        self.__identifier = identifier
        self.__name = name
        self.__coordinates = vector
        for coord in self.__coordinates:
            if not isinstance(coord, int) and not isinstance(coord, np.float32) and not isinstance(coord, float):
                raise ValueError("Impossible to create object: %s not possible." % coord)
        try:
            self.__radius = my_atomic_radii[name]
        except KeyError:
            self.__radius = 1
        self.__residue_id = residue_id
        self.__chain = chain
    ##--METHODS--#
    def get_identifer(self):
        return self.__identifier
    def get_name(self):
        return self.__name
    def get_coords(self):
        return self.__coordinates.get_coords()
    def get_n(self, n):
        return self.__coordinates.get_n(n)
    def get_radius(self):
        return self.__radius
    def get_residue_id(self):
        return self.__residue_id
    def get_chain(self):
        return self.__chain

###---ASA CALCULUS---###-FULLY ADAPTED
         
#Generates points in a sphere using 'Golden Section Spiral' algorithm
def generate_sphere_points(n):
    points = []
    inc = math.pi*(3 - math.sqrt(5))
    offset = 2/float(n)
    for k in range(int(n)):
        y = k*offset - 1 + (offset/2)
        r = math.sqrt(1 - y*y)
        phi = k*inc
        points.append(vector([math.cos(phi)*r, y, math.sin(phi)*r]))
    return points

def make_boxes(a, d_max):
    b = defaultdict(list) # space divided into boxes
    for i in range(len(a)):
        atom = a[i]
        box_coor = tuple(int(math.floor(x / d_max)) for x in atom._atom__coordinates)
        b[box_coor].append(i)
    return b

def add_bond(a, a1, a2, conn, d_max):
    '''
    If distance between atoms a1 and a2 is less than d_max (neighboring atoms),
    add atoms a1 and a2 in adjacency list conn to each other
    '''
    atom1 = a[a1]
    atom2 = a[a2]
    m = vector(atom1._atom__coordinates - atom2._atom__coordinates)
    if m._vector__magnitude2 <= d_max * d_max:
        conn[a1].append(a2)
        conn[a2].append(a1)

def neighbor_atoms(b, box):
    '''
    Returns list of atoms from half of neighbouring boxes of the box
    another half is accounted when symmetric (opposite) boxes considered
    '''
    na = [] # list for neighboring atoms
    x, y, z = box # coordinates of the box
    # top layer consisting of 9 boxes
    if (x + 1, y + 1, z +1) in b: na.extend(b[(x + 1, y + 1, z +1)])
    if (x, y + 1, z +1) in b: na.extend(b[(x, y + 1, z +1)])
    if (x + 1, y, z +1) in b: na.extend(b[(x + 1, y, z +1)])
    if (x, y, z +1) in b: na.extend(b[(x, y, z +1)])
    if (x - 1, y + 1, z +1) in b: na.extend(b[(x - 1, y + 1, z +1)])
    if (x + 1, y - 1, z +1) in b: na.extend(b[(x + 1, y - 1, z +1)])
    if (x, y - 1, z +1) in b: na.extend(b[(x, y - 1, z +1)])
    if (x - 1, y, z +1) in b: na.extend(b[(x - 1, y, z +1)])
    if (x - 1, y - 1, z +1) in b: na.extend(b[(x - 1, y - 1, z +1)])
    # half of the middle layer excluding the box itself (4 boxes)
    if (x + 1, y + 1, z) in b: na.extend(b[(x + 1, y + 1, z)])
    if (x, y + 1, z) in b: na.extend(b[(x, y + 1, z)])
    if (x + 1, y, z) in b: na.extend(b[(x + 1, y, z)])
    if (x + 1, y - 1, z) in b: na.extend(b[(x + 1, y - 1, z)])
    return na

def adjacency_list(a, d_max):
    '''
    Returns adjacency list from coordinate file
    in O(len(a)) time
    '''
    b = make_boxes(a, d_max) # put atoms into the boxes with dmax length side
    conn = [[] for i in range(len(a))] # list of bond lengths each atom implicated
    for box in b:
        lb = len(b[box])
        for i in range(lb):
            a1 = b[box][i]
            # check possible connections inside the box
            for j in range(i+1, lb):
                a2 = b[box][j]
                add_bond(a, a1, a2, conn, d_max)
            # check connections with atoms from neighbouring boxes
            na = neighbor_atoms(b, box) # list of such atoms
            for a2 in na:
                add_bond(a, a1, a2, conn, d_max)
    return conn

def find_neighbor_indices(atoms, indices, probe, k):
    """
    Returns list of indices of atoms within probe distance to atom k. 
    """
    neighbor_indices = []
    atom_k = atoms[k]
    radius = atom_k._atom__radius + probe + probe
    for i in indices:
        if i == k: continue
        atom_i = atoms[i]
        m = vector(atom_k._atom__coordinates - atom_i._atom__coordinates)
        dist2 = m._vector__magnitude2 # ToAn
        if dist2 < (radius + atom_i._atom__radius) ** 2: # ToAn
            neighbor_indices.append(i)
    return neighbor_indices

def calculate_asa(atom_list, probe, n_sphere_point=960):
    """
    Returns the accessible-surface areas of the atoms, by rolling a
    ball with probe radius over the atoms with their radius
    defined.
    """
    sphere_points = generate_sphere_points(n_sphere_point)
    
    const = 4.0 * math.pi / len(sphere_points)
    areas = defaultdict(dict)
    neighbor_list = adjacency_list(atom_list, 2 * (probe + max(atom_list, key=lambda p: p._atom__radius)._atom__radius))
    for i, atom_i in enumerate(atom_list):

        neighbor_indices = [neig for neig in neighbor_list[i]]
        neighbor_indices = find_neighbor_indices(atom_list, neighbor_indices, probe, i)
        n_neighbor = len(neighbor_indices)
        j_closest_neighbor = 0
        radius = probe + atom_i._atom__radius

        n_accessible_point = 0
        for point in sphere_points:
            is_accessible = True
            surface_point = vector(point)
            test_point = atom_i._atom__coordinates + vector(surface_point.scale(radius))
            cycled_indices = list(range(j_closest_neighbor, n_neighbor))
            cycled_indices.extend(list(range(j_closest_neighbor)))
      
            for j in cycled_indices:
                atom_j = atom_list[neighbor_indices[j]]
                r = atom_j._atom__radius + probe
                test_vector = vector(atom_j._atom__coordinates - test_point)
                diff2 = test_vector._vector__magnitude2
                if diff2 < r*r:
                    j_closest_neighbor = j
                    is_accessible = False
                    break
            if is_accessible:
                n_accessible_point += 1
    
        area = const*n_accessible_point*radius*radius
        try:
            areas[atom_i._atom__chain][atom_i._atom__residue_id] = areas[atom_i._atom__chain][atom_i._atom__residue_id] + area
        except:
            areas[atom_i._atom__chain][atom_i._atom__residue_id] = area

    return areas

##--NEW FUNCTIONS--##--CUSTOM
#Function to generate a list of atom objects from a molecule object.
def generate_atom_list(molecule_o):
    atom_list = []
    sets = [set(molecule_o.get('serial', 'not water'))]
    for se in sets:
        for serial in se:
            ret = "serial " + str(int(serial))
            resid = molecule_o.get('resid', ret)[0]
            chain = molecule_o.get('chain', ret)[0]
            name = molecule_o.get('name', ret)[0]
            coords = molecule_o.get('coords', ret)
            atom_o = atom(serial, name, vector([coords[0], coords[1], coords[2]]), resid, chain)
            atom_list.append(atom_o)
    return atom_list

#Function to normalize values of a dictionary of dictionares.
def normalize_areas(dictionary):
    l=list(dictionary.values())
    s=list()
    s2=list()
    for l2 in l:
        s.append(max(l2.values()))
    max_size = max(s)
    for chain in dictionary.keys():
        for residue in dictionary[chain].keys():
            dictionary[chain][residue] = dictionary[chain][residue]/max_size
    return dictionary
