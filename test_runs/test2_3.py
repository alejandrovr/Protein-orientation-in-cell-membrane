from htmd import *
from main_functions import apply_main, dummy_leaflet

molecule_id='1auo'
opm_molecule_path='./examples/'+molecule_id+'.pdb'
conditions = ['resname HIS and (name CE1 or name CD2 or name ND1)',
'(resname TRP or resname PHE or resname TYR) and (name CE2 or name CD2 or name CD1)']
f_structure_a_n = 3
pmin=[-100,-100]
pmax=[100,100]
mem_window = 30
max_it = None
asa = False
solvation = False
superpose = False
my_normalized_hydrophobicies = {'PHE': 1, 'ILE': 0.99, 'TRP': 0.97, 'LEU': 0.97, 'VAL': 0.76, 'MET': 0.74,
                               'TYR': 0.63, 'CYS': 0.49, 'ALA': 0.41, 'THR': 0.13, 'HIS': 0.08, 'GLY': 0,
                                'SER': -0.05, 'GLN': -0.1, 'ARG': -0.14, 'LYS': -0.23, 'ASN': -0.28, 'GLU': -0.31,
                               'PRO': -0.46, 'ASP': -0.55}

mol=Molecule(molecule_id)
for element in conditions:
    try:
        filtering = filtering + " or (" + element + ")"
    except NameError:
        if len(conditions) > 1:
            element = "(" + element + ")"
        filtering = element
if asa == True:
    import asa_config
    import asa_module
    atoms = asa_module.generate_atom_list(mol)
    areas_dict = asa_module.calculate_asa(atoms, asa_config.distance)
    scale_dict = asa_module.normalize_areas(areas_dict)
elif asa == False:
    scale_dict = None

(mol, z_1, z_2) = apply_main(mol, filtering, f_structure_a_n, my_normalized_hydrophobicies, mem_window, scale_dict, solvation, max_it, 1)

mem1 = dummy_leaflet(pmin,pmax,z_1)
mem2 = dummy_leaflet(pmin,pmax,z_2)
allmol=mol.copy()
allmol.append(mem1)
allmol.append(mem2)
allmol.center()
if opm_molecule_path != None:
    opm=Molecule(opm_molecule_path)
    opm.center()
    if superpose != True:
        opm.moveBy([50,0,0])
    allmol.append(opm)
allmol.reps.add(sel='protein', style='NewCartoon',color=0)
allmol.reps.add(sel='resname DUM', style='Licorice')
