#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
This script computes formation enthalpies (kJ/mol) of molecules provided
as input in the form of SMILES strings or MDL .mol files.
It does not work for compounds such as octalene due to an unsuitable 
definition of bond aromaticity in RDKit.
If the formation enthalpy of a specific entry cannot be calculated, the input
entry is output instead of the value of the formation enthalpy.
"""

import sys, os, fileinput
from collections import Counter
from rdkit import Chem

# atomic heats of formation kJ/mol
HOF = dict(C=717, H=218, Br=112, Cl=121, F=79,I=107, N=473, O=249, S=279)

# additive energetic contributions kJ/mol
BDE = {
    'Br-C': -261.0,
    'C#C': -828.7,
    'C#N': -860.9,
    'C-C': -432.4,
    'C-Cl': -340.1,
    'C-F': -440.5,
    'C-H': -415.5,
    'C-I': -200.9,
    'C-N': -376.4,
    'C-O': -454.5,
    'C-S': -338.1,
    'C=C': -652.2,
    'C=N': -614.7,
    'C=O': -787.4,
    'C=S': -537.6,
    'C~C': -565.2,
    'C~N': -493.1,
    'C~O': -496.2,
    'C~S': -400.9,
    'F-N': -229.7,
    'H-N': -390.8,
    'H-O': -457.7,
    'H-S': -364.5,
    'N-N': -245.2,
    'N-O': -267.9,
    'N=N': -448.6,
    'N=O': -559.0,
    'N~N': -376.0,
    'N~O': -355.0,
    'O-O': -252.1,
    'O-S': -323.4,
    'O=S': -501.6,
    'S-S': -284.3,
    'C..C': 30.2,
    'C..Cl': 12.8,
    'C..F': 8.2,
    'C..H': 14.9,
    'C..N': 23.4,
    'C..O': 37.3,
    'C..S': 23.2,
    'Cl..Cl': 15.9,
    'Cl..N': 25.4,
    'F..F': -34.4,
    'F..H': -8.5,
    'F..N': 15.1,
    'H..N': 17.7,
    'H..O': 31.6,
    'H..S': 12.7,
    'N..N': 16.1,
    'N..O': 8.2,
    'O..O': 18.1,
    'O..S': 35.3,
    'S..S': 14.1,
    '=C=': 31.5,
    'C3or4nitros': 61.2,
    'N->O': -60.5,
    'N3': -72.7,
    'NO2': 212.6,
    'R3': 189.7,
    'R4': 94.8,
    'R4a': 425.2,
    'R5': 11.3,
    'R5a': 53.4,
    'cage4-4-4': -12.5,
    'cage6-6-6': -7.0,
    'fused_aro_aro': 127.1,
    'nAtin3aromaticRings': -148.6,
    'nCyclicAromAmides': 239.2
    }

# ignored geminal interactions
IGNORED = set(['H..H', 'N..S', 'Br..F', 'Br..Cl', 'Cl..F', 'I..I', 'Br..Br',
               'F..I', 'Br..H', 'Cl..H', 'Br..C', 'Br..O', 'H..I', 'C..I',
               'Cl..O', 'F..O', 'I..O', 'Br..N', 'I..N'])

MAX_RING_SIZE = 5  # ring strain is neglected if ring size > 5

# N->O bond is considered in this script as a correction
USED_ATOM_COR = ('N3', 'NO2', 'C3or4nitros', '=C=', 'N->O')

# One-character codes for bond orders
BOC = {Chem.rdchem.BondType.SINGLE: '-',
       Chem.rdchem.BondType.DOUBLE: '=',
       Chem.rdchem.BondType.TRIPLE: '#',
       Chem.rdchem.BondType.AROMATIC: '~'
      }

def get_bond(mol, a, b):
    """ select a static bond in mol.b from two atoms a and b """
    for bond in mol.b:
        if set([a, b]) == set([bond.a1, bond.a2]):
            return bond
    return None

class Ring:
    def __init__(self, mol, atom_indices):
        self.type = None                              # type of ring eg. R4 or R6a
        self.atoms = [mol.a[i] for i in atom_indices] # list of atoms in ring
        self.size = len(atom_indices)
        self.aromatic = get_is_aromatic_ring(mol,self.atoms)
        self.type = 'R'+str(self.size)+('a'*int(self.aromatic))

def are_fused_rings(r1, r2):
    intersec = set(r1.atoms) & set(r2.atoms)
    if len(intersec) != 2:
        return False
    a1, a2 = list(intersec)
    if a2 not in a1.nn:
        return False
    return True

def get_is_aromatic_ring(mol, ring_atoms):
    """ a ring is considered as aromatic if all its bonds are aromatic """
    for i in range(len(ring_atoms)):
        a1, a2 = ring_atoms[i], ring_atoms[i-1]
        b = get_bond(mol, a1, a2)
        if b is None:
            print('no bond between two subsequent atoms in ring')
        if not b.GetIsAromatic():
            return False
    return True

def nAtin3aromaticRings(mol):
    rings = [r for r in mol.rings if r.aromatic]
    nr = len(rings)
    if nr < 3:
        return 0
    n = 0
    for k in range(2, nr):   # loop over all triplets of aromatic rings
        for j in range(1, k):
            for i in range(j):
                r1, r2, r3 = rings[i], rings[j], rings[k]
                intersec = set(r1.atoms) & set(r2.atoms) & set(r3.atoms)
                if len(intersec) == 1:
                    n += 1
    return n

def GetCageCorrections(mol):
    """ return dict countering corrections: cage4-4-4 and cage6-6-6 """
    dic = dict()
    rings = [r for r in mol.rings if not r.aromatic]
    nr = len(rings)
    if nr < 3:
        return dic
    for k in range(2, nr):   # loop over all triplets of aromatic rings
        for j in range(1, k):
            for i in range(j):
                r1, r2, r3 = rings[i], rings[j], rings[k]
                intersec = set(r1.atoms) & set(r2.atoms) & set(r3.atoms)
                if len(intersec) == 1:
                    ringsizes = set(len(r.atoms) for r in (r1, r2, r3))
                    if ringsizes in ({4}, {6}):
                        kk = 'cage'+'-'.join(3*[str(ringsizes.pop())])
                        dic[kk] = dic.get(kk, 0) + 1
    return dic

def getSMARTS(mol, smarts):
    """ return all substructures matching a SMARTS string
    => list of tuples of indexes
    """
    substructure = Chem.MolFromSmarts(smarts)
    return mol.GetSubstructMatches(substructure)

def get_SMARTS_atoms(mol, smarts):
    """ return set of atomic indexes that are
    the first atom in specified substructure
    """
    return set([t[0] for t in getSMARTS(mol, smarts)])

def pairs(items):
    """ return all pairs of items """
    for i in range(1, len(items)):
        for j in range(i):
            yield items[i], items[j]
    return

def sel(atoms, crit, index=None):
    """ selection routine:
    select bonds of an atom from specified bond order: bonds = sel(atom, '~')
    select atoms from a list or a bond or its neighbors: atoms = sel(atom, 'C4-0')
    if index is not None: return list[index] instead of list
    """
    if hasattr(atoms, 'GetSymbol'):  # this is in fact an Atom
        l = [b for b in atoms.bonds if b.order == crit]
        if len(l) > 0:
            if index is None:
                return l
            return l[index]
        atoms = atoms.nn
    elif hasattr(atoms, 'a1'):         # this is in fact a Bond
        atoms = [atoms.a1, atoms.a2]
    l = [a for a in atoms if crit in (a.symbol, a.type, a.arh)]
    if index is None:
        return l
    if len(l) <= index:
        return None
    return l[index]

def nel(atoms, crit):  # return number of atom satisfying criterion 'crit'
    return len(sel(atoms, crit))

def get_mol(molinput):
    """ return RDkit molecule from input string 'molinput' with main attributes:
        a = list of 'static' atoms with their own attributes
        b = list of 'static' bonds with their own attributes
        dica = Counter of all atoms (ie. elements) in the molecule including Hs
        dicb = Counter of all covalent bonds including those involving Hs
        dicc = dict of other features: geminal-ring-structural corrections
    """

    # pylint: disable=too-many-return-statements
    # Required in this case by the number of possible outcomes.

    if os.path.isfile(molinput):
        mol = Chem.MolFromMolFile(molinput)
    else:
        mol = Chem.MolFromSmiles(molinput)
    if mol is None:
        return 'error reading molecule'
    mol = Chem.AddHs(mol)
    netot = sum([a.GetAtomicNum() for a in mol.GetAtoms()])
    if netot % 2 == 1:
        return 'radical'
    mol.name = molinput
    mol.a = [a for a in mol.GetAtoms()]
    if len(mol.a) < 3:
        return 'only_%u_atoms_in_molecule' %len(mol.a)
    mol.b = [b for b in mol.GetBonds()]

    for a in mol.a:
        a.cor = set()
        a.symbol = a.GetSymbol()
        if a.symbol not in HOF:
            return 'no_parameters_for_%s' %a.symbol
        a.ii = [x.GetIdx() for x in a.GetNeighbors()]
        a.nn = [mol.a[i] for i in a.ii]
        a.nH = len([x for x in a.nn if x.GetSymbol() == 'H'])
        a.totdeg = a.GetTotalDegree()
        a.type = a.symbol+str(a.totdeg)
        a.arh = a.type + '-' + str(a.nH)
        a.aromatic = a.GetIsAromatic()
        a.bonds = list()
    for b in mol.b:
        b.cor = set()
        i1, i2 = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        b.a1, b.a2 = mol.a[i1], mol.a[i2]
        b.a1.bonds.append(b)
        b.a2.bonds.append(b)
        s1, s2 = b.a1.symbol, b.a2.symbol
        if s1 > s2:
            b.a1, b.a2 = b.a2, b.a1
            s1, s2 = s2, s1
        b.order = BOC[b.GetBondType()]
        b.s1, b.s2 = s1, s2
        b.type = s1 + b.order + s2

    # Ensure correct bond orders in nitro NO2 and azide N3 groups
    for a in mol.a:
        if a.type == 'N3' and nel(a, 'O1-0') == 2: # N in nitro
            for O1 in sel(a, 'O1-0'):
                bo = get_bond(mol, a, O1)
                bo.order = '='
                bo.type = bo.s1 + '=' + bo.s2
        elif a.arh == 'N2-0':        # azide -N3 or isocyanide -NC
            Nterm = sel(a, 'N1-0', 0)
            Nfirst = sel(a, 'N2', 0)
            if Nterm and Nfirst:
                bN3_1, bN3_2 = get_bond(mol, a, Nterm), get_bond(mol, a, Nfirst)
                bN3_1.order = bN3_2.order = '='
    
    # count atoms and bonds
    mol.dica = Counter(a.symbol for a in mol.a)
    mol.dicb = Counter(b.type for b in mol.b)

    # add atomic and bond structural corrections
    mol.dicc = dict()
    for cor in USED_ATOM_COR:
        mol.dicc[cor] = len([a for a in mol.a if cor in a.cor])
    for cor in set.union(*(b.cor for b in mol.b)):
        mol.dicc[cor] = len([b for b in mol.b if cor in b.cor])

    # add geminal interactions to dic
    for a in mol.a:
        for p in pairs(a.nn):
            a1, a2 = p
            if a1 in a2.nn:  # a1 bonded to a2
                continue     # 3-membered ring
            if a1.symbol > a2.symbol:
                a1, a2 = a2, a1
            geminal = a1.symbol + '..' + a2.symbol
            mol.dicc[geminal] = mol.dicc.get(geminal,0) + 1

    # add ring corrections to dic
    ri = mol.GetRingInfo()
    mol.rings = [Ring(mol, r) for r in ri.AtomRings()]
    for ring in mol.rings: # + mol.ringtypes:
        if ring.size > MAX_RING_SIZE:
            continue
        mol.dicc[ring.type] = mol.dicc.get(ring.type, 0) + 1

    # fused rings
    for r1, r2 in pairs(mol.rings):
        if are_fused_rings(r1, r2) and r1.aromatic and r2.aromatic:
            mol.dicc['fused_aro_aro'] = mol.dicc.get('fused_aro_aro', 0) + 1

    # cumulenic C
    mol.dicc['=C='] = len(get_SMARTS_atoms(mol, 'C(=*)=*'))

    # azide
    mol.dicc['N3'] = len(get_SMARTS_atoms(mol, '[NX2]~[NX2]~[NX1]'))
    #[$(*-[NX2-]-[NX2+]#[NX1]),$(*-[NX2]=[NX2+]=[NX1-])]

    # amides involved in aromatic rings
    mol.dicc['nCyclicAromAmides'] = len(get_SMARTS_atoms(mol, '[cR](=[OX1])[#7]'))

    # nitro groups
    NO2 = '[$([NX3](=O)=O),$([NX3+](=O)[O-])]'
    mol.dicc['NO2'] = len(get_SMARTS_atoms(mol, NO2))

    # C atoms bonded to 3 or 4 nitros
    s = 'C(%s)(%s)%s' %(NO2, NO2, NO2)
    mol.dicc['C3or4nitros'] = len(get_SMARTS_atoms(mol, s))

    # N-oxide: not azoxy, nitro, nitroso, nor nitrate.
    NOx = '[$([#7+][OX1-]),$([#7v5]=[OX1]);!$([#7](~[O])~[O]);!$([#7]=[#7])]'
    NOx = '[$([#7X3]~[OX1]);!$([NX3](~[OX1])~[OX1])]'
    mol.dicc['N->O'] = len(get_SMARTS_atoms(mol, NOx))

    # atoms belonging to 3 aromatic rings
    mol.dicc['nAtin3aromaticRings'] = nAtin3aromaticRings(mol)
    # other cage corrections
    cage = GetCageCorrections(mol)
    mol.dicc.update(cage)
    return mol

#### MAIN PROGRAM #####################################################

for molinput in fileinput.input():
    print ("*************************************")
    print ("Molecule is ", molinput)
    molec = get_mol(molinput.rstrip())
    if molec is None:
        print('None')
        continue
    try:
        # are there missing parameters ?
        missing = set(molec.dica) - set(HOF)  # elements beyond scope of the method
        missing |= (set(molec.dicb) - IGNORED) - set(BDE)
        missing |= (set(molec.dicc) - IGNORED) - set(BDE)
        if missing:
            txt_list = ' '.join(sorted(list(missing)))
            print('None. missing: %s' %(txt_list))
            continue
        hof = sum([HOF[k]*molec.dica[k] for k in molec.dica])
        hof += sum([BDE[k]*molec.dicb[k] for k in molec.dicb])
        hof += sum([BDE[k]*molec.dicc[k] for k in molec.dicc if k not in IGNORED])
        print ("Enthalpy formation (kJ/mol)")
        print('%.1f' %hof)
    except AttributeError:
        print('%s' %molec.name)
