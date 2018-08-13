"""
Library for computing chemical complexity of structures based on chemical framework RDKit.
Contains functions for calculation of:
Bertz complexity index
Whitlock complexity index
Barone complexity index
Total walk count
SMCM - synthetic and molecular complexity 
"""

# Copyright (c) 2018 Milan Vorsilak
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from rdkit import Chem
from rdkit.Chem import GraphDescriptors
import math


def BertzCT(mol, cutoff = 100, dMat = None, forceDMat = 1):
    """
    A topological index meant to quantify "complexity" of molecules.

    Consists of a sum of two terms, one representing the complexity
    of the bonding, the other representing the complexity of the
    distribution of heteroatoms.

    Original article: S. H. Bertz, J. Am. Chem. Soc., vol 103, 3599-3601 (1981)   
    
    Wrapper of GraphDescriptors.BertzCT
    """
    return GraphDescriptors.BertzCT(mol, cutoff, dMat, forceDMat)
    
def WhitlockCT(m, ringval=4, unsatval=2, heteroval=1, chiralval=2, protectingGroups=False):
    """
    A chemically intuitive measure (metric) for molecular complexity is described.

    Default keyword parameter values are taken from the original article. 
    In the artical, the masking of protecting groups is suggested but is not yet implemented in the code (#TODO).

    Original article: H. W. Whitlock, J. Org. Chem., 1998, 63, 7982-7989.
    Benzyls, fenyls, etc. are not treated at all.
    
    >>> mol=Chem.MolFromSmiles("CC(=O)O[C@@H]1C2=C(C)[C@H](C[C@@](O)([C@@H](OC(=O)C)[C@@H]3[C@@]4(CO[C@@H]4C[C@H](O)[C@@]3(C)C1=O)OC(C)=O)C2(C)C)OC(=O)[C@H](O)[C@@H](NC(=O)C)C")
    >>> WhitlockCT(mol)
    67
    
    >>> mol=Chem.MolFromSmiles("C1(C)=CC(C)(C)C2(CCC3)CCC(C)C123")
    >>> WhitlockCT(mol)
    20
    
    >>> mol=Chem.MolFromSmiles("O=C1C=CC(=O)C=C1")
    >>> WhitlockCT(mol)
    14
    """
    # possibly TODO phenyls and others are considered as protecting group thus could be removed (keyword argument removePhenyl)
    mol = Chem.AddHs(m)
    n_bonds = 0
    n_unsat = 0
    n_hetero = 0
    n_chiral = len(Chem.FindMolChiralCenters(mol,includeUnassigned=True))
    n_arom = 0

    bonds = mol.GetBonds()
    atoms = mol.GetAtoms()
    
    for bond in bonds:
        n_bonds+=1
        btype=bond.GetBondType()
        if btype==Chem.rdchem.BondType.DOUBLE:
            n_unsat+=1
            n_bonds+=1
        elif btype==Chem.rdchem.BondType.TRIPLE:
            n_unsat+=2
            n_bonds+=2
        elif btype==Chem.rdchem.BondType.AROMATIC:
            n_arom+=1
    for atom in atoms:
        if atom.GetAtomicNum()!=1 and atom.GetAtomicNum()!=6:
            n_hetero+=1
    
    n_rings = len(bonds)-len(atoms)+len(Chem.GetMolFrags(mol))
    
    complexity = ringval*n_rings+unsatval*n_unsat+heteroval*n_hetero+chiralval*n_chiral
    return complexity


def BaroneCT(m,chiral=False):
    """
    Original article: R. Barone and M. Chanon, J. Chem. Inf. Comput. Sci., 2001, 41 (2), pp 269–272
    Chiral extension: Qi Huang, Lin-LiLi, Sheng-Yong Yang, J. Mol. Graph. Model. 2010, 28 (8), pp 775–787

    Parameter values are hardcoded as in the articles.

    >>> mol=Chem.MolFromSmiles("CC1=CCCC1=O")
    >>> BaroneCT(mol)
    135
    >>> mol=Chem.MolFromSmiles("O=C1CCC2(C)C(=O)CCCC2=C1")
    >>> BaroneCT(mol)
    270
    >>> BaroneCT(mol,True)
    290
    
    >>> mol=Chem.MolFromSmiles("CC1(C)CCCC2(C)C3CCC(C2O)C13")
    >>> BaroneCT(mol)
    288
    
    
    #mol=Chem.MolFromSmiles("c1ccccc1-c1ccc(OCCCCCC(=O)C(=O)OC)cc1")
    
    >>> mol=Chem.MolFromSmiles("O=NCNCCCCc1cccc(OC2=COC=C2)c1")
    >>> BaroneCT(mol,True)
    363
    
    >>> mol=Chem.MolFromSmiles("c1ccccc1C(=O)C1=NC=C(C=N1)CCCCCCNC(=O)C(N)N=O")
    >>> BaroneCT(mol,True)
    512
    
    >>> mol=Chem.MolFromSmiles("N1CNCC1SCC=CC=C(CC(C)=O)Cc1cccc(C)c1")
    >>> BaroneCT(mol,True)
    419
    """
    cmpx = 0
    rinfo = m.GetRingInfo()
    for aring in rinfo.AtomRings():
        cmpx += len(aring)*6
    for atom in m.GetAtoms():
        deg = atom.GetExplicitValence()
        if deg == 1:
            cmpx += 3
        elif deg == 2:
            cmpx += 6
        elif deg == 3:
            cmpx += 12
        else:
            cmpx += 24
        if atom.GetAtomicNum() != 6:
            cmpx += 6
        else:
            cmpx +=3
    if chiral:
        cmpx += 20 * len(Chem.FindMolChiralCenters(m,includeUnassigned=True))        
    return cmpx





_SMCM_ENs = {5:0.851,6:1,7:1.149,8:1.297,9:1.446,15:1.086,16:1.235,17:1.384,35:1.244,53:1.103}

"""
Table of bond
atom    atom    single    double    triple    aromatic
C        C      1.000     0.500     0.333     0.667
C        N      0.857     0.429     0.286     0.571
C        O      0.750     0.375
C        F      0.667
C        P      0.400     
C        S      0.375     0.188               0.250
C        Cl     0.353
C        Br     0.171
C        I      0.113
N        N      0.735     0.367               0.490
N        O      0.643     0.321               0.423
O        S      0.281     0.141
"""
_SMCM_BNs = {Chem.rdchem.BondType.SINGLE:{6:{6:1.0,7:0.857,8:0.75,9:0.667,15:0.4,16:0.375,17:0.353,35:0.171,53:0.113},7:{7:0.735,8:0.643},8:{16:0.281}},
       Chem.rdchem.BondType.DOUBLE:{6:{6:0.5,7:0.429,8:0.375,16:0.188,},7:{7:0.367,8:0.321},8:{16:0.141}},
       Chem.rdchem.BondType.TRIPLE:{6:{6:0.333,7:0.286}},
       Chem.rdchem.BondType.AROMATIC:{6:{6:0.667,7:0.571,16:0.25},7:{7:0.49,8:0.423}}}

def _SMCM_GetBondScore(btype,a1,a2):
    ats = sorted([a1,a2])
    s = 0
    try:
        s = _SMCM_BNs[btype][ats[0]][ats[1]]
    except KeyError:
        pass
    return s
    
    


#SMARTS patterns for SMCM metric which lower structure complexity

_SMCM_SMARTSs = [Chem.MolFromSmarts(patt) for patt in (
        "[$([C!D1]-@C(=O)-@[N!a])]",    # peptides, cyclic
        "[$([C!D1]-!@C(=O)-!@[N!a])]",    # peptides, acyclic
        "[$([#6]=,:[*R!C!c]=,:[#6R]~@[*R!C!c])]",    # C attached to Hetero in ring
        "[$([A!C!H](-!@[#6])-!@[#6])]",    # hetero attached to 2 carbons not in ring
        "[$([#7]=[C!$(*a)]-!@[N!H0])]",    # amidine, guanidine
        "[$([#8]=[#6H0]-!@[N!H2])!$(NC[!#6])]",    # nonterminal amide
        "[#8]=[#6](-!@[NH1,NH0])-!@[N!H2]",    # nonterminal urea
        "[$(O=[CD3]([#6])[#6])!$([!#6]=CC=O)]",    # ketone, not diketo
        "[$([#16X2v2]([#6])[#6])]",    # thio ether
        "[$([#8X2H0][#6]=[#8D1])!$(O~C(~O)~O)]",    # carboxyl ester, not carbonate
        "[$([#8X2v2]([#6])[#6])]",    # ether oxygen
        "C1:c-@C-@N-@C-@C-@O-@C-@c:c-1",    # SMARTS from ICCB's Diversity-Oriented
        "[N,O,C]C(=O)C1=C[C@H](*)C[C@H](O)O1",    # Synthesis approach19
        "C[C@H]1O[C@H]~3O[C@H]~C~2~C~C[C@H]([C@@H]~1)[C@@H]~23",    #  
        "C[C@@H]1O[C@@H]~2O[C@H][C@@H][C@@H]~3~C~C~C~1[C@@H]~23",    #  
        "C[C@H]2C[C@]14~C~C~C~C[C@H]1Oc3cccc(CN2)c34",    #  
        "[#6][C@H]1C[C@@H]([#6])O[C@@H](-a)O1",    #  
        "C2=CC[C@@H]1C(=O)~*~*~C(=O)[C@@H]1[C@@H]2-a",    #  
        "C2=CC[C@@H]1C(=O)~*~C(=O)[C@@H]1[C@@H]2-a",    #  
        "a-[CH,CH2;R0;0*]",    # from ref 17
        "[R;0*]-[CH2R0,NHR0,OR0;0*]-[R]",    #  
        "*-[CD3H,ND2;R0;0*](-a)-a",    #  
        "[a]-&!@[a]",    #  
        "[NR;0*]-[CD3R0;0*](=O)-[R]",    #  
        "[NR;0*]-[CD2R0;0*]-[R]",    #  
        "[NR;0*]-[CD2R0;0*]-[CD2,CD3,OD2,ND2,ND3,aD2,aD3]",    #  
        "a-[NHR0]-[CR0;0*](=O)-[OR0,NR0,0*]",    #  
        "[CR,NR]=[CR]-&!@[a]",    #  
        "[$([#6](-!@[#7])-!@[#7!H0])]",    # NHCN not in ring
        "[$([#6](@[#7])@[#7!H0])]",    # NHCN in ring
        "[$([A!C!H](-@[#6])-@[#6])]",    # hetero attached to 2 carbons in ring
        "[$([#6](=[#8])-[#6]=[!#6])]",    # diketo, keto-Ryl
        "[$([#16X3v4,#16X4v6]([#6])[A])]",    # hypervalent sulfur
        "[$(C(~O)(~O)~O)]"    # carbonate
    )]

def SMCM(m):
    """
    Original article: TK Allu, TI Oprea, J. Chem. Inf. Model. 2005, 45(5), pp. 1237-1243.
    >>> mol=Chem.MolFromSmiles("CN(C)CCOCCOCCN(C)C")
    >>> SMCM(mol)
    17.034
    >>> mol1=Chem.MolFromSmiles("OC(=O)CC(C)(O)CC(O)=O")
    >>> SMCM(mol1)
    21.485
    """
    ats = m.GetAtoms()
    bonds = m.GetBonds()
    a_score = 0
    c=0
    for a in ats:
        if a.GetAtomicNum()==6:
            c+=1
        a_score += _SMCM_ENs[a.GetAtomicNum()]
        
    b_score = 0
    for b in bonds:
        b_score += _SMCM_GetBondScore(b.GetBondType(),b.GetBeginAtom().GetAtomicNum(),b.GetEndAtom().GetAtomicNum())
        
    chiral_a = len(Chem.FindMolChiralCenters(m,includeUnassigned=True))*2
#    mot_score = sum([len(m.GetSubstructMatches(patt)) for patt in SMARTSs])
    mot_score = 0
    for patt in _SMCM_SMARTSs:
        mot_score += len(m.GetSubstructMatches(patt))
            
    # TODO add ring score

    score = a_score + b_score + chiral_a - mot_score 
    
    return score

def _AWC(k,atom,table,neighbors):
    if not table[k].has_key(atom):
        table[k][atom] = 0
        for i in neighbors[atom]:
            table[k][atom] += _AWC(k-1,i,table,neighbors)
    return table[k][atom]
    
def TWC(m):
    """
    Original article: Gerta Rucker and Christoph Rucker, J. Chem. Inf. Comput. Sci. 1993, 33, 683-695
    twc = 1/2 sum(k=1..n-1,sum(i=atoms,awc(k,i)))
    >>> mol=Chem.MolFromSmiles("CC(C)C(C)C(C)CC(C)C")
    >>> TWC(mol)
    21784
    
    >>> mol=Chem.MolFromSmiles("C[Pt](C)(C)(C)(C)CCCC")
    >>> TWC(mol)
    21784
    
    >>> mol=Chem.MolFromSmiles("CC(C)C(C)CC(C)CCCC")
    >>> TWC(mol)
    40145
    
    >>> mol=Chem.MolFromSmiles("CC(C)CCC(C)C(C)CCC")
    >>> TWC(mol)
    40145
    
    /*  not working
        mol=Chem.MolFromSmiles("CC(C)C(C)C(C)(CC)CCC")
        TWC(mol)
        69926
    */
    
    >>> mol=Chem.MolFromSmiles("CCC(C)C(C(C)C)C(C)CC")
    >>> TWC(mol)
    69926
    
    >>> mol=Chem.MolFromSmiles("CCC(C)CC(C)CCCCC")
    >>> TWC(mol)
    31474

    >>> mol=Chem.MolFromSmiles("CCC(C)CCCC(CC)CC")
    >>> TWC(mol)
    31474
    """
    twc = 0
    table = []
    neighbors = []
    atoms = m.GetAtoms()
    for k in xrange(len(atoms)):
        table.append({})
    for atom in atoms:
        neighbors.append([i.GetIdx() for i in atom.GetNeighbors()])
        table[0][atom.GetIdx()]=len(neighbors[atom.GetIdx()])
    
    for k in xrange(len(atoms)-1): 
        for i in atoms:
            twc += _AWC(k,i.GetIdx(),table,neighbors)
    twc /= 2
    try:
        return math.log10(twc)
    except ValueError as e:
        return float("nan")

    
