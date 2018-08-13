#!/usr/bin/env python

"""
Library for generating hard-to-synthesize structures based on Molpher algorithm (and RDKit toolkit).
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

import nonpher.complex_lib as cmplx

from rdkit import Chem
from rdkit.Chem import Descriptors

from molpher.core.morphing.operators import *
from molpher.core.morphing import Molpher
from molpher.core import MolpherMol


class ComplexityFilter:
    """
    3 limit class: 99th percentil, 999th permille, max for each (Bertz, Whitlock, Chiral Barone, SMCM)
    complexity metrics and for corresponding MW from ZINC database (2012).
    
    MW bins starts with up to 150Da, then each following bin contains 50Da heavier structures till 600Da and in the last
    bin there are structures heavier than 600Da.
    """
    _table = [
       # 99th percentil
       [[257, 394, 525, 679, 864, 1036, 1261, 1435, 1573, 1826, 2051],          # Bertz
        [14, 17, 20, 23, 26, 30, 34, 38, 41, 51, 67],                           # Whitlock
        [234, 309, 384, 462, 540, 619, 717, 794, 866, 1002, 1174],              # Barone
        [21.7, 28.8, 35.0, 40.2, 44.0, 48.8, 54.9, 60.5, 67.0, 81.7, 107.6],    # SMCM
       ],
       # 999th permille
       [[390, 520, 673, 845, 1043, 1230, 1456, 1644, 1784, 2052, 2540],          # Bertz
        [19, 22, 25, 28, 32, 36, 42, 47, 51, 64, 81],                            # Whitlock
        [281, 360, 438, 522, 607, 695, 795, 884, 968, 1085, 1408],               # Barone
        [21.7, 33.7, 40.1, 46.2, 49.5, 55.2, 63.5, 70.0, 81.3, 100.7, 128.9],    # SMCM
       ],
       # maximum value
       [[487, 782, 1078, 1358, 1675, 1957, 2417, 2654, 2675, 2714, 3210],        # Bertz
        [32, 46, 52, 59, 61, 64, 76, 70, 80, 78, 100],                           # Whitlock
        [418, 566, 702, 822, 1152, 995, 2250, 1500, 3636, 1445, 4668],           # Barone
        [36.0, 50.8, 60.1, 68.1, 79.4, 68.3, 96.9, 100.1, 109.8, 119.2, 156.1],  # SMCM
       ],
    ]

    def __init__(self, limit="999", nOfCmpxs=1):
        """
        Default complexity limit is set on at least 1 exceeding metric on the 999th permille level
        """
        if limit == "99":
            self._limit_index = 0
        elif limit == "999":
            self._limit_index = 1
        elif limit == "max":
            self._limit_index = 2
        else:
            raise ValueError("Invalid value %s" % (str(limit)))

        self._nOfCmpxs = nOfCmpxs

    def isTooComplex(self, rdmol):
        """
        Return if structure is too complex        
        """
        mol = rdmol
        mw = Descriptors.MolWt(mol)
        # bin index
        b = int(mw / 50) - 2
        if b < 0:
            b = 0
        elif b > 10:
            b = 10
        # # of exceeding thresholds
        cmps = [int(cmplx.BertzCT(mol) > ComplexityFilter._table[self._limit_index][0][b]),
                int(cmplx.WhitlockCT(mol) > ComplexityFilter._table[self._limit_index][1][b]),
                int(cmplx.BaroneCT(mol, chiral=True) > ComplexityFilter._table[self._limit_index][2][b]),
                int(cmplx.SMCM(mol) > ComplexityFilter._table[self._limit_index][3][b]),
               ]
        return sum(cmps) >= self._nOfCmpxs


def nonpher(rdmol, max_steps=30):
    """
    Generator of random morphs    
    """
    mol = MolpherMol(rdmol)
    step = 0
    operators = [
        AddAtom(),
        AddBond(),
        ContractBond(),
        InterlayAtom(),
        MutateAtom(),
        RemoveBond(),
        RerouteBond(),
        RemoveAtom(),
    ]
    molpher = Molpher(mol, operators, attempts=1)

    while step < max_steps:
        molpher()
        new_morphs = molpher.getMorphs()
        if len(new_morphs) > 0:
            molpher.reset(new_morphs[0])
            step += 1
            yield new_morphs[0].asRDMol()


def complex_nonpher(smiles, max_steps=30, complexity_filter=ComplexityFilter()):
    """
    Function returning morph which exceeding desired complexity, if no such is generated, returning None

    keyword argument complexity_filter has to be object with method isTooComplex accepting rdkit.Mol and returning
    True or False if the molecule is too complex
    """
    if not complexity_filter.isTooComplex(Chem.MolFromSmiles(smiles)):
        rdstructs = nonpher(smiles, max_steps)
        for s in rdstructs:
            if complexity_filter.isTooComplex(s):
                return Chem.MolToSmiles(s)
    return None


if __name__ == "__main__":
    import argparse
    import sys
    parser = argparse.ArgumentParser(description="Who knows")
    parser.add_argument('-H', "--header", action='store_true')
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout)
    args = parser.parse_args()
    inp = args.infile
    with args.outfile as out:
        if args.header:
            header = inp.readline()
            h_spls = header.split(",")
            out.write("%s,%s,morph_smiles\n" % h_spls[:2])
        for line in inp:
            spls = line.strip().split(",")
            idx, smi = spls[:2]
            morph = complex_nonpher(smi)
            if morph:
                out.write("%s,%s,%s\n" % (idx, smi, morph))