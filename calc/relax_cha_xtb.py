#!/usr/bin/env python3
"""Constrained GFN2-xTB relaxation of the active site in the real-CHA alpha-O
cluster: framework (Si/Al/caps) FROZEN at crystallographic positions, only
Fe + oxo + the 6 ring O relaxed. Pre-optimizes the Fe=O coordination before the
CASSCF/NEVPT2 spin-state QM. Atom order: 0=Fe, 1=oxo, 2-7=ring T, 8-13=ring O."""
import numpy as np
from ase.io import read, write
from ase.constraints import FixAtoms
from ase.optimize import BFGS
from tblite.ase import TBLite

atoms = read("calc/alpha_o_cha.xyz")
free = {0, 1, 8, 9, 10, 11, 12, 13}          # Fe, oxo, 6 ring O
atoms.set_constraint(FixAtoms(indices=[i for i in range(len(atoms)) if i not in free]))
for kw in (dict(method="GFN2-xTB", charge=0.0, multiplicity=5, verbosity=0),
           dict(method="GFN2-xTB", charge=0.0, verbosity=0),
           dict(method="GFN2-xTB")):
    try:
        atoms.calc = TBLite(**kw); atoms.get_potential_energy(); print("xTB kwargs:", kw); break
    except Exception as e:
        print("xTB init retry:", type(e).__name__, str(e)[:60])

opt = BFGS(atoms, logfile="calc/relax_cha.log")
opt.run(fmax=0.10, steps=250)

p = atoms.get_positions()
feO = sorted(np.linalg.norm(p[8:14] - p[0], axis=1))
print("Fe-oxo   : %.3f A" % np.linalg.norm(p[1] - p[0]))
print("Fe-ringO : ", [round(x, 2) for x in feO])
write("calc/alpha_o_cha_relaxed.xyz", atoms)
print("wrote calc/alpha_o_cha_relaxed.xyz")
