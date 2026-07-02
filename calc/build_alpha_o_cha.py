#!/usr/bin/env python3
"""
Build an alpha-O cluster from the REAL IZA CHA framework (via ZSE / ase).

Real crystallographic framework geometry -> a 6-ring (the alpha-site window) is
cut out, dangling Si-O bonds are silanol-capped (Si-O-H), two non-adjacent ring
Si are replaced by Al (Loewenstein), and an Fe(IV)=O is placed along the ring
normal. Atom order: Fe (0), oxo O (1), then the framework — so the v2 pipeline's
atom-indexed AVAS ['0 Fe 3d','1 O 2p'] selects Fe 3d + oxo 2p.

This is the standard literature alpha-O model ([FeO]2+ in a 2-Al 6-ring of CHA,
PNAS 2018 10.1073/pnas.1721717115) but built on REAL framework coordinates — the
O-donor geometry around Fe is crystallographic, not hand-drawn.

Writes calc/alpha_o_cha.xyz and prints diagnostics (atom count, Fe-O distances,
min pairwise distance overlap check). NO QM here.
"""
import numpy as np
from zse.collections import framework
from zse.rings import get_rings

SIO = 1.95   # Si-O bond cutoff (Angstrom)
OH = 0.96
FE_OXO = 1.62


def main():
    cha = framework("CHA")
    sym = cha.get_chemical_symbols()
    Si = [i for i, s in enumerate(sym) if s == "Si"]
    sizes, rings_list, _, big = get_rings(cha, Si[0])
    P = big.get_positions(); S = big.get_chemical_symbols()
    # pick a SPATIALLY COMPACT 12-atom (6T) ring (others are PBC-broken in `big`)
    ring = None
    for r in rings_list:
        if len(r) == 12 and sum(S[j] == "Si" for j in r) == 6:
            pos = P[r]
            mx = max(np.linalg.norm(pos[a] - pos[b]) for a in range(12) for b in range(a + 1, 12))
            if mx < 7.0:
                ring = r; break
    if ring is None:                                # fall back: unwrap the first 6-ring
        ring = next(r for r in rings_list if len(r) == 12 and sum(S[j] == "Si" for j in r) == 6)
    ring0 = ring[0]

    def uw(i):                                       # minimum-image position relative to ring[0]
        return P[ring0] + big.get_distance(ring0, int(i), mic=True, vector=True)

    def nbr(j, elem):                                # mic neighbor search
        d = big.get_distances(int(j), list(range(len(big))), mic=True)
        return [k for k in np.where(d < SIO)[0] if S[k] == elem and k != j]

    ring_Si = [j for j in ring if S[j] == "Si"]
    ring_O = [j for j in ring if S[j] == "O"]
    assert len(ring_Si) == 6 and len(ring_O) == 6, (len(ring_Si), len(ring_O))
    ring_O_set = set(ring_O)

    # Loewenstein: 2 Al = ring-Si pair with max (unwrapped) separation -> non-adjacent
    uwSi = {j: uw(j) for j in ring_Si}
    best = (None, -1)
    for a in range(6):
        for b in range(a + 1, 6):
            dd = np.linalg.norm(uwSi[ring_Si[a]] - uwSi[ring_Si[b]])
            if dd > best[1]:
                best = ((a, b), dd)
    al_idx = set(best[0])

    atoms = []   # (symbol, xyz)
    for n, si in enumerate(ring_Si):
        atoms.append(("Al" if n in al_idx else "Si", uw(si)))
    for o in ring_O:
        atoms.append(("O", uw(o)))
    # external O of each ring Si -> silanol cap (keep O, add H toward its outer Si)
    ext_O = []
    for si in ring_Si:
        for o in nbr(si, "O"):
            if o not in ring_O_set and o not in ext_O:
                ext_O.append(o)
    for o in ext_O:
        atoms.append(("O", uw(o)))
        outer = [k for k in nbr(o, "Si") if k not in ring_Si]
        if outer:
            v = big.get_distance(int(o), int(outer[0]), mic=True, vector=True)
            atoms.append(("H", uw(o) + OH * v / np.linalg.norm(v)))

    # Fe + oxo along the ring normal (centroid of UNWRAPPED ring O), toward +normal
    ringO_pos = np.array([uw(o) for o in ring_O])
    centroid = ringO_pos.mean(axis=0)
    # normal via SVD of centered ring-O positions
    _, _, vh = np.linalg.svd(ringO_pos - centroid)
    normal = vh[2] / np.linalg.norm(vh[2])
    fe = centroid + 0.8 * normal
    oxo = fe + FE_OXO * normal
    cluster = [("Fe", fe), ("O", oxo)] + atoms

    # diagnostics
    xyz = np.array([c[1] for c in cluster]); syms = [c[0] for c in cluster]
    from collections import Counter
    dmin = np.inf
    for i in range(len(xyz)):
        dd = np.linalg.norm(xyz[i + 1:] - xyz[i], axis=1) if i + 1 < len(xyz) else np.array([np.inf])
        if len(dd): dmin = min(dmin, dd.min())
    feO = sorted(np.linalg.norm(ringO_pos - fe, axis=1))
    print("n atoms :", len(cluster), dict(Counter(syms)))
    print("Al sites:", sorted(al_idx), " sep=%.2f A" % best[1])
    print("Fe-oxo  : %.3f A" % np.linalg.norm(oxo - fe))
    print("Fe-(ring O) nearest 3: ", [round(x, 2) for x in feO[:3]])
    print("min pairwise dist: %.3f A %s" % (dmin, "(OK)" if dmin > 0.8 else "(*** OVERLAP ***)"))

    with open("calc/alpha_o_cha.xyz", "w") as f:
        f.write(f"{len(cluster)}\nalpha-O cluster from real IZA CHA 6-ring (Fe,oxo first)\n")
        for s, p in cluster:
            f.write(f"{s} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")
    print("wrote calc/alpha_o_cha.xyz")


if __name__ == "__main__":
    main()
