import numpy as np


def lattice(symmetry='fcc', dims=(3, 3, 3)):
    assert symmetry in ['fcc', 'bcc', 'sc'], 'Invalid symmetry chosen'
    ylen, xlen, zlen = dims
    atoms = []

    for x in range(xlen):
        for y in range(ylen):
            for z in range(zlen):
                atoms.append([x, y, z])
                if symmetry == 'fcc':
                    atoms.append([x + .5, y + .5, z])
                    atoms.append([x, y + .5, z + .5])
                    atoms.append([x + .5, y, z + .5])
                elif symmetry == 'bcc':
                    atoms.append([x + .5, y + .5, z + .5])
                elif symmetry == 'sc':
                    pass
    return {'atoms': np.array(atoms)}


def cut_lattice(lattice: dict, norm=(1, 0, 0)):
    assert 'atoms' in lattice
    norm = np.asarray(norm)

    unit_norm = norm / np.linalg.norm(norm)

    lattice['dists'] = lattice['atoms'] @ unit_norm
    lattice['projs'] = lattice['atoms'] - np.outer(lattice['dists'], unit_norm)

    ordered_dists = np.sort(lattice['dists'])
    print(ordered_dists)
    maxd = np.max(np.abs(ordered_dists[:-1] - ordered_dists[1:]))

    lattice['mod_dists'] = lattice['dists'] % (maxd * 2)
    cut_idxs = []

    for i, proj in enumerate(lattice['projs']):
        j = 0
        while j < len(cut_idxs):
            cut_idx = cut_idxs[j]
            cut_proj = lattice['projs'][cut_idx]
            print(proj - cut_proj)
            if np.sum((proj - cut_proj)**2) < 1:
                if lattice['dists'][i] < lattice['dists'][j]:
                    cut_idxs.remove(cut_idx)
                    cut_idxs.append(i)
                    break
            j += 1
        else:
            cut_idxs.append(i)

    lattice['cut_idxs'] = np.array(cut_idxs)
    return lattice
