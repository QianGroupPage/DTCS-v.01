import numpy as np
import cv2
import mayavi.mlab as plt3d
import time

from surface_crn.connectivity.cut_gen import lattice, cut_lattice

def connectivity(lattice: dict):
    pass


def reparamaterize(points, cut, img_size, margin=10):
    basis_0 = np.asarray(cut, dtype=np.float)
    print(np.linalg.norm(basis_0))
    basis_0 /= np.linalg.norm(basis_0)
    basis_1 = np.asarray([0, 0, 1], dtype=np.float)
    basis_1 -= basis_0 * np.dot(basis_0, basis_1)
    basis_2 = np.asarray([0, 1, 0], dtype=np.float)
    basis_2 -= basis_0 * np.dot(basis_0, basis_2) + basis_1 * np.dot(basis_1, basis_2)

    print('Bases:', basis_0, basis_1, basis_2)
    # print('Dots:', np.dot(basis_0, basis_1), np.dot(basis_1, basis_2), np.dot(basis_0, basis_2))
    new_points = points @ np.vstack([basis_1, basis_2]).T
    min, max = np.min(new_points, axis=0), np.max(new_points, axis=0)
    new_points -= min
    new_points *= (np.asarray(img_size) - margin*2) / (max - min)
    new_points += margin
    print(np.min(new_points, axis=0), np.max(new_points, axis=0))
    return new_points.astype(np.int)


def draw_cut(atoms, cut, img_size=(400, 400, 3)):
    lat = cut_lattice(atoms, norm=cut)
    cut_atoms = lat['atoms'][lat['cut_idxs']]
    cut_dists = lat['mod_dists'][lat['cut_idxs']]
    rp_atoms = reparamaterize(cut_atoms, cut, img_size[:1])

    dmax, dmin = np.max(cut_dists), np.min(cut_dists)
    img = np.zeros(img_size, dtype=np.uint8)
    for i, a in enumerate(rp_atoms):
        x, y = a
        d = (cut_dists[i] - dmin) * 200 / (dmax - dmin) + 55
        img = cv2.circle(img, (x, y), 10, (d, d, d), -1)

    # print(atoms.shape, pruned_atoms.shape, dists.shape)
    plt3d.clf()
    time.sleep(.1)
    plt3d.clf()
    plt3d.figure('pruned_atoms')
    plt3d.points3d(*rp_atoms.T, cut_dists.T*100, cut_dists.T, scale_factor=40, scale_mode='none')
    # plt3d.figure('projs')
    # plt3d.points3d(*lat['projs'].T)
    return img


if __name__ == '__main__':
    def na(x):
        pass

    cv2.namedWindow('img')
    cv2.createTrackbar('X', 'img', 1, 5, na)
    cv2.createTrackbar('Y', 'img', 0, 5, na)
    cv2.createTrackbar('Z', 'img', 0, 5, na)
    last_cut = (0, 0, 0)

    atoms = lattice(symmetry='fcc', dims=(3, 3, 3))
    while True:
        cut = cv2.getTrackbarPos('X', 'img'),\
              cv2.getTrackbarPos('Y', 'img'),\
              cv2.getTrackbarPos('Z', 'img')
        if cut != last_cut:
            img = draw_cut(atoms, cut)
            last_cut = cut
            cv2.imshow('img', img)

        k = cv2.waitKey(1) & 0xFF
        if k == ord('q'):
            cv2.destroyAllWindows()
        elif k == ord('s'):
            cv2.imwrite('connectivity.png', img)
            cv2.destroyAllWindows()
