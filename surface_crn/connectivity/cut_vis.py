from ase.visualize import view
import ase
from ase.build import fcc111, bcc111
import numpy as np
import cv2
import mayavi.mlab as plt3d
import time

# slabA = fcc111('Al', size=(5,5,5))
# slabB = bcc111('Fe', size=(5,5,5))
#
# view(slabA)
# view(slabB)

xlen, ylen, zlen = 3, 3, 3

def compute(x,y,z):
    unit = np.array([x,y,z])

    atoms = []
    for x in range(xlen):
        for y in range(ylen):
            for z in range(zlen):
                i = x + y*xlen + z*xlen*ylen
                atoms.append([x, y, z])
                atoms.append([x + .5, y + .5, z])
                atoms.append([x, y + .5, z + .5])
                atoms.append([x + .5, y, z + .5])
    atoms = np.array(atoms)

    unit_norm = unit / np.linalg.norm(unit)

    dists = atoms @ unit_norm
    projs = atoms - np.outer(dists, unit_norm)

    ordered_dists = np.sort(dists)
    maxd = np.max(np.abs(ordered_dists[:-1] - ordered_dists[1:]))

    mod_dists = dists % (maxd)
    pruned_idxs = []
    def pcopy(i, proj):
        for pruned_idx in pruned_idxs:
            if np.sum(np.abs(proj - projs[pruned_idx])) < .5:
                return
        pruned_idxs.append(i)

    for i, p in enumerate(projs):
        pcopy(i, p)

    pruned_idxs = np.array(pruned_idxs)

    dmax, dmin = np.max(dists[pruned_idxs]), np.min(dists[pruned_idxs])
    imgx, imgy = 400, 400
    img = np.zeros((imgx, imgy, 3), dtype=np.uint8)
    for i, a in enumerate(atoms[pruned_idxs]):
        x, y, _ = a
        x = ((x + .5) * imgx) / xlen
        y = ((y + .5) * imgy) / ylen
        x, y = int(x), int(y)
        d = (dists[pruned_idxs][i] - dmin) * 200 / (dmax - dmin) + 55
        img = cv2.circle(img, (x, y), 10, (d, d, d), -1)

    # print(atoms.shape, pruned_atoms.shape, dists.shape)
    plt3d.clf()
    time.sleep(.1)
    plt3d.clf()
    plt3d.figure('pruned_atoms')
    plt3d.points3d(*atoms[pruned_idxs].T)
    plt3d.figure('projs')
    plt3d.points3d(*projs.T)
    return img


def na(x): pass
cv2.namedWindow('img')
cv2.createTrackbar('X', 'img', 1, 5, na)
cv2.createTrackbar('Y', 'img', 0, 5, na)
cv2.createTrackbar('Z', 'img', 0, 5, na)
last_x, last_y, last_z = 0, 0, 0
while True:
    x = cv2.getTrackbarPos('X', 'img')
    y = cv2.getTrackbarPos('Y', 'img')
    z = cv2.getTrackbarPos('Z', 'img')
    if x != last_x or y != last_y or z != last_z:
        img = compute(x, y, z)
        last_x, last_y, last_z = x, y, z
    cv2.imshow('img', img)
    k = cv2.waitKey(1) & 0xFF
    if k == ord('q'):
        break
    # k = cv2.waitKey(0)
    # if k == 27:         # wait for ESC key to exit
    #     cv2.destroyAllWindows()
    # elif k == ord('s'): # wait for 's' key to save and exit
    #     cv2.imwrite('shapes.png', img)
    #     cv2.destroyAllWindows()
