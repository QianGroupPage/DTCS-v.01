import numpy as np
import cv2
import mayavi.mlab as plt3d
import time
import ase
from ase.visualize import view
from ase.build import fcc100, fcc110, fcc111
from ase.build import bcc100, bcc110, bcc111
import matplotlib.pyplot as plt


def connectivity(positions: np.ndarray, max_dist=None, deg=2, max_deg=3):
    position_dict = {deg - 1: positions}

    n = positions.shape[0]
    positions_matrix = np.expand_dims(positions, axis=0).repeat(n, axis=0)
    dists = positions_matrix - positions_matrix.swapaxes(0, 1)
    dists[..., 2] = positions_matrix[..., 2] + positions_matrix.swapaxes(0, 1)[..., 2]
    dists = np.sqrt(np.sum(dists**2, axis=2))
    dist_bins = np.sort(np.unique(dists))

    # print(dist_bins)
    max_dist = max_dist if max_dist else dist_bins[2] + .01
    if max_dist < dist_bins[1]:
        return {}, position_dict
    # print(dists.astype(np.float16))
    # print(edges.astype(np.float16))

    edges = []
    for i in range(n):
        for j in range(i):
            if dists[i, j] <= max_dist:
                edges.append([i, j])
    edges = np.asarray(edges, dtype=np.int)

    edge_dict = {deg: edges}
    if np.any(edges):
        edge_positions = np.sum(positions[edges], axis=1) / 2
        edge_positions[..., 2] = np.max(positions[edges][..., 2], axis=1)
        if deg < max_deg:
            e, p = connectivity(edge_positions, max_dist * .8, deg + 1, max_deg)
            edge_dict.update(e)
            position_dict.update(p)
        else:
            position_dict[deg] = edge_positions
        return edge_dict, position_dict
    return {}, {}


def get_remap(positions, new_max: np.ndarray, padding: np.ndarray):
    max = np.max(positions, axis=0)
    min = np.min(positions, axis=0)
    def remap(positions: np.ndarray):
        positions = positions.copy()
        positions -= min
        positions *= (new_max - 2 * padding) / (max - min)
        positions += padding
        return positions.astype(np.int)
    return remap


def get_colormap(min: int, max: int):
    scale = 255. / (max - min + 1)
    def colormap(v):
        v -= min
        v *= scale
        v = int(v)
        color = [int(x * 255.) for x in plt.get_cmap('hsv')(v)]
        return tuple(color[:-1])
    return colormap


def draw_cut(atoms: ase.Atoms, max_dist=None, img_size=(800, 800), padding=20):
    positions = atoms.get_positions().astype(np.float)
    remap = get_remap(positions, new_max=np.asarray([img_size[0], img_size[1], 255]),
                      padding=np.asarray([padding, padding, 100]))

    img = np.zeros((img_size[0], img_size[1], 3), dtype=np.uint8)
    edge_dict, position_dict = connectivity(atoms.get_positions(), max_dist=max_dist)
    edge_sets = edge_dict.keys()
    if edge_sets:
        max_edge, min_edge = max(edge_sets), min(edge_sets)
        cmap = get_colormap(min_edge, max_edge)
        # TODO: refactor this hunk of code
        for edge_set in range(max_edge, min_edge - 1, -1):
            edges = edge_dict[edge_set]
            remapped = remap(position_dict[edge_set-1])[..., :-1]
            for i, edge in enumerate(edges):
                x0, y0 = remapped[edge[0]]
                x1, y1 = remapped[edge[1]]
                img = cv2.line(img, (x0, img_size[1] - y0), (x1, img_size[1] - y1), color=cmap(edge_set), thickness=edge_set)

                x, y = remap(position_dict[edge_set][i])[..., :-1]
                img = cv2.circle(img, center=(x, img_size[1] - y), radius=edge_set*2, color=cmap(edge_set+.1), thickness=-1)

    for i, p in enumerate(remap(positions)):
        x, y, d = p
        d = int(d)
        img = cv2.circle(img, center=(x, img_size[1] - y), radius=10, color=(d,d,d), thickness=-1)

    return img


if __name__ == '__main__':
    def na(x):
        pass

    cv2.namedWindow('img')
    cv2.createTrackbar('X', 'img', 0, 100, na)

    size=(5,5,1)
    atoms = fcc111('Al', size=size)
    # atoms = bcc111('Fe', size=size)
    view(atoms)


    last_d = None
    while True:
        d = cv2.getTrackbarPos('X', 'img')
        if d != last_d:
            img = draw_cut(atoms, max_dist=.5+d/30.)
            last_d = d
        cv2.imshow('img', img)

        k = cv2.waitKey(1) & 0xFF
        if k == ord('q'):
            cv2.destroyAllWindows()
        elif k == ord('s'):
            cv2.imwrite('connectivity.png', img)
            cv2.destroyAllWindows()
