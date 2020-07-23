import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
from ase.io.vasp import read_vasp


def show_triangulation(file_path="", points=None, project_down="z", title="Triangulation", display_zoom=2):
    """
    :param file_path: a string for the path to a POSCAR file;
    :param points: a numpy array containing x, y, z locations for points;
    :param project_down: the direction from which you are looking into the surface
    :param title: title for the figure;
    :param display_zoom: zoom level for the figure, default to 2. Increase to enlarge the figure.
    :return:
    """
    if file_path:
        points = poscar_to_positions(file_path)
        if points == []:
            raise Exception(f"No valid points in file {file_path}")
    elif points is None:
        raise Exception("Please provide an input file either as a POSCAR file or " +
                        "as a numpy array of x, y, z locations")

    if project_down == "z":
        indices = range(1, 3)
    elif project_down == "y":
        indices = [1, 3]
    elif project_down == "x":
        indices = [2, 3]

    R = np.array(((1, 0), (0, 1.1)))
    R_inv = np.linalg.inv(R)

    if len(points) == 0:
        points = pd.read_csv(file_path).to_numpy()[:, 1:]
    points = points[:, [i - 1 for i in indices]]
    points = de_dup(points)

    points = points.dot(R)
    tri = Delaunay(points)
    points = points.dot(R_inv)
    plt.triplot(points[:, 0], points[:, 1], tri.simplices)
    fig = plt.gcf()
    w, h = fig.get_size_inches()
    fig.set_size_inches(w * display_zoom, h * display_zoom)
    plt.gca().set_aspect('equal', adjustable='box')

    plt.plot(points[:, 0], points[:, 1], 'o')

    plt.title(title)
    plt.show()
    return tri, points


def de_dup(arr):
    # https://stackoverflow.com/questions/31097247/remove-duplicate-rows-of-a-numpy-array
    new_array = [tuple(row) for row in arr]
    return np.array([[e for e in r] for r in set(new_array)])


def poscar_to_positions(file_path):
    cell = read_vasp(file_path)
    return cell.get_positions()
