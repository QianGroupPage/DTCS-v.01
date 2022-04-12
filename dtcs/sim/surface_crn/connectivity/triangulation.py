import numpy as np
import pandas as pd
from scipy.spatial import Delaunay
from ase.io.vasp import read_vasp
from ase.visualize import view


def show_triangulation(file_path="", points=None, project_down="z", title="Triangulation", display_zoom=2,
                       distort_factor=1.1, ignore_threhold=2):
    """
    :param file_path: a string for the path to a POSCAR file;
    :param points: a numpy array containing x, y, z locations for points;
    :param project_down: the direction from which you are looking into the surface
    :param title: title for the figure;
    :param display_zoom: zoom level for the figure, default to 2. Increase to enlarge the figure.
    :return:
    """
    # TODO: refactor the code below.
    if file_path:
        points = read_vasp(file_path).get_positions()
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
    else:
        raise Exception(f"{project_down} is not a possible value for parameter project_down. project down must be" +
                        " one of \'x\', \'y\', \'z\'.")

    R = np.array(((1, 0), (0, distort_factor)))
    R_inv = np.linalg.inv(R)

    if len(points) == 0:
        points = pd.read_csv(file_path).to_numpy()[:, 1:]
    points = de_dup_threshold(points, threshold=ignore_threhold)
    points = points[:, [i - 1 for i in indices]]
    points = de_dup(points)

    points = points.dot(R)
    tri = Delaunay(points)
    points = points.dot(R_inv)
    # plt.triplot(points[:, 0], points[:, 1], tri.simplices)
    # fig = plt.gcf()
    # w, h = fig.get_size_inches()
    # fig.set_size_inches(w * display_zoom, h * display_zoom)
    # plt.gca().set_aspect('equal', adjustable='box')
    #
    # plt.plot(points[:, 0], points[:, 1], 'o')
    #
    # plt.title(title)
    # plt.show()
    return tri, points


def de_dup_threshold(points, threshold=1, project_down="z"):
    """
    Drop all points except for all the points within a certain threshold from the top layer.
    :param points:
    :param threshold: in Angstroms; if an atom is further away from the top layer than this number,
    drop the atom.
    :return: all the points that remain.
    """
    if project_down == "z":
        index = 2
    elif project_down == "y":
        index = 1
    elif project_down == "x":
        index = 0

    heights = set(h for h in points[:, index])
    return points[points[:, index] > max(heights) - threshold]


def de_dup(arr):
    # https://stackoverflow.com/questions/31097247/remove-duplicate-rows-of-a-numpy-array
    new_array = [tuple(row) for row in arr]
    return np.array([[e for e in r] for r in set(new_array)])


def poscar_to_positions(file_path, supercell_dimensions=1):
    if isinstance(supercell_dimensions, tuple) or isinstance(supercell_dimensions, list):
        if len(supercell_dimensions) > 3:
            raise Exception()
        elif len(supercell_dimensions) == 2:
            supercell_dimensions = tuple(supercell_dimensions + supercell_dimensions[-1])
        elif len(supercell_dimensions) == 1:
            supercell_dimensions = tuple(supercell_dimensions * 3)
    else:
        supercell_dimensions = tuple(supercell_dimensions for _ in range(3))
    atoms = read_vasp(file_path) * supercell_dimensions

    # view(atoms)
    # TODO: plus atoms.cell into the atoms obbject used in voronoi/create_supercell
    return atoms.get_positions(), atoms


def grid_size(points):
    """
    :param points: a 2-D array of x, y locations
    :return: the number of unique x locations and the number of unique y locations.
    """
    return np.unique(points[:, 0]).size, np.unique(points[:, 1]).size
