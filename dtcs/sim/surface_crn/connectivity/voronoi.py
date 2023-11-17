import numpy as np
from scipy.spatial import Voronoi
import copy
from ase import Atoms


def produce_voronoi(points):
    return Voronoi(points)


def voronoi_infinite_regions(vor):
    """
    Produce a version of vor.ridge_points, where each infinite point is given
    an index, attached to the end of a modified version of vor.vertices.

    :param vor: a Voronoi object returned by Scipy
    :return: modified version of vor.ridge_points, and vor.vertices.
    """
    center = vor.points.mean(axis=0)
    ptp_bound = vor.points.ptp(axis=0)

    finite_segments = []
    infinite_segments = []
    far_points = []

    updated_vertices = vor.vertices.copy()
    # TODO: this should be ridge_vertices.
    updated_ridge_points = copy.deepcopy(vor.ridge_points)
    for pointidx, simplex in zip(vor.ridge_points, vor.ridge_vertices):
        simplex = np.asarray(simplex)
        new_simplex = np.copy(simplex)
        if np.all(simplex >= 0):
            finite_segments.append(simplex)
        else:
            i = simplex[simplex >= 0][0]  # finite end Voronoi vertex
            infinite_index = [simplex < 0][0]

            t = vor.points[pointidx[1]] - vor.points[pointidx[0]]  # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[pointidx].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            if (vor.furthest_site):
                direction = -direction
            far_point = vor.vertices[i] + direction * ptp_bound.max()

            infinite_segments.append([i, len(far_points)])
            far_points.append(far_point)

            new_simplex = new_simplex.tolist()
            new_simplex[infinite_index] = len(updated_vertices)

            # TODO:
            updated_vertices.append(far_point)
        updated_ridge_points.append(new_simplex)
    return updated_ridge_points, updated_vertices


def voronoi_plot_2d(vor, ax=None, **kw):
    """
    A customized version of voronoi plotting based on Scipy.

    Plot the given Voronoi diagram in 2-D
    Parameters
    ----------
    vor : scipy.spatial.Voronoi instance
        Diagram to plot
    ax : matplotlib.axes.Axes instance, optional
        Axes to plot on
    show_points: bool, optional
        Add the Voronoi points to the plot.
    show_vertices : bool, optional
        Add the Voronoi vertices to the plot.
    line_colors : string, optional
        Specifies the line color for polygon boundaries
    line_width : float, optional
        Specifies the line width for polygon boundaries
    line_alpha: float, optional
        Specifies the line alpha for polygon boundaries
    point_size: float, optional
        Specifies the size of points
    Returns
    -------
    fig : matplotlib.figure.Figure instance
        Figure for the plot
    See Also
    --------
    Voronoi
    Notes
    -----
    Requires Matplotlib.
    Examples
    --------
    Set of point:
    >>> import matplotlib.pyplot as plt
    >>> points = np.random.rand(10,2) #random
    Voronoi diagram of the points:
    >>> from scipy.spatial import Voronoi, voronoi_plot_2d
    >>> vor = Voronoi(points)
    using `voronoi_plot_2d` for visualisation:
    >>> fig = voronoi_plot_2d(vor)
    using `voronoi_plot_2d` for visualisation with enhancements:
    >>> fig = voronoi_plot_2d(vor, show_vertices=False, line_colors='orange',
    ...                 line_width=2, line_alpha=0.6, point_size=2)
    >>> plt.show()
    """
    from matplotlib.collections import LineCollection

    if vor.points.shape[1] != 2:
        raise ValueError("Voronoi diagram is not 2-D")

    line_colors = kw.get('line_colors', 'k')
    line_width = kw.get('line_width', 1.0)
    line_alpha = kw.get('line_alpha', 1.0)

    center = vor.points.mean(axis=0)
    ptp_bound = vor.points.ptp(axis=0)

    finite_segments = []
    infinite_segments = []
    for pointidx, simplex in zip(vor.ridge_points, vor.ridge_vertices):
        simplex = np.asarray(simplex)
        if np.all(simplex >= 0):
            finite_segments.append(vor.vertices[simplex])
        else:
            i = simplex[simplex >= 0][0]  # finite end Voronoi vertex

            t = vor.points[pointidx[1]] - vor.points[pointidx[0]]  # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[pointidx].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            if (vor.furthest_site):
                direction = -direction
            far_point = vor.vertices[i] + direction * ptp_bound.max()

            # ax.plot(far_point[0], far_point[1], 'o', markersize=5 * point_size, color="orange")
            infinite_segments.append([vor.vertices[i], far_point])

    ax.add_collection(LineCollection(finite_segments,
                                     colors=line_colors,
                                     lw=line_width,
                                     alpha=line_alpha,
                                     linestyle='solid'))
    ax.add_collection(LineCollection(infinite_segments,
                                     colors=line_colors,
                                     lw=line_width,
                                     alpha=line_alpha,
                                     linestyle='solid'))

    _adjust_bounds(ax, vor.points)

    return ax.figure


def _adjust_bounds(ax, points):
    margin = 0.1 * points.ptp(axis=0)
    xy_min = points.min(axis=0) - margin
    xy_max = points.max(axis=0) + margin
    ax.set_xlim(xy_min[0], xy_max[0])
    ax.set_ylim(xy_min[1], xy_max[1])


def voronoi_finite_polygons_2d(vor, radius=None):
    """
    Based on https://stackoverflow.com/questions/20515554/colorize-voronoi-diagram, with
    major modifications.

    Reconstruct infinite voronoi regions in a 2D diagram to finite
    regions.

    Parameters
    ----------
    vor : Voronoi
        Input diagram
    radius : float, optional
        Distance to 'points at infinity'.

    Returns
    -------
    regions : list of tuples
        Indices of vertices in each revised Voronoi regions.
    vertices : list of tuples
        Coordinates for revised Voronoi vertices. Same as coordinates
        of input vertices, with 'points at infinity' appended to the
        end.

    """
    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")

    new_regions = [[]]
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)
    if radius is None:
        radius = vor.points.ptp().max()
        ptp_bound = vor.points.ptp(axis=0)

    # Construct a map containing all ridges for a given point
    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        # Deal with the triangles with only 1 finite point
        # In this case, (v1, -1) may pair with two other ridge points.
        if p1 in all_ridges and (v1, v2) in all_ridges[p1]:
            all_ridges.setdefault(p1, {})[(v2, v1)] = p2
        else:
            all_ridges.setdefault(p1, {})[(v1, v2)] = p2

        if p2 in all_ridges and (v1, v2) in all_ridges[p2]:
            all_ridges.setdefault(p2, {})[(v2, v1)] = p1
        else:
            all_ridges.setdefault(p2, {})[(v1, v2)] = p1

    # Reconstruct infinite regions
    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]

        if all(v >= 0 for v in vertices):
            # finite region
            new_regions.append(vertices)
            continue

        # reconstruct a non-finite region
        ridges = all_ridges[p1]
        new_region_vertices = []

        for i in range(len(vertices)):
            if vertices[i] >= 0:
                new_region_vertices.append(vertices[i])
                continue

            # Make the infinite region into a triangle, one of whose vertex would be its head.
            if len(vertices) == 2:
                next_vertex_index = prev_vertex_index = 0 if i == 1 else 1
            else:
                next_vertex_index = 0 if i + 1 == len(vertices) else i + 1
                prev_vertex_index = len(vertices) - 1 if i == 0 else i - 1

            for v1_index, v2_index in [[prev_vertex_index, i], [i, next_vertex_index]]:
                v1, v2 = vertices[v1_index], vertices[v2_index]
                if v1 >= 0 and v2 >= 0:
                    continue

                # print(len(vertices))
                # print(vertices)
                # This gets helps on the issue of two edges of an infinity triangle sharing names.
                if (v1, v2) in ridges:
                    p2 = ridges[(v1, v2)]
                    # print(f"p2, normal choice ridges[({v1}, {v2})]", p2)
                else:
                    p2 = ridges[(v2, v1)]
                    # print(f"p2, choose other ridges[({v2}, {v1})]", p2)

                if v1 == -1:
                    v_near = v2
                else:
                    v_near = v1

                # Compute the missing endpoint of an infinite ridge
                t = vor.points[p2] - vor.points[p1]     # tangent
                t /= np.linalg.norm(t)
                n = np.array([-t[1], t[0]])  # normal

                midpoint = vor.points[[p1, p2]].mean(axis=0)
                direction = np.sign(np.dot(midpoint - center, n)) * n

                if (vor.furthest_site):
                    direction = -direction
                far_point = vor.vertices[v_near] + direction * ptp_bound.max()

                new_region_vertices.append(len(new_vertices))
                new_vertices.append(far_point)
        new_regions.append(new_region_vertices)
    return new_regions, np.asarray(new_vertices)


def fold_numbers(vor):
    """
    Return the fold numbers of the intersection nodes.
    """

    counts = [0 for _ in vor.vertices]

    for v1_index, v2_index in vor.ridge_vertices:
        if v1_index >= 0:
            counts[v1_index] += 1
        if v2_index >= 0:
            counts[v2_index] += 1

    #  TODO: account for edge nodes
    return [2, max(np.unique(counts))]

def voronoi_include_edge(points):
    """
    Produce a Voronoi object, where the regions at the edge are presumed to be same as the regions in the center.
    :param points:
    :return:
    """
    super_cell_positions = (Atoms("top", points) * (3, 3, 3)).get_positions()
    vor = Voronoi(super_cell_positions)
    included_points = [i for i in len(vor.points) if vor.points[i] in points]

    new_ridges = []
    new_ridge_points = []
    for i, vertices in enumerate(vor.ridge_vertices):
        vertex_1, vertex_2 = vertices
        if vertex_1 in included_points or vertex_2 in included_points:
            new_ridges.append(vertices)
            new_ridge_points.append(vor.ridge_points[i])

    # if a vertex is not adjacent to any point, don't include it.
    new_regions = []
    for i, region_index in enumerate(vor.point_region):
        if i in included_points:
            new_regions.append(vor.regions[region_index])


def create_supercell(points, cell=None, supercell_factor=(2, 2, 1)):
    points = points.tolist()
    points = [tuple(row + [0]) for row in points]

    atoms_positive = Atoms("O" * len(points), positions=points, cell=cell)
    pos_positions = (atoms_positive * (2, 2, 1)).get_positions()
    atoms_neg = Atoms("O" * len(pos_positions), positions=pos_positions, cell=cell*(-1))

    positive_points = [tuple(row[:-1]) for row in (atoms_positive * supercell_factor).get_positions()]
    negative_points = [tuple(row[:-1]) for row in (atoms_neg * supercell_factor).get_positions()]
    new_points = np.unique(np.array(list(set(tuple(row) for row in list(set(negative_points + positive_points))))),
                           axis=0)
    new_points = np.array(list(set(tuple(round(r, 3) for r in row) for row in new_points)))
    return new_points


class VoronoiGraph:
    def __init__(self, points):
        """
        A 2-d array of points.
        :param points:
        """
        vor = Voronoi(points)

        # Each dictionary attribute shall have immutable keys.
        self.points_dict = {i: p for i, p in enumerate(vor.points)}

        self.points = VoronoiGraph.index_dict_to_arrays(self.points_dict)

        self.vertices_dict = {i: v for i, v in enumerate(vor.vertices)}
        self.vertices = VoronoiGraph.index_dict_to_arrays(self.vertices_dict)

        self.regions_dict = {i: r for i, r in enumerate(vor.regions)}
        self.regions = VoronoiGraph.index_dict_to_arrays(self.regions_dict)
        self.point_region_dict = {i: region_index for i, region_index in enumerate(vor.point_region)}
        self.point_region = VoronoiGraph.index_dict_to_arrays(self.point_region_dict)

        # A dictionary of the points connected to each vertex.
        self.vertex_points_dict = self.compute_vertex_points()

        self.ridge_points = vor.ridge_points
        self.ridge_vertices = vor.ridge_vertices
        self.furthest_site = vor.furthest_site

    @staticmethod
    def index_dict_to_arrays(index_dict):
        """
        Update self.points according to points_dict.
        """
        arr = []
        for k in sorted(index_dict):
            arr.append(index_dict[k])
        import dtcs
        dtcs.debug = arr
        return np.array(arr, dtype=object)

    def compute_vertex_points(self):
        """
        Compute the list of points adjacent to each vertex.
        """
        new_vertex_points_dict = {}
        for point_index, region_index in self.point_region_dict.items():
            for vertex_index in self.regions_dict[region_index]:
                if vertex_index not in new_vertex_points_dict:
                    new_vertex_points_dict[vertex_index] = []
                new_vertex_points_dict[vertex_index].append(point_index)
        return new_vertex_points_dict

    def delete_vertex_point(self, vertex_index, point_index):
        """
        Delete a pair of vertex_index and its adjacent point_index from vertex_point_dict.
        If the last point_index of any vertex_index is deleted, delete the vertex_index from the dictionary.
        """
        # -1 won't count as a valid vertex index
        if vertex_index == -1:
            return
        point_indices = self.vertex_points_dict[vertex_index]
        self.vertex_points_dict[vertex_index] = [i for i in point_indices if i != point_index]
        if self.vertex_points_dict[vertex_index] == []:
            #  TODO: this should be popping out everything?
            self.vertex_points_dict.pop(vertex_index)
            self.vertices_dict.pop(vertex_index)

    def filter_points(self, points):
        """
        Delete points whose coordinates match any row in 2-d array points.
        """
        points = [tuple(round(r, 3) for r in row) for row in points]
        new_ridge_vertices = []
        new_ridge_points = []
        for i, vertices in enumerate(self.ridge_vertices):
            ridge_point_1, ridge_point_2 = self.ridge_points[i]
            if tuple(self.points_dict[ridge_point_1]) in points or tuple(self.points[ridge_point_2]) in points:
                new_ridge_vertices.append(vertices)

                # TODO: one of the points could be outside of any valid region.
                new_ridge_points.append(self.ridge_points[i])

        new_points_dict = {}
        for point_index, point in self.points_dict.items():
            if tuple(point) in points:
                new_points_dict[point_index] = point
            else:
                region_index = self.point_region_dict.pop(point_index)

                # Get rid of the vertex-point relationships based on this point
                if region_index in self.regions_dict:
                    popped_region = self.regions_dict.pop(region_index)
                    for vertex_index in popped_region:
                        self.delete_vertex_point(vertex_index, point_index)

        self.ridge_vertices = new_ridge_vertices
        self.ridge_points = new_ridge_points
        self.points_dict = new_points_dict
        self.reindex()

    def delete_points(self, points):
        """
        Delete points whose coordinates match any row in 2-d array points.
        """
        print("start deleting the points")
        points = [tuple(row) for row in points]
        new_ridge_vertices = []
        new_ridge_points = []
        for i, vertices in enumerate(self.ridge_vertices):
            ridge_point_1, ridge_point_2 = self.ridge_points[i]
            print(f'deleting vertex number {i}')
            if tuple(self.points_dict[ridge_point_1]) in points or tuple(self.points[ridge_point_2]) in points:
                new_ridge_vertices.append(vertices)
                new_ridge_points.append(self.ridge_points[i])

        print("finish calculating the ridge points")
        new_points_dict = {}
        for point_index, point in self.points_dict.items():
            if tuple(point) in points:
                region_index = self.point_region_dict.pop(point_index)

                popped_region = self.regions_dict.pop(region_index)
                for vertex_index in popped_region:
                    self.delete_vertex_point(vertex_index, point_index)
            else:
                new_points_dict[point_index] = point
        # self.reindex()
        print("points deletion finished")

    def reindex(self):
        """
        Give sequential indices to all points and vertices, intended for use after removals within the VoronoiGraph.
        """
        point_indexing_dict = {original_index: new_index for new_index,
                               original_index in enumerate(sorted(self.points_dict))}

        vertex_indexing_dict = {original_index: new_index for new_index,
                                original_index in enumerate(sorted(self.vertices_dict))}

        region_indexing_dict = {original_index: new_index for new_index,
                                original_index in enumerate(sorted(self.regions_dict))}

        for region in self.regions_dict.values():
            for i in range(len(region)):
                region[i] = vertex_indexing_dict[region[i]] if region[i] != -1 else -1

        for point_index in self.point_region_dict:
            if self.point_region_dict[point_index] in region_indexing_dict:
                self.point_region_dict[point_index] = region_indexing_dict[self.point_region_dict[point_index]]
            else:
                print(f"{self.point_region_dict[point_index]} not in region_indexing_dict")

        for i in range(len(self.ridge_points)):
            new_points = []
            for point_index in self.ridge_points[i]:
                if point_index in point_indexing_dict:
                    new_points.append(point_indexing_dict[point_index])
                else:
                    new_points.append(-1)
            self.ridge_points[i] = new_points

        for i in range(len(self.ridge_vertices)):
            new_vertices = []
            for vertex_index in self.ridge_vertices[i]:
                if vertex_index == -1:
                    new_vertices.append(-1)
                else:
                    new_vertices.append(vertex_indexing_dict[vertex_index])

            # TODO: fix ridge indexing error.
            self.ridge_vertices[i] = new_vertices

        new_vertex_points_dict = {}
        for vertex_index, points in self.vertex_points_dict.items():
            if vertex_index != -1:
                new_points = []
                for point_index in points:
                    if point_index in point_indexing_dict:
                        new_points.append(point_index)

                new_vertex_points_dict[vertex_indexing_dict[vertex_index]] = new_points
        self.vertex_points_dict = new_vertex_points_dict

        new_points_dict = {}
        for new_point_index, old_point_index in enumerate(sorted(self.points_dict)):
            new_points_dict[new_point_index] = self.points_dict[old_point_index]
        self.points_dict = new_points_dict

        new_vertices_dict = {}
        for new_vertex_index, old_vertex_index in enumerate(sorted(self.vertices_dict)):
            new_vertices_dict[new_vertex_index] = self.vertices_dict[old_vertex_index]
        self.vertices_dict = new_vertices_dict

        new_regions_dict = {}
        for new_region_index, old_region_index in enumerate(sorted(self.regions_dict)):
            new_regions_dict[new_region_index] = self.regions_dict[old_region_index]
        self.regions_dict = new_regions_dict

        self.points = VoronoiGraph.index_dict_to_arrays(self.points_dict)
        self.vertices = VoronoiGraph.index_dict_to_arrays(self.vertices_dict)
        self.regions = VoronoiGraph.index_dict_to_arrays(self.regions_dict)
        self.point_region = VoronoiGraph.index_dict_to_arrays(self.point_region_dict)




















