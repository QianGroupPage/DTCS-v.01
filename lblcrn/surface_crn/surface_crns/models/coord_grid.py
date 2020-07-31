from lblcrn.surface_crn.surface_crns.base.node import Node
from lblcrn.surface_crn.connectivity.triangulation import show_triangulation, poscar_to_positions
from lblcrn.surface_crn.connectivity.neighbors import neighbors_dict, voronoi_neighbors_dict
from lblcrn.surface_crn.connectivity.voronoi import voronoi_plot_2d, voronoi_finite_polygons_2d
from collections import Counter
from ase import Atoms, Atom
from scipy.spatial import Voronoi
import numpy as np
import json


class CoordGrid(object):
    """
    Representation of a CRN on a 2D finite mesh grid, where we're given each node's exact coordinates.

    Only allows reactions between directly-adjacent locations, and every edge has equal weight.
    """
    def __init__(self, points, distort_factor=1.1, ignore_threhold=1):
        """
        :param points: a numpy array where each row is a coordinate for a top site.
        """
        self.nodes = []
        self.nodes_by_site = {}
        #  TODO: currently this does not use distort factor; so factor it out.
        tri, points = show_triangulation(points=points, distort_factor=distort_factor,
                                         ignore_threhold=ignore_threhold)
        # Below is the old coding
        # self._populate_grid_delaunay(neighbors_dict(tri, points))
        self.vor = None
        self._populate_grid_voronoi(points)
        self._number_grid()

    @staticmethod
    def from_poscar(poscar_path, distort_factor=1.1, ignore_threhold=1):
        """
        Transform a POSCAR file into a Coordinate Grid.
        :param distort_factor:
        :param ignore_threhold:
        :return:
        """
        return CoordGrid(poscar_to_positions(poscar_path), distort_factor=distort_factor,
                         ignore_threhold=ignore_threhold)

    def _populate_grid_delaunay(self, neighbors_dict):
        """
        Populate the nodes list of the grid using Delaunay Triangulation.

        :param neighbors_dict: a dictionary where each key is a default name, for instance, "Top";
        The value would be a dictionary where each key, a tuple of cartesian coordinates,
        corresponds to a list of of the key's neighbors, also a tuple or cartesian coordinates.
        For example: {"Top": {(1, 1): [(0,1), (1,0)]}, "3F": {(0.5, 0.5): [(1,1), (0, 1)]}
        """
        nodes = {}
        for default_name in neighbors_dict.keys():
            self.nodes_by_site[default_name] = []
            for loc in neighbors_dict[default_name]:
                nodes[loc] = Node(position=loc, state=default_name)
                self.nodes.append(nodes[loc])
                self.nodes_by_site[default_name].append(nodes[loc])

        for default_name, d in neighbors_dict.items():
            for loc, neighbors in d.items():
                nodes[loc].neighbors = [(1, nodes[n_loc]) for n_loc in neighbors]

    def _populate_grid_voronoi(self, points):
        """
        Populate the grid based on a set of coordinates in points.

        :param points:
        """
        vor = Voronoi(points)
        self.vor = vor
        neighbors_dict = voronoi_neighbors_dict(vor)

        nodes = {}
        #  TODO: currently self.nodes_by_site["Top"][i] corresponds to points[i] and vor.points[i]
        #  TODO: the above is too implicit.
        self.nodes_by_site["Top"] = []
        for loc in points:
            loc = tuple(c for c in loc)
            new_node = Node(position=loc, state="Top")
            nodes[loc] = new_node

            # Nodes_by_site["top"] would be a list with same indexing as points
            self.nodes_by_site["Top"].append(new_node)
            self.nodes.append(new_node)

        # TODO: bridge sites.
        # vor.ridge_points  # There should be a 2-site between each of these.
        self.nodes_by_site["Bridge"] = []
        for i, vertex_indices in enumerate(self.vor.ridge_vertices):
            v1_index, v2_index = vertex_indices
            if v1_index == -1 or v2_index == -1:
                p1_index, p2_index = self.vor.ridge_points[i]
                # TODO: make sure this confusing use of "None" doesn't show up elsewhere.
                # self.nodes_by_site["Bridge"].append(None)
                # # self.nodes.append(None)
                loc = tuple(np.mean(self.vor.points[[p1_index, p2_index]], axis=0))
            else:
                loc = tuple(np.mean(self.vor.vertices[[v1_index, v2_index]], axis=0))
            new_node = Node(position=loc, state="Bridge")
            nodes[loc] = new_node
            self.nodes_by_site["Bridge"].append(new_node)
            self.nodes.append(new_node)

        #  TODO: completely fix connectivity
        #  Add the bridge site on a boundary line going to infinity
        # for i, point_indices in enumerate(self.vor.ridge_points):
        #     p1_index, p2_index = point_indices
        #     if any([vertex_index == -1 for vertex_index in self.vor.ridge_vertices[i]]):
        #         loc = tuple(np.mean(self.vor.points[[p1_index, p2_index]], axis=0))
        #         new_node = Node(position=loc, state="Bridge")
        #         nodes[loc] = new_node
        #         self.nodes_by_site["Bridge"].append(new_node)
        #         self.nodes.append(new_node)

        self.nodes_by_site["Intersection"] = []
        for loc in vor.vertices:
            loc = tuple(c for c in loc)
            new_node = Node(position=loc, state="Intersection")
            nodes[loc] = new_node
            self.nodes_by_site["Intersection"].append(new_node)
            self.nodes.append(new_node)

        # TODO: Make new nodes' neighbors
        for i, individual_neighbors_dict in enumerate(neighbors_dict["Top"]):
            loc = tuple(c for c in points[i])
            current_node = nodes[loc]

            for top_index in individual_neighbors_dict["Top"]:
                loc = tuple(c for c in points[top_index])
                current_node.neighbors.append((nodes[loc], 1))

            for intersection_index in individual_neighbors_dict["Intersection"]:
                loc = tuple(c for c in vor.vertices[intersection_index])
                current_node.neighbors.append((nodes[loc], 1))

            for bridge_index in individual_neighbors_dict["Bridge"]:
                bridge_neighbor_node = self.nodes_by_site["Bridge"][bridge_index]
                current_node.neighbors.append((bridge_neighbor_node, 1))

        for i, individual_neighbors_dict in enumerate(neighbors_dict["Intersection"]):
            loc = tuple(c for c in vor.vertices[i])
            current_node = nodes[loc]

            for top_index in individual_neighbors_dict["Top"]:
                loc = tuple(c for c in points[top_index])
                current_node.neighbors.append((nodes[loc], 1))

        for bridge_index, individual_neighbors_dict in enumerate(neighbors_dict["Bridge"]):
            current_node = self.nodes_by_site["Bridge"][bridge_index]
            # Handle the case where a node with an infinite neighbor is None.
            if not current_node:
                continue

            for top_index in individual_neighbors_dict["Top"]:
                current_node.neighbors.append((self.nodes_by_site["Top"][top_index], 1))

    def to_atoms(self):
        """
        Output an ASE atoms object representing the states, and the node locations of the current grid.
        """
        # map sites to "default" atoms for the purposes of using ASE's built-in information on specific atom types.
        site_to_atom = {"Top": "Ag", "Intersection": "Cl", "Bridge": "O"}
        # TODO: test and handle the case when n.state is a molecule.
        return Atoms([Atom(site_to_atom[n.state] if n.state in site_to_atom else n.state,
                     [c for c in n.position] + [0]) for n in self])

    def species_count(self):
        """
        :return: count of each species in the grid
        """
        return Counter([str(n.state) for n in self])

    def _number_grid(self):
        """
        Give each node in the grid a unique identifying number based on its location.
        :return:
        """
        site_names = ["Top", "Bridge", "Intersection"]
        site_abbreviations = {s: s[0] for s in site_names}
        for site_name in site_names:
            # TODO: this is not clean
            nodes = [n for n in self.nodes_by_site[site_name] if n is not None]
            nodes = sorted(nodes, key=lambda n: (n.position[1], n.position[0]))

            for i, n in enumerate(nodes):
                n.node_id = site_abbreviations[site_name] + str(i)

    def print_connectivity_dictionary(self):
        """
        Produce and then pretty print the connectivity dictionary
        """
        print(json.dumps(self.connectivity_dictionary, sort_keys=True, indent=4))

    def voronoi_pic(self, return_fig=False, show_node_number=False):
        """
        Produce a Voronoi figure for the current graph.
        """
        import matplotlib.pyplot as plt
        vor = self.vor

        ax = plt.gca()
        ax.set_aspect('equal', adjustable='box')
        fig = voronoi_plot_2d(vor, ax=ax, point_size=2, show_points=False)

        ax.plot(vor.vertices[:, 0], vor.vertices[:, 1], 'o', markersize=8 * 2, color="Green")

        #  TODO: fill the infinite regions
        regions, vertices = voronoi_finite_polygons_2d(self.vor)
        for point_index in range(len(vor.point_region)):
            region = regions[vor.point_region[point_index]]
            if -1 not in region:
                polygon = [vertices[i] for i in region]
                ax.fill(*zip(*polygon), color="#CCCCCC", alpha=0.8)

                # polygon = np.array(polygon)
                # ax.plot(polygon[:, 0], polygon[:, 1], 'o', markersize=2)
        #
        # for i, vertex in enumerate(self.vor.vertices):
        #     ax.text(vertex[0], vertex[1], i)

        for node in self.nodes_by_site["Bridge"]:
            if node:
                ax.plot(node.position[0], node.position[1], 'o', markersize=2 * 4, color="red")
                # ax.text(node.position[0], node.position[1], 39)

        if show_node_number:
            for node in self.nodes:
                ax.text(node.position[0], node.position[1], node.node_id, ha='center', va='center')

        fig.set_figheight(8)
        fig.set_figwidth(12)

        # plt.axis('off')  # TODO: include this in the standard mode.
        # return fig
        # ax.plot(vertices[:, 0], vertices[:, 1], 'o', markersize=3 * 2, color="Red")
        # ax.set_xlim(-10, 20)
        # ax.set_ylim(-15, 30)
        if return_fig:
            return fig

    @property
    def num_nodes(self):
        return sum([1 for _ in self])

    @property
    def top_nodes(self):
        return self.nodes_by_site["Top"].copy()

    @property
    def connectivity_dictionary(self):
        """
        Compute a dictionary representing the connectivity graph.
        """
        d = {}
        for site_name, nodes in self.nodes_by_site.items():
            d[site_name] = {}
            for node in nodes:
                if node is not None:
                    # TODO: weight should come second
                    d[site_name][node.node_id] = sorted([neighbor_node.node_id for neighbor_node, _ in node.neighbors
                                                         if neighbor_node is not None])
        return d

    def __iter__(self):
        return iter(self.nodes)
