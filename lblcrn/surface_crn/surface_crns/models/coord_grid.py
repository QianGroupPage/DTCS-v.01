from lblcrn.surface_crn.surface_crns.base.node import Node
from lblcrn.surface_crn.connectivity.triangulation import show_triangulation, poscar_to_positions
from lblcrn.surface_crn.connectivity.neighbors import voronoi_neighbors_dict
from lblcrn.surface_crn.connectivity.voronoi import voronoi_plot_2d, voronoi_finite_polygons_2d, \
    VoronoiGraph, create_supercell
from lblcrn.common.num_to_word import num2word
from collections import Counter
from ase import Atoms, Atom
from scipy.spatial import Voronoi
import numpy as np
import json
import math
import random


class CoordGrid(object):
    """
    Representation of a CRN on a 2D finite mesh grid, where we're given each node's exact coordinates.

    Only allows reactions between directly-adjacent locations, and every edge has equal weight.
    """
    def __init__(self, points=[], vor=None, distort_factor=1.1, ignore_threhold=1, supercell_dimensions=1):
        """
        :param points: a numpy array where each row is a coordinate for a top site.
        """
        # Provide three consistently equivalent ways to refer to the nodes.
        self.nodes = []
        self.nodes_by_site = {}
        self.id_to_node = {}

        # map id to node for all nodes, including disabled nodes.
        self.id_to_all_nodes = {}
        # A dictionary from the x coordinate of the row to the list of nodes that belong to the row, only for
        #  square grids.
        self.rows_dict = {}

        if vor is None:
            #  TODO: currently this does not use distort factor; so factor it out.
            tri, points = show_triangulation(points=points, distort_factor=distort_factor,
                                             ignore_threhold=ignore_threhold)
            # TODO: convert points into ASE cell, apply supercell, and then cut off redundant atoms, vertices,
            #  and ridges.
            # Below is the old coding
            # self._populate_grid_delaunay(neighbors_dict(tri, points))
            self.vor = Voronoi(points)
        else:
            self.vor = vor
            points = vor.points
        self._populate_grid_voronoi(points)
        self._number_grid()

        min_x, max_x = math.inf, -math.inf
        min_y, max_y = math.inf, -math.inf
        for node in self.nodes:
            node_x, node_y = node.position
            min_x = min(min_x, node_x)
            min_y = min(min_y, node_y)

            max_x = max(max_x, node_x)
            max_y = max(max_y, node_y)
        self.xrange = (min_x, max_x)
        self.yrange = (min_y, max_y)

        self.rows_dict = self.__compute_rows()

    @staticmethod
    def from_poscar(poscar_path, distort_factor=1.1, ignore_threhold=2, supercell_dimensions=1):
        """
        Transform a POSCAR file into a Coordinate Grid.
        :param distort_factor:
        :param ignore_threhold:
        :return:
        """
        # super_positions = poscar_to_positions(poscar_path, supercell_dimensions * 3)
        allowed_positions, atoms = poscar_to_positions(poscar_path, supercell_dimensions * 1)
        # print("allowed_positions from 'poscar_to_positions'", allowed_positions)

        tri, allowed_positions = show_triangulation(points=allowed_positions, distort_factor=distort_factor,
                                         ignore_threhold=ignore_threhold, project_down="z")
        # points_to_delete = [point for point in super_positions if point not in allowed_positions]

        vor = VoronoiGraph(create_supercell(allowed_positions, cell=atoms.cell))
        # print("allowed positions", [[round(c, 2) for c in p] for p in allowed_positions])
        vor.filter_points(allowed_positions)
        # print("final points", vor.points)
        new_grid = CoordGrid(vor=vor, distort_factor=distort_factor,
                         ignore_threhold=ignore_threhold, supercell_dimensions=supercell_dimensions)
        new_grid.vor = vor
        return new_grid

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
        neighbors_dict = voronoi_neighbors_dict(self.vor)

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
        self.nodes_by_site["Bridge"] = []
        for i, vertex_indices in enumerate(self.vor.ridge_vertices):
            v1_index, v2_index = vertex_indices
            if v1_index == -1 or v2_index == -1:
                p1_index, p2_index = self.vor.ridge_points[i]
                loc = tuple(np.mean(self.vor.points[[p1_index, p2_index]], axis=0))
            else:
                loc = tuple(np.mean(self.vor.vertices[[v1_index, v2_index]], axis=0))
            new_node = Node(position=loc, state="Bridge")
            nodes[loc] = new_node
            self.nodes_by_site["Bridge"].append(new_node)
            self.nodes.append(new_node)

        self.nodes_by_site["Intersection"] = []
        for loc in self.vor.vertices:
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
                loc = tuple(c for c in self.vor.vertices[intersection_index])
                current_node.neighbors.append((nodes[loc], 1))

            for bridge_index in individual_neighbors_dict["Bridge"]:
                bridge_neighbor_node = self.nodes_by_site["Bridge"][bridge_index]
                current_node.neighbors.append((bridge_neighbor_node, 1))

        for i, individual_neighbors_dict in enumerate(neighbors_dict["Intersection"]):
            loc = tuple(c for c in self.vor.vertices[i])
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

    def __disable_node(self, node_id):
        """
        Temporarily prevent a node from showing in the simulations.
        """
        self.id_to_node.pop(node_id)
        self.nodes = list(self.id_to_node.values())

    def __reenable_node(self, node_id):
        """
        Re-enable a disabled node.
        """
        node = self.id_to_all_nodes[node_id]
        self.id_to_node[node_id] = node
        self.nodes.append(node)

    def __compute_rows(self):
        """
        Compute the lists of nodes with the same x-axis values. Assemble the lists into a dictionary
        indexed by x-values.
        """
        # Map a tuple of coordinates, where each number is rounded to 2 digits after decimal, to
        # the corresponding node object.
        position_tuple_to_node = {}
        for node in self.nodes_by_site["Top"]:
            position_tuple_to_node[tuple(round(c, 2) for c in node.position)] = node
        x_positions = set(k[0] for k in position_tuple_to_node.keys())
        rows = {x_position: [] for x_position in x_positions}
        for x, y in position_tuple_to_node:
            rows[x].append(position_tuple_to_node[x, y])
        return rows

    def disable_alternate_rows(self):
        # TODO: only top sites should be counted, then associated other sites get filtered out.
        sorted_x_vals = sorted(self.rows_dict.keys())
        affected_neighbors = {}
        for i in range(0, len(sorted_x_vals), 2):
            for node in self.rows_dict[sorted_x_vals[i]]:
                self.__disable_node(node.node_id)

                for neighbor_node in node.neighbors:
                    if neighbor_node in affected_neighbors:
                        affected_neighbors[neighbor_node] = 1
                    affected_neighbors[neighbor_node] += 1

                    # disabled_node should not be in node.neighbors.

        # 1-degree recursive disabling of neighboring nodes, this needs to be recursive
        for neighbor_node, affected_times in affected_neighbors.items():
            if len(neighbor_node.neighbors) == affected_times:
                self.__disable_node(neighbor_node.node_id)
        return

    def enable_alternate_rows(self):
        sorted_x_vals = sorted(self.rows_dict.keys())
        affected_neighbors = {}
        for i in range(0, len(sorted_x_vals), 2):
            for node in self.rows_dict[sorted_x_vals[i]]:
                self.__reenable_node(node.node_id)

                for neighbor_node in node.neighbors:
                    if neighbor_node in affected_neighbors:
                        affected_neighbors[neighbor_node] = 1
                    affected_neighbors[neighbor_node] += 1

            #     TODO: the re-enabled nodes must be put back into its neighbor's list of nodes.
            # 1-degree recursive disabling of neighboring nodes, this needs to be recursive
            for neighbor_node, affected_times in affected_neighbors.items():
                if len(neighbor_node.neighbors) == affected_times:
                    self.__reenable_node(neighbor_node.node_id)
        return

    def _number_grid(self):
        """
        Give each node in the grid a unique identifying number based on its location.

        This gives the initial numbering, before any node is disabled.
        :return:
        """
        site_names = ["Top", "Bridge", "Intersection"]
        site_abbreviations = {s: s[0] for s in site_names}
        self.id_to_node = {}
        for site_name in site_names:
            # TODO: this is not clean
            nodes = [n for n in self.nodes_by_site[site_name] if n is not None]
            nodes = sorted(nodes, key=lambda n: (n.position[1], n.position[0]))

            for i, n in enumerate(nodes):
                n.node_id = site_abbreviations[site_name] + str(i)
                self.id_to_node[n.node_id] = n
        self.id_to_all_nodes = self.id_to_node.copy()

    def _get_node(self, node_id):
        if not self.id_to_node:
            self._number_grid()
        return self.id_to_node[node_id]

    def clear_timestamps(self):
        """
        Set the timestamps of all nodes in the grid to 0.
        """
        for node in self:
            node.timestamp = 0
        return

    def get_global_state(self):
        """
        Get the global state of nodes as a dictionary from node_id to strings.
        """
        state_dict = {}
        for node in self:
            state_dict[node.node_id] = node.state
        return state_dict

    def set_global_state(self, state_dict):
        """
        Set the states of nodes using a dictionary from node ID to node object.
        Also resets timestamps.
        """
        for node_id, node_state in state_dict.items():
            self._get_node(node_id).state = node_state
        self.clear_timestamps()

    def set_initial_concentrations(self, initial_concentrations, random_seed=30):
        """
        Set the initial concentrations of all species following the dictionary
        of initial_concentrations.
        """
        random.seed(random_seed)
        for fold_name, d in initial_concentrations.items():
            nodes = self.nodes_by_fold_name[fold_name]
            species = []
            for species_name, species_concentration in d.items():
                species.extend([species_name] * species_concentration)
            indices = random.sample(range(len(nodes)), len(species))
            chosen = [nodes[i] for i in indices]
            for i, n in enumerate(chosen):
                n.state = species[i]

    def set_default_species(self, default_species_dict):
        """
        Set the default states of the grid according to default_species_dict.
        """
        for site_name, species_name in default_species_dict.items():
            for node in self.nodes_by_fold_name[site_name]:
                node.state = species_name

    def print_connectivity_dictionary(self, node=None):
        """
        Produce and then pretty print the connectivity dictionary
        """
        if node is not None:
            d = self.connectivity_dictionary
            res = {}
            for type, type_d in d.items():
                res[type] = {}
                for node_id, neighbors in type_d.items():
                    if node in neighbors or node_id == node:
                        res[type][node_id] = neighbors.copy()
            print(json.dumps(res, sort_keys=True, indent=4))
        else:
            print(json.dumps(self.connectivity_dictionary, sort_keys=True, indent=4))

    def voronoi_pic(self, ax=None, return_fig=False, show_node_number=False, color_index=None, set_fig_size=True,
                    fig_type="voronoi"):

        """
        Produce a Voronoi figure for the current graph.
        """
        # TODO: save the polygons, the ridge lines, and the atoms, for fast drawing on a PyGame screen.
        vor = self.vor

        if not ax:
            import matplotlib.pyplot as plt
            ax = plt.gca()
        #     TODO: this may be screwing up the storage

        if fig_type == "colored_dots":
            fig = plt.gcf()
            # Colored-dots graph
            for node in self:
                if node is None:
                    continue
                else:
                    if color_index is None:
                        ax.plot(node.position[0], node.position[1], 'o', markersize=1, color="green")
                    else:
                        ax.plot(node.position[0], node.position[1], 'o', markersize=1, color=color_index[node.state])
        else:
            # ax.set_aspect('equal', adjustable='box')
            fig = voronoi_plot_2d(vor, ax=ax, point_size=2, show_points=False)

            for node in self.nodes_by_site["Intersection"]:
                if node:
                    # TODO: add node color
                    if color_index:
                        state = node.state
                        if state in self.fold_names_dict:
                            state = self.fold_names_dict[state]
                        species_color = color_index[state]
                    else:
                        species_color = "green"
                    ax.plot(node.position[0], node.position[1], 'o', markersize=8 * 2, color=species_color)

            # ax.plot(vor.vertices[:, 0], vor.vertices[:, 1], 'o', markersize=8 * 2, color="Green")

            #  TODO: fill the infinite regions
            regions, vertices = voronoi_finite_polygons_2d(self.vor)
            for point_index in range(len(vor.point_region)):
                region = regions[vor.point_region[point_index]]
                if -1 not in region:
                    polygon = [vertices[i] for i in region]
                    if color_index:
                        state = self.nodes_by_site["Top"][point_index].state
                        if state in self.fold_names_dict:
                            state = self.fold_names_dict[state]
                        species_color = color_index[state]
                    else:
                        species_color = "#CCCCCC"
                    if isinstance(species_color, list):
                        species_color = tuple(species_color)
                    # print(species_color)
                    ax.fill(*zip(*polygon), color=species_color, alpha=0.8)

                    # polygon = np.array(polygon)
                    # ax.plot(polygon[:, 0], polygon[:, 1], 'o', markersize=2)
            #
            # for i, vertex in enumerate(self.vor.vertices):
            #     ax.text(vertex[0], vertex[1], i)

            for node in self.nodes_by_site["Bridge"]:
                if node:
                    # TODO: add node color
                    if color_index:
                        state = node.state
                        if state in self.fold_names_dict:
                            state = self.fold_names_dict[state]
                        species_color = color_index[state]
                    else:
                        species_color = "red"
                    ax.plot(node.position[0], node.position[1], 'o', markersize=2 * 4, color=species_color)
                    # ax.text(node.position[0], node.position[1], 39)

        if show_node_number:
            for node in self.nodes:
                ax.text(node.position[0], node.position[1], node.node_id, ha='center', va='center')

        if set_fig_size:
            fig.set_figheight(8)
            fig.set_figwidth(12)

        # plt.axis('off')  # TODO: include this in the standard mode.
        # return fig
        # ax.plot(vertices[:, 0], vertices[:, 1], 'o', markersize=3 * 2, color="Red")
        ax.set_xlim(self.xrange[0] - 1, self.xrange[1] + 1)
        ax.set_ylim(self.yrange[0] - 1, self.yrange[1] + 1)
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
    def fold_names(self):
        """
        Build a list consisting of all site names named after the site intersection numbers.

        For instance, ["top", "twofold", "threefold", "fourfold"]
        :return:
        """
        res = []
        for key, nodes in self.nodes_by_site.items():
            if key == "Top":
                continue
            elif nodes:
                res.append(max([len(node.neighbors) for node in nodes]))
        return ["top"] + [f"{num2word(num)}fold" for num in res if num > 1]

    @property
    def fold_names_dict(self):
        d = {}
        for k, v in self.nodes_by_site.items():
            if k == "Top":
                d["Top"] = "top"
            elif k == "Bridge":
                d["Bridge"] = "twofold"
            else:
                intersection_name = [n for n in self.fold_names if n not in ("top", "twofold")][0]
                d["Intersection"] = intersection_name
        return d

    @property
    def nodes_by_fold_name(self):
        fold_names = self.fold_names
        d = {}
        for k, v in self.nodes_by_site.items():
            if k == "Top":
                d["top"] = v
            elif k == "Bridge":
                d["twofold"] = v
            else:
                intersection_name = [n for n in fold_names if n not in ("top", "twofold")][0]
                d[intersection_name] = v
        return d

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

    @property
    def size_2d(self):
        """
        The best approximate width and height of this grid, as if it were a rectangle.
        """
        return (self.x_size, self.y_size)

    # TODO: following is used to compute size, verify it is correct
    @property
    def x_size(self):
        return len(set(node.position[0] for node in self.nodes))

    @property
    def y_size(self):
        return len(set(node.position[1] for node in self.nodes))

    def __iter__(self):
        return iter(self.nodes)
