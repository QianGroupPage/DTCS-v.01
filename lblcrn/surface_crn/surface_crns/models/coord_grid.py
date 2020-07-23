from lblcrn.surface_crn.surface_crns.base.node import Node
from lblcrn.surface_crn.connectivity.triangulation import show_triangulation
from lblcrn.surface_crn.connectivity.neighbors import neighbors_dict
from collections import Counter


class CoordGrid(object):
    """
    Representation of a CRN on a 2D finite mesh grid, where we're given each node's exact coordinates.

    Only allows reactions between directly-adjacent locations, and every edge has equal weight.
    """
    def __init__(self, points):
        """
        :param points: a numpy array where each row is a coordinate for a top site.
        """
        self.nodes = []
        tri, points = show_triangulation(points=points)
        self.populate_grid(neighbors_dict(tri, points))

    def populate_grid(self, neighbors_dict):
        """
        Populate the nodes list of the grid.

        :param neighbors_dict: a dictionary where each key is a default name, for instance, "Top";
        The value would be a dictionary where each key, a tuple of cartesian coordinates,
        corresponds to a list of of the key's neighbors, also a tuple or cartesian coordinates.
        For example: {"Top": {(1, 1): [(0,1), (1,0)]}, "3F": {(0.5, 0.5): [(1,1), (0, 1)]}
        :return:
        """
        nodes = {}
        for default_name in neighbors_dict.keys():
            for loc in neighbors_dict[default_name]:
                nodes[loc] = Node(position=loc, state=default_name)
                self.nodes.append(nodes[loc])

        for default_name in neighbors_dict.keys():
            for loc, neighbors in neighbors_dict[default_name].items():
                nodes[loc].neighbors = [(1, nodes[n_loc]) for n_loc in neighbors_dict]

    def species_count(self):
        """
        :return: count of each species in the grid
        """
        return Counter([str(n.state) for n in self])

    @property
    def num_nodes(self):
        return sum([1 for _ in self])

    def __iter__(self):
        return iter(self.nodes)
