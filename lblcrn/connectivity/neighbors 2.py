import numpy as np
from itertools import permutations
from lblcrn.surface_crn.connectivity.voronoi import voronoi_infinite_regions


def neighbors_dict(tri, points):
    """
    :param tri: a tri object created by scipy.spatial.Delaunay;
    :param points: a numpy array where each row is a location of a point;
    :return: a dictionary from a location (tuple) to a list of its neighbors.
    """
    top_sites = [tuple(c for c in point) for point in points]
    top_sites_neighbors = {}

    two_sites = []
    two_sites_neighbors = {}

    three_sites = []
    three_sites_neighbors = {}

    for i, simplice in enumerate(tri.simplices):
        three_site = triangle_centroid(points[simplice])
        three_sites.append(tuple(c for c in three_site))
        three_site_key = f"3F{i}"

        for j, site_index in enumerate(simplice):
            if site_index not in top_sites_neighbors:
                top_sites_neighbors[site_index] = []
            top_sites_neighbors[site_index].extend(simplice[:j])
            if j < len(simplice) - 1:
                top_sites_neighbors[site_index].extend(simplice[j + 1:])
            top_sites_neighbors[site_index].append(three_site_key)

            if three_site_key not in three_sites_neighbors:
                three_sites_neighbors[three_site_key] = []
            three_sites_neighbors[three_site_key].append(site_index)

        # Add two sites: not compeletely implemented.
        # for index_1, index_2 in permutations(simplice):
        #     two_site_key = f"2F{len(two_sites)}"

    res = {"Top": {}, "3F": {}}
    for key, val in top_sites_neighbors.items():
        res["Top"][top_sites[key]] = []
        for n_index in val:
            if str(n_index).startswith("3F"):
                res["Top"][top_sites[key]].append(three_sites[int(str(n_index)[2:])])
            else:
                res["Top"][top_sites[key]].append(top_sites[n_index])

    for key, val in three_sites_neighbors.items():
        res["3F"][three_sites[int(key[2:])]] = []
        for n_index in val:
            res["3F"][three_sites[int(key[2:])]].append(top_sites[n_index])
    return res


def triangle_centroid(loc_array):
    """
    :param loc_array: a numpy array where each row is the coordinates of a points
    :return: the geometrical centroid
    """
    return np.mean(loc_array, axis=0)


def voronoi_neighbors_dict(vor):
    top_neighbors_list = [{"Top": [], "Intersection": [], "Bridge": []} for _ in vor.points]

    intersection_neighbors_list = [{"Top": []} for _ in vor.vertices]
    bridge_neighbors_list = [{"Top": []} for _ in vor.ridge_points]

    regions = vor.regions
    for point_index, region_index in enumerate(vor.point_region):
        for vertex_index in regions[region_index]:
            if vertex_index != -1:
                top_neighbors_list[point_index]["Intersection"].append(vertex_index)
                intersection_neighbors_list[vertex_index]["Top"].append(point_index)

    for i, indices in enumerate(vor.ridge_points):
        n1_index, n2_index = indices

        # TODO: see if -1 any of n1_index and n2_index will be -1.
        top_neighbors_list[n1_index]["Top"].append(n2_index)
        top_neighbors_list[n2_index]["Top"].append(n1_index)

        # TODO: to support points that used to be infinite on the edge, use updated ridge_vertices else where
        #  one on one correspondence between bridge and ridge (a pair of ridge points and ridger vertices), of
        #  course still hold.
        top_neighbors_list[n1_index]["Bridge"].append(i)
        top_neighbors_list[n2_index]["Bridge"].append(i)
        bridge_neighbors_list[i]["Top"].extend(indices)
    return {"Top": top_neighbors_list, "Intersection": intersection_neighbors_list, "Bridge": bridge_neighbors_list}

