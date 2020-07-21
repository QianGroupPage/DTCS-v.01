FLOAT_REGEX = "-?[0-9]+\.*[0-9]+"

INT_REGEX = "-?[1-9][0-9]*"

SPECIES_REGEX = "([a-zA-Z][0-9]*)+"

# The index (starting from 1) where the list of species is written at
# The following line would contain the quantities of each species
SPECIES_TYPE_LINE_INDEX = 6

# Each line following a line containing any keyword in this list represents
# the coordinates for 1 column.
START_OF_COORDS = ["Cartesian", "Direct"]
DIRECT_REGEX = "Direct"

import re, sys, os
import numpy as np
import pandas as pd


def is_species_list(line):
    """
    Return true if a line represents the list of species in the POSCAR file.
    """
    species_regex = re.compile(SPECIES_REGEX)
    return all([species_regex.match(w) for w in line.split()])

def direct_to_cartesian(scale_factor, lattice_vectors_transposed, l):
    """
    Lattice_vectors_transposed and l have to be Python arrays;
    Each row in l is the direct coordinates of an atom.
    """
    coord_weights = np.sum(np.array(lattice_vectors_transposed), axis=0)
    results = []
    for r in l:
        results.append([scale_factor * coord_weights[i] * r[i]for i in range(len(r))])
    return results

def is_direct(coordinate_flag):
    """
    Check if coordinate_flag, a string, indicates direct coordinate.
    """
    direct_regex = re.compile(DIRECT_REGEX)
    return direct_regex.fullmatch(coordinate_flag)



def transform_to_csv(filename):
    """
    Parse POSCAR for a list of coordinates of all atoms in the file,
    and save to a Python file with the same name.

    :param filename: name of the file to transform.
    :return: None
    """
    f = open(filename, "r")
    name = os.path.splitext(filename)[0]
    new_coordinates = open("%s_coords.csv" % name, "w")

    float_regex = re.compile(FLOAT_REGEX)
    int_regex = re.compile(INT_REGEX)

    # Construct a regex based on "(Cartesian | Direct)".
    start_of_coords_regex = re.compile("(%s)" % "|".join(map(re.escape, START_OF_COORDS)))
    is_coords, is_species_quantity, found_species_list, c = False, False, False, 0

    lattice_vectors_transposed = []
    scale_factor, coord_type = None, None
    l = []

    species = None

    for line in f:
        words = line.split()
        c += 1

        if c == 2:
            if len(words) != 1:
                raise  ValueError("The 2nd line in your Poscar file indicates the scale factor. The line has to "
                       + "contain only 1 number. \n The erroneous line is shown as below: \n"
                         + "{}\n".format(line))
            else:
                scale_factor = float(words[0])
        if 3 <= c <= 5:
            if len(words) != 3:
                raise  ValueError("The 2nd line in your Poscar file indicates the scale factor. The line has to "
                       + "contain exactly 3 numbers separated by whitespace. \n The erroneous line is shown as below: \n"
                         + "{}\n".format(line))
            else:
                v_t = [float(w) for w in words if float_regex.fullmatch(w)]
                lattice_vectors_transposed.append(v_t)

        if is_coords and len(words) != 0:
            l.append([float(w) for w in words if float_regex.fullmatch(w)])
        elif is_coords and len(words) == 0:
            # An empty line after a bunch of coordinates implies a different session
            is_coords = False
        elif start_of_coords_regex.match(line):
            is_coords = True
            coord_type = words[0]
        elif is_species_list(line) and not found_species_list:
            species = words
            is_species_quantity = True
            found_species_list = True
        elif is_species_quantity:
            is_species_quantity = False
            quantity = [int(w) for w in words if int_regex.fullmatch(w)]

    if is_direct(coord_type):
        l = direct_to_cartesian(scale_factor, lattice_vectors_transposed, l)

    # Take care of certain Poscar files where the header indicates only one species name
    # and where's no following line indicating the number of this species
    if len(species) == 1:
        quantity = [len(l)]

    if species is not None:
        j = 0
        for i in range(len(species)):
            for _ in range(quantity[i]):
                l[j].insert(0, species[i])
                j += 1

    results = pd.DataFrame(np.array(l))
    results.to_csv(new_coordinates)

    f.close()
    new_coordinates.close()


def df_to_poscar(filename, df, meta_data, coord_system="Direct"):
    """
    Save a Pandas dataframe as a POSCAR file.

    :param filename: the name for the new POSCAR file.
    :param df: a Pandas dataframe, where each row represents an atom.
    :param meta_data: the first 5 lines of a POSCAR file.
    :param coord_system: either "Direct" or "Cartesian"
    :return:
    """
    f = open(filename, "w")
    f.write(meta_data)


    # Separate different species by 3 spaces.
    names, counts = [], []

    prev, curr = None, None
    # If the following remains zero after going through the column of names,
    # we know that the input dataframe is empty.
    count = 0

    for name in df.iloc[:, 0]:
        prev, curr = curr, name

        if prev is not None and prev == curr:
            count = count + 1
        elif prev is not None:
            names.append(prev)
            counts.append(count)
            # Count the current one.
            count = 1
        else:
            # Count 1 if prev is none and current is not None.
            count = 1

    f.write("{}{}/n".format(" "*3, (" "*3).join(names)))
    f.write("{}{}/n".format(" "*4, (" "*5).join(counts)))

    f.write("{}/n".format(coord_system))

    # Each coordinate value has 2 spaces ahead of it, including the first coordinate value.
    for (index_label, row_series) in df.iloc[:, 1:].iterrows():
        f.write("{}{}/n".format(" " * 2, (" " * 2).join(row_series.values)))

    # Attach a blank line towards the end of the file.
    f.write("/n")
    f.close()

def grab_meta_data(filename):
    """
    Take the meta_data section of a POSCAR file.

    We define this section of a POSCAR file as the first 5 lines, the header line, the scale factor,
    and the next 3 lines for a unit cell.
    :return: a string representing the meta_data.
    """
    with open(filename, "r") as f:
        return "".join(f.readlines[:5])




if __name__ == "__main__":

  for filename in sys.argv[1:]:
    transform_to_csv(filename)
