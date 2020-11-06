"""
Suggest fitting parameters based on data in the RawExperiment class.
"""

import numpy as np
import pandas as pd


def _read_suggestions_table(table_path):
    """
    :param table_path: path to the table of suggested values;
    :return: dataframe of suggestions.
    """
    suggestions_df = pd.read_excel(io=table_path)
    suggestions_df['line_shape_params'] = suggestions_df['line_shape_params'].apply(_read_string_tuples)
    return suggestions_df


def _read_string_tuples(s):
    """
    :param s: a string representing a tuple, in the format of `(1.1, 2,  3 ) '.
    :return: a tuple.
    """
    if isinstance(s, str):
        s = s.strip()[1:-1]
        return tuple(map(float, [c.strip() for c in s.split(',')]))
    elif np.isnan(s):
        return ()


suggestions_df = \
    _read_suggestions_table('/Users/ye/Desktop/lbl-crn/lblcrn/_resources/metals.xlsx')


def process_species_name(species_name):
    """
    :param species_name: a name of a species. The first section shall be the atom/molecule name, the last section shall
    be the energy level. The middle section is the orbital information; only the first quantum number is allowed;
    :return: a list containing each section of the species name.
    """
    species_name = species_name.rstrip()
    allowed_separators = ["-", " "]

    for separator in allowed_separators:
        species_name = species_name.replace(separator, "_")

    name_sections = [name_section.strip() for name_section in species_name.split("_")]
    return name_sections


def suggest_fitting_params(species_name):
    """
    :param species_name: name of a species; the orbital must not have a spin number;
    :return: a list of dictionaries of fitting parameters, each dictionary corresponding to one peak.
    """
    element, orbital, _ = process_species_name(species_name)
    if element in [s for s in suggestions_df["element"].tolist()]:
        species_df = suggestions_df.loc[suggestions_df["element"] == element]

        if not species_df.empty:
            orbital_df = species_df[species_df.apply(lambda row: row["orbital"].startswith(orbital), axis=1)]

            if orbital_df.empty:
                print(f"Element {element} Orbital {orbital} is not in database.")
                return

            main_df = orbital_df.loc[orbital_df["orbital"] == orbital]
            if not main_df.empty:
                main_dict = orbital_df.iloc[0].to_dict()
                if len(orbital_df.index) == 1:
                    suggestions = [main_dict]
                # if there's another row, take it as the lower binding energy peak
                elif len(orbital_df.index) == 2:
                    low_be = species_df.loc[species_df["orbital"] != orbital].iloc[0].to_dict()

                    high_be = low_be.copy()

                    for key in ['fwhm', 'fwhm_error', 'line_shape', 'line_shape_params']:
                        if np.isnan(low_be[key]) or not low_be[key]:
                            low_be[key] = main_dict[key]
                        high_be[key] = main_dict[key]

                    high_be["be"] = low_be["be"] + main_df.iloc[0]["split"]
                    high_be["be_error"] = low_be["be_error"] + main_df.iloc[0]["split_error"]

                    suggestions = [low_be, high_be]
            else:
                if len(species_df.index) == 2:
                    suggestions = [species_df.iloc[i].to_dict() for i in range(2)]
            return pd.DataFrame.from_records(suggestions)
    else:
        print(f"Element {element} not in database.")







