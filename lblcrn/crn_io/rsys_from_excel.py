from lblcrn.crn_sym import SpeciesManager, Rxn, RevRxn, Conc, RxnSystem
from lblcrn.model_input.state_variable import T, P
from lblcrn.tables import read_line_delimited_excel
import numpy as np
import pandas as pd


def build_bulk_crn(species, rxns, initial_concentration):
    sm = SpeciesManager()
    species_name_to_sympy = {s: (sm.sp(s) if s != "*" else sm.sp("S")) for s in species['Species']}
    rxn_sys_input = []
    for _, row in rxns.iterrows():
        reactants = [species_name_to_sympy[r] for r in [row.iloc[0], row.iloc[1]] if not pd.isnull(r)]
        products = [species_name_to_sympy[r] for r in [row.iloc[2], row.iloc[3]] if not pd.isnull(r)]
        forward_rate = row['Rate(forward)/s']
        backward_rate = row['Rate(back)/s']
        if backward_rate == 0 or np.isclose(backward_rate, 0):
            if len(reactants) == 1 and len(products) == 1:
                rxn = Rxn(reactants[0], products[0], forward_rate)
            elif len(reactants) == 2 and len(products) == 1:
                rxn = Rxn(reactants[0] + reactants[1], products[0], forward_rate)
            elif len(reactants) == 1 and len(products) == 2:
                rxn = Rxn(reactants[0], products[0] + products[1], forward_rate)
            else:
                rxn = Rxn(reactants[0] + reactants[1], products[0] + products[1], forward_rate)
        else:
            if len(reactants) == 1 and len(products) == 1:
                rxn = RevRxn(reactants[0], products[0], forward_rate, backward_rate)
            elif len(reactants) == 2 and len(products) == 1:
                rxn = RevRxn(reactants[0] + reactants[1], products[0], forward_rate, backward_rate)
            elif len(reactants) == 1 and len(products) == 2:
                rxn = RevRxn(reactants[0], products[0] + products[1], forward_rate, backward_rate)
            else:
                rxn = RevRxn(reactants[0] + reactants[1], products[0] + products[1], forward_rate, backward_rate)
        rxn_sys_input.append(rxn)

    rxn_sys_input.append(T(36))
    for _, row in initial_concentration.iterrows():
        rxn_sys_input.append(Conc(species_name_to_sympy[row.iloc[0]], row.iloc[1]))

    rxn_sys_input.append(sm)
    return RxnSystem(*rxn_sys_input)


def excel_to_crn(io):
    input_dfs = read_line_delimited_excel(io)
    rsys = build_bulk_crn(*input_dfs[:3])
    tmax = input_dfs[3].iloc[0][0]
    return rsys, tmax
