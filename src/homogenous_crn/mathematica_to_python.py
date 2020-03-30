syntactical_dict = {
    ';': '',  # delete all Mathematica line separators
    '[': '(',
    ']': ')',
    '{': '[',
    '}': ']',
    ' ,': ', ',
    'rxn': 'crn.rxn',
    'conc': 'crn.conc',
    'SimulateRxnsys': 'crn.simulate_rxnsys',
    # Additional syntax for the H2O input file
    'AccountingForm[y2[tmax] /. sol]\nAccountingForm[x3[tmax] /. sol]\nAccountingForm[x4[tmax] /. sol]\nAccountingForm[x5[tmax] /. sol]\nAccountingForm[x6[tmax] /. sol]':
        'results.append([T, P, sol[-1, 5], sol[-1, 0], sol[-1, 1], sol[-1, 2], sol[-1, 3]])',
    "T=": "T = ",
    "P=": "P = "
}

PRE_LOGUE = """
import crn_sym as crn
import csv
results = []
"""

POST_LOGUE = """
with open("H2O_CRN_output_results.csv","w+") as my_csv:
    csvWriter = csv.writer(my_csv, delimiter=',')
    csvWriter.writerow(crn.species_in_sys(rsys))
    csvWriter.writerows(results)
"""

import re, sys, os
from .util import multiple_replace

def transform_to_mathematica(filename):
    """
    Parse Mathematica notebook to Python and save to a Python file with the same name.

    :param filename: name of the file ro transform.
    :return: None
    """
    f = open(filename, "r")
    name = os.path.splitext(filename)[0]
    new_file = open("{}.py".format(name), "w")

    text = multiple_replace(syntactical_dict, f.read())

    # Process P = AccountingForm(0.00001) into P = 0.00001
    regex = re.compile("AccountingForm\(([0-9]+\.*[0-9]+)\)")
    text = regex.sub(lambda mo: mo.group(1), text)

    # Process bad Python format ,12 into , 12
    regex = re.compile(",([0-9])")
    text = regex.sub(lambda mo: ", %s" % mo.group(1), text)

    # Place all formulas into strings, as specified in Ye's CRN transation
    regex = re.compile("([a-zA-Z][0-9] \+ [a-zA-Z][0-9])")
    text = regex.sub(lambda mo: "\"%s\"" % mo.group(1), text)

    # Then quote the remaining singleton symbols
    regex = re.compile("([a-zA-Z][0-9])(?![(? \+)|\"])")
    text = regex.sub(lambda mo: "\"%s\"" % mo.group(1), text)

    new_file.writelines(PRE_LOGUE)
    new_file.writelines(text)
    new_file.writelines(POST_LOGUE)

    f.close()
    new_file.close()


if __name__ == "__main__":

  for filename in sys.argv[1:]:
      if filename.endswith('.nb'):
          transform_to_mathematica(filename)
      else:
          print "%s is not a Mathematica notebook" % filename
