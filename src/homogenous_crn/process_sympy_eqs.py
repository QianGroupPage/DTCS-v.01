from typing import Callable, Mapping, Collection, Sequence, Mapping, List, Tuple, Collection
import sympy as sym
from sympy.utilities.lambdify import lambdastr
from .util import multiple_replace
from .species import Species

def process_sympy_eqs(eqs: Sequence[sym.Eq] , species_ordering: Sequence[str] = None) -> Tuple[Callable, Mapping, Mapping]:
    # Process a list of eqs consisting of eqs of the following types.

    # 1. Initial concentrations:
    # Eq(x, 0.001),
    # Eq(y, 0.0001),

    # 2. Functions
    # Eq(x, x + 5 * y)
    # Eq(y, -3 * x + 6 * y)

    # 3. Schedules
    # The right hand side has to be some key-value pairs, or tuples
    # Eq(x, [(0: 0.001, 1: 0.007, 100: 0.0007)])
    # Eq(y, {0: 4, 70: 10, 110: -2})
    # Eq(z, [(50, 10), (40, -9)])

    # Development hint:
    # Definition, methods, and docs for Sympy's eq objects are all in sympy/core/relational.
    # you'll see that EQ is a sublcass of relational.

    functions = {}  #  A dictionary that maps variable names to a tuple of (Sympy symbol, mathematical
    # formulas as a string).
    initial_values = {}  #  A dictionary that maps variable names to numbers.
    schedules = {}  #  A dictionary that maps variable names to a list of (time, concentration) tuples.

    all_rhs_symbols = []
    all_lhs_symbols = []

    for eq in eqs:
        # free_symbols in sympy/core/basic returns a set.
        rhs_symbols = list(eq.rhs.free_symbols)
        lhs_symbols = list(eq.lhs.free_symbols)

        all_rhs_symbols.extend(rhs_symbols)
        all_lhs_symbols.extend(lhs_symbols)

        if is_function(eq):
            # Lamdify is in sympy.utilities.lambdify
            # The regular format for a python lambda is lambda x, y, z: x + 5.
            # This assumes no ':' in any variable name. Python variable names have to
            # be a combination of uppercase and lowercase letters and the underscore character _,
            # and a digit as long as it's not the first character.
            lambda_expr = get_formula_from_lambstr(lambdastr(rhs_symbols, eq.rhs))
            functions[repr(lhs_symbols[0])] = (rhs_symbols, lambda_expr)
        elif is_initial_value(eq):
            #TODO: solve a system of equations for initial values
            initial_values[repr(lhs_symbols[0])] = sym.core.sympify(eq.rhs).evalf
        elif is_schedule(eq):
            schedules[repr(lhs_symbols[0])] = to_tuple(eq.rhs)
        else:
            raise ValueError("Your equation is in the wrong format.")

    all_rhs_symbols = [repr(s) for s in list(set(all_rhs_symbols)) if s != 't']
    all_rhs_symbols.sort()
    all_lhs_symbols = [repr(s) for s in list(set(all_lhs_symbols)) if s != 't']
    all_lhs_symbols.sort()
    assert len(all_lhs_symbols) == len(all_rhs_symbols), "Symbols with no default values are not supported."

    if species_ordering is None:
        all_symbols = all_rhs_symbols
    else:
        assert len(all_rhs_symbols) == len(species_ordering), "The optional argument is incomplete"
        all_symbols = species_ordering

    syntactical_dict = {}
    for i in range(len(all_symbols)):
        syntactical_dict[all_symbols[i]] = "y[{}]".format(str(i))

    return_values = []
    for i in range(len(all_symbols)):
        replaced_lambda_expr = multiple_replace(syntactical_dict, functions[all_symbols[i]][1][1])
        return_values.append(replaced_lambda_expr)

    derivative_function_str = "lambda t, y: [{}]".format(','.join(return_values))
    derivative_function = eval(derivative_function_str, None, None)

    return derivative_function, initial_values, schedules


def rxns_to_python_derivative_function(rxns):
    """
    Return a Python ode function corresponding to the reaction system interpretable.
    """
    # A dictionary that maps species names to its index
    species_ordering = rxns.symbol_index

    # A list of ode formulas following the expressions.
    odes = rxns.get_ode_expressions()

    syntactical_dict = {}
    for species, i in species_ordering.items():
        syntactical_dict[str(species)] = "y[{}]".format(str(i))

    return_values = [None for _ in range(len(rxns.symbol_index))]
    for species, i in species_ordering.items():
        symbols = list(odes[i].free_symbols)
        lambda_expr = get_formula_from_lambstr(lambdastr(symbols, odes[i]))

        replaced_lambda_expr = multiple_replace(syntactical_dict, lambda_expr[1])
        return_values[i] = replaced_lambda_expr

    derivative_function_str = "lambda t, y: [{}]".format(','.join(return_values))
    derivative_function = eval(derivative_function_str, None, None)

    return derivative_function

def rxns_to_substances(rxns):
    substances = [None for _ in range(len(rxns.symbol_index))]
    for species, i in rxns.symbol_index.items():
        substances[i] = Species(str(species), None)
    return substances

def rxns_to_initial_values(rxns):
     symbol_index = rxns.symbol_index
     init_vals = [0.0 for _ in range(len(symbol_index))]
     for conc_eq in rxns.conc_eqs:
         init_vals[symbol_index[conc_eq.species]] = sym.sympify(conc_eq.expression).evalf()
     return init_vals

# def rxns_to_schedule(rxns):
    # i
    # for schedule in rxns.schedule:



###### Helper methods. ######
def get_formula_from_lambstr(lambstr: str):
    #  Here's lambdastr's return statement "lambda %s: (%s)" % (args, expr)"
    #  This function takes the returned string, grabs the expr, and return it.

    # The regular format for a python lambda is lambda x, y, z: x + 5.
    # This assumes no ':' in any variable name. Python variable names have to
    # be a combination of uppercase and lowercase letters and the underscore character _,
    # and a digit as long as it's not the first character.
    return lambstr.split(':', 1)


def to_tuple(collections: Collection) -> List[Tuple]:
    if isinstance(collections, list):
        return collections
    elif isinstance(collections, dict):
        result = []
        result = []
        for k, v in collections.items():
            result.append((k, v))
        return result
    else:
        raise ValueError("Schedule has to either a dictionary or a list of tuples.")


###### Criteria for separating Sympy equations. ######
def is_function(eq: sym.Eq) -> bool:
    return expr_is_one_symbol(eq.lhs) and not is_schedule(eq)


def is_initial_value(eq: sym.Eq) -> bool:
    # eq(a b c, 0.001) sets all initial values to 0.0001
    return eq.rhs.is_number


def is_schedule(eq: sym.Eq) -> bool:
    return expr_is_one_symbol(eq.lhs) and expr_is_time_seq(eq.rhs)


def expr_is_time_seq(expr: sym.Expr) -> bool:
    # If the right side is a dictionary or a list, it represents the sequence.
    return isinstance(expr, list) or isinstance(expr, dict)


def expr_is_one_symbol(expr: sym.Expr) -> bool:
    # Atoms is in sympy/core/ basic.py
    # atoms return a set ... as of now; if their code changes, we may have to follow suit.
    return len(expr.atoms()) == 1 and len(expr.free_symbols) == 1
