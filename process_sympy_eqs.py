def process_sympy_eqs(eqs):
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

    for eq in eqs:
        if is_function(eq):
            # free_symbols in sympy/core/basic returns a set.
            symbols = list(eq.rhs.free_symbols)
            # Lamdify is in sympy.utilities.lambdify
            # The regular format for a python lambda is lambda x, y, z: x + 5.
            # This assumes no ':' in any variable name. Python variable names have to
            # be a combination of uppercase and lowercase letters and the underscore character _,
            # and a digit as long as it's not the first character.
            lambda_expr = get_formula_from_lambstr(lambdastr(symbols, eq.rhs))

            functions[eq.lhs.free_symbols[0]] = (symbols, lambda_expr)
        elif is_initial_value(eq):
            #TODO: solve a system of equations for initial values
            #TODO: grab the numerical values
            initial_values[eq.lhs.free_symbols[0]] = sympify(eq.rhs).evalf
        elif is_schedule(eq):
            schedules[eq.lhs.free_symbols[0]] = to_tuple(eq.rhs)
        else:
            raise ValueError("Your equation is in the wrong format.")

        return functions, initial_values, schedules


###### Helper methods. ######
def get_formula_from_lambstr(lambstr):
    #  Here's lambdastr's return statement "lambda %s: (%s)" % (args, expr)"
    #  This function takes the returned string, grabs the expr, and return it.

    # The regular format for a python lambda is lambda x, y, z: x + 5.
    # This assumes no ':' in any variable name. Python variable names have to
    # be a combination of uppercase and lowercase letters and the underscore character _,
    # and a digit as long as it's not the first character.
    return lambstr.split(':', 1)


def to_tuple(collections):
    if isinstance(collections, list):
        return collections
    elif isinstance(collections, dict):
        result = []
        for k, v in collections.items():
            result.append((k, v))
        return result
    else:
        raise ValueError("Schedule has to either a dictionary or a list of tuples.")


###### Criteria for separating Sympy equations. ######
def is_function(eq):
    return expr_is_one_symbol(eq.lhs) and not is_schedule(eq)


def is_initial_value(eq):
    # eq(a b c, 0.001) sets all initial values to 0.0001
    return eq.rhs.is_number


def is_schedule(eq):
    return expr_is_one_symbol(eq.lhs) and expr_is_time_seq(eq.rhs)


def expr_is_time_seq(expr):
    raise NotImplementedError


def expr_is_one_symbol(expr):
    # Atoms is in sympy/core/ basic.py
    # atoms return a set ... as of now; if their code changes, we may have to follow suit.
    return len(atoms(expr)) == 1 and len(expr.free_symbols) == 1
