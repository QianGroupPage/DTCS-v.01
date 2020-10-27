.. Examples

========
Examples
========

In-Depth Example
----------------

This will be a step-by step explanation of setting up a basic predator-prey scenario.
First, import the library::

    from lblcrn import *

Then make your species manager. This is the object that keeps track of all
the information relating to species. See the docs here (TODO)::
    
    sm = SpeciesManager()

Using the species manager, you can declare as many species as desired as follows::

    x1 = sm.sp('x', {0:2})
    x2 = sm.sp('y', {0:1})


Then, given these species, arbitrary reaction systems can be created using the :code:`Rxn`,
:code:`Conc` and other classes detailed here (TODO)::

    rsys = RxnSystem(
        Rxn(x1 + x2, 2 * x2, 1.5),
        Rxn(x1, 2 * x1, 1),
        Rxn(x2, 1, 1),
        sm
    )

Simply specify the reaction system, runtime, and other optional parameters to simulate the system::

    s = simulate_crn(rsys, time=40, max_step=1e-2)

Given the simulated solution, a variety of different visualizations and analyses can be performed.
For example, the time series data can be plotted with the plot command::

    s.plot()

.. image:: _static/img/predator_prey_time_series.png
    :width: 400
    :alt: Predator-Prey Time Series data


Reaction Conditions
-------------------

.. code-block::

    from lblcrn import *
    import sympy as sym

    sm = SpeciesManager()
    x1 = sm.sp('sinewave', Orbital('1s', 1))
    x2 = sm.sp('lagged', Orbital('1s', 2))
    x3 = sm.sp('laggy chem', Orbital('1s', 2))
    x4 = sm.sp('chemical', Orbital('1s', 2))

    rsys = RxnSystem(
        sm,

        Rxn(x1, x4, k=0.01),
        Rxn(x2, x4, k=0.01),
        Rxn(x3, x4, k=0.1),

        ConcEq(x1, -sym.sin(T)),
        Conc(x2, 0.5),
        ConcDiffEq(x2, x1),
        Term(x3, x1),
    )

    solution = simulate(rsys, time_max=45, max_step=0.01)
