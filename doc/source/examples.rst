.. Examples

========
Examples
========

Bulk CRN Predator Prey
----------------------

This will be a step-by step explanation of setting up a basic predator-prey scenario.
First, import the library::

    from lblcrn.crn_sym import *
    from lblcrn.experiments.simulate *

Then make your species manager. This is the object that keeps track of all
the information relating to species. See the docs here (TODO)::
    
    sm = SpeciesManager()

Using the species manager, you can declare as many species as desired as follows::

    x1 = sm.sp('x')
    x2 = sm.sp('y')


Then, given these species, arbitrary reaction systems can be created using the :code:`Rxn`,
:code:`Conc` and other classes detailed here (TODO)::

    rsys = RxnSystem(
        Rxn(x1 + x2, 2 * x2, 1.5),
        Rxn(x1, 2 * x1, 1),
        Rxn(x2, None, 1),
        sm
    )

Simply specify the reaction system, runtime, and other optional parameters to simulate the system::

    s = simulate(rsys, time=40, max_step=1e-2)

Given the simulated solution, a variety of different visualizations and analyses can be performed.
For example, the time series data can be plotted with the plot command::

    s.plot()

.. image:: _static/img/predator_prey_time_series.png
    :width: 400
    :alt: Predator-Prey Time Series data


Bulk CRN Something Else
-----------------------
