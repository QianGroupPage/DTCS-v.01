.. LBL CRN docs master file

LBL Chemical Reaction Simulator
===============================

.. toctree::
   :name: mastertoc
   :maxdepth: 2
   :hidden:

   examples
   contact
   license
   api_docs/index

.. TODO: Introduction

How to Install
--------------
If you are using Anaconda::

    conda install git
    conda install git+https://github.com/rithvikp/lbl-crn

If you are using normal Python::

    pip install git+https://github.com/rithvikp/lbl-crn

Once it's installed (it might take a while), run ``python -m lblcrn --version``
to check if everything installed.

Walkthrough
-----------
This will walk you through the workflow of a predator-prey system. Import::

    from lblcrn.homogenous_crn import *

First, you need to make a SpeciesManager to keep track of all your species.
Then, you can create some species::

    sm = SpeciesManager()
    prey = sm.sp('rabbit', Orbital('1s', 535.0), {0:2})
    pred = sm.sp('fox', Orbital('1s', 535.5), {0:1})

Then, you define your reaction system::

    rsys = RxnSystem(
        Rxn(prey + pred, 2 * pred, 1.5),
        Rxn(prey, 2 * prey, 1),
        Rxn(pred, 1, 1),
        sm
    )

Then, you solve the system. This simulates the system for ``time``::

    solution = solve(rsys, time=45, max_step=0.01)
    solution.process()
    solution.basic_plot()

And then to plot the gaussian:

.. code-block::

    solution.plot_gaussian(envelope=True)

You can download this example `here <https://github.com/rithvikp/lbl-crn/blob/master/examples/predator_prey.ipynb>`

.. TODO Citation Information

Indices and Tables
------------------

- :ref:`genindex`
- :ref:`modindex`
- :ref:`search`




