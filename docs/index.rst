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

    conda install pip
    conda install git
    pip install git+https://github.com/rithvikp/lbl-crn

If you are using normal Python::

    pip install git+https://github.com/rithvikp/lbl-crn

Once it's installed (it might take a while), run ``python -m lblcrn --version``
to check if everything installed.

API Rundown
-----------

For exhaustive documentation, see the full :ref:`API Docs <api_docs>`

- Help Utilities
    - ``lblcrn_help()``
    - ``lblcrn_docs()``
    - ``lblcrn_echo_on()`` and ``lblcrn_echo_off()``
    - ``lblcrn_examples.load()``
- Monty JSON save/load utilities
- ``RxnSystem``, which takes as input:
    - one ``SpeciesManager``
    - any amount of: ``Rxn``, ``RevRxn``, ``Schedule``, ``Conc``, ``ConcEq``, ``ConcDiffEq``
- ``simulate_crn()``, which takes a ``RxnSystem`` and creates a ``CRNTimeSeries``
- ``simulate_xps()``, which takes a ``RxnSystem`` and creates an ``XPSExperiment``
- ``CRNTimeSeries``
    - ``.df``, a pandas DataFrame
    - ``.at()``, to look at one time
    - ``.xps_with()`` to make an ``XPSExperiment``
    - ``.plot()``
- XPSExperiment
    - ``.df``, a pandas DataFrame
    - ``.set_experimental()``
    - ``.set_gas_interval()``
    - ``.set_scale_factor()``
    - ``.plot()``

Walkthrough
-----------
This will walk you through the workflow of a predator-prey system. Import::

    from lblcrn import *

First, you need to make a SpeciesManager to keep track of all your species.
Then, you can create some species::

    sm = SpeciesManager()
    prey = sm.sp('rabbit', Orbital('1s', 1.0))
    pred = sm.sp('fox', Orbital('1s', 2))

Then, you define your reaction system::

    rsys = RxnSystem(
        sm,

        Rxn(prey + pred, 2 * pred, 2),
        Rxn(prey, 2 * prey, 1),
        Rxn(pred, None, 1),

        Conc(pred, 1),
        Conc(prey, 1),
    )

Then, you solve the system. This simulates the system for ``time``::

    time_series = simulate_crn(rsys, time_max=45, max_step=0.01)
    time_series.plot()

And then to plot some x-ray spectroscopy gaussians::

    time_series.xps_with(t=20).plot()

You can access this example through::

    lblcrn_examples.load('predator-prey')

.. TODO Citation Information

Indices and Tables
------------------

- :ref:`genindex`
- :ref:`modindex`
- :ref:`search`




