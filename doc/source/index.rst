.. LBL CRN docs master file

LBL Chemical Reaction Simulator
===============================

.. toctree::
   :name: mastertoc
   :maxdepth: 2
   :hidden:

   bulk_examples
   contact
   license

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

For exhaustive documentation, see the full API Docs TODO

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

.. TODO Citation Information

Indices and Tables
------------------

- :ref:`genindex`
- :ref:`modindex`
- :ref:`search`




