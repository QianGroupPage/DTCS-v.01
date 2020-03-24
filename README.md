Surface Chemical Reaction Networks with Jupyter Notebook
============================================================
This repository is a fork of Samuel Clamson's Surface CRN code as
described below.

To use the Jupyter notebook, run `jupyter notebook` in your terminal and open `SCRN Pipeline`. Follow
the directions in the notebook.

This Surface CRN can simulate a 10 x 10 surface with fast speed. It's a good idea to use 
a grid of this size to develop a new surface or test your ideas. A 20 x 20 surface is recommended for your final results.


Simulation of Chemical Reaction Networks (CRNs) on a surface
============================================================

This package is for simulating surface CRNs, as defined in Qian & Winfree 2014 ("Parallel and Scalable Computation and Spatial Dynamics with DNA-Based Chemical Reaction Networks on a Surface").


This package provides:

* A core simulator for simulating surface CRNs on arbitrary graphs.
* Models and views for square-grid and hex-grid surface CRNs.
* A manifest-reading system to simplify setting up new surface CRNs.
* A pygame-powered GUI for visualizing square-grid and hex-grid surface CRNs.

**If you are using surface_crns on OSX, read the note for you in "Prerequisites".**

Prerequisites
=============
The GUI simulator requires pygame.

This simulator is only tested in Python 3. It will likely not run in Python 2.

**Important note for OSX users**: There is a known problem with Pygame mouse click registration on OSX with certain installations of Python. In Python 3, if you use virtualenv or a Conda installation of Python, Pygame will only register mouse movement when the mouse button is held down. You can fix this by turning off virtualenv or (if you use Conda) by using pythonw instead of python to run your simulations. If, for whatever reason, you cannot do either of these things, make sure to move your mouse *during* each click, ending the click on whatever button you want to press.

**Important note for OSX 10.14 users**: The current public stable version of pygame (1.9.6) does not work in OSX 10.14 in some verisons of Python 3. The development version of pygame does work on OSX 10.14. If you are running OSX 10.14 and you don't get a window when you run ``SurfaceCRNQueueSimulator``, try downloading the pygame repository at https://github.com/pygame/pygame, navigate into the downloaded folder, and run ``python setup.py install``.

Installation
============

Install with the following (probably requires sudo): ``pip install git+git://github.com/sclamons/surface_crns.git@master``

To update, run (also with superuser privilege): ``pip install --upgrade --no-deps git+git://github.com/sclamons/surface_crns.git@master``

If you are using OSX and you get an error while installing surface\_crns
or pygame, you may need to install several SDL libraries used by pygame.
Install them with brew (<https://brew.sh/>) with the command:

`brew install sdl sdl_image sdl_mixer sdl_ttf portmidi`

What's a surface CRN?
=====================

In brief, a surface CRN is a chemical reaction network where individual molecules are tethered to positions on a surface; or, alternatively, a surface CRN is an asynchronous cellular automata with transition rules that resemble those of unimolecular and bimolecular chemical reactions.

A surface CRN consists of a number of sites with states on an arbitrary lattice (often a grid, but possibly any graph), along with transition rules of the form ``A -> B`` or ``A + B -> C + D``, with a real-valued reaction rate for each transition rule. States change stochastically according to the transition rules for the surface CRN, with expected time specified by the reaction rate for each reaction.

For example, a site with state ``A`` can undergo ``A -> B`` to instantaneously switch to state ``B``.

Bimolecular reactions (of the form ``A + B -> C + D``) can occur only if *both* of the species before the arrow are present and next to each other on the lattice (note that there is no "directionality" inherent to the surface, so an ``A`` next to a ``B`` is equivalent to a ``B`` next to an ``A``). When a bimolecular reaction occurs, the first reactant instantaneously becomes the first product, and the second reactant simultaneously and becomes the second product, so the states ``A B`` would become ``C D`` (or, equivalently, ``B A`` would become ``D C``).

What can you do with a surface CRN? Some examples are shown in Qian and Winfree 2014; See what you can make!

How do I use the simulator?
===========================

The fastest way to start using the simulator is to use our hosted version of the simulator at http://centrosome.caltech.edu/Surface_CRN_Simulator/srv/.

To use the simulator directly, yourself, first install this package, as described under "Installation". This will install an executable script SurfaceCRNQueueSimulator, which you can run directly from the command line with ``SurfaceCRNQueueSimulator -m <manifest_file>``, where ``<manifest_file>`` is the name of a valid manifest file (see below). We have provided several example manifests in the ``examples`` folder. Please note that some of these manifests use information from other files in the ``examples`` folder. The easiest way to use these examples is to copy the entire folder.

You can also run the simulator from a Python script. To run a surface CRN specified in a file "example_manifest.txt", you would run (in your Python terminal or script)

``from surface_crns import SurfaceCRNQueueSimulator
SurfaceCRNQueueSimulator.simulate_surface_crn("example_manifest.txt")``

If you are one of the users affected by the Pygame mouse click registration bug mentioned in "Prerequisites", you can run SurfaceCRNQueueSimulator with pythonw. First, find the full name and location of SurfaceCRNQueueSimulator on your machine by running ``which SurfaceCRNQueueSimulator``. Then run ``pythonw <simulator_name> -m <manifest_file>``, where ``<simulator_name>`` is the output of the ``which`` command you ran.

What is a "manifest file"?
==========================

A manifest file is a text file that defines a surface CRN. It includes the rules of the surface CRN, an initial condition, a colormap to control how the surface CRN will be displayed, and a number of run options.

How do I write a manifest file?
===============================
You can find manifest files for several example surface CRNs in the "examples" folder and its subfolders (.txt files). There are also a number of examples available at http://www.dna.caltech.edu/Surface_CRN_Simulator/srv/.

A surface chemical reaction network (sCRN) on a square grid can be specified by a manifest file consisting of four parts:

* Transition rules
* Initial state
* Colormap
* General Settings

Any line in a manifest beginning with a '#' will be treated as a comment and ignored.

## Transition Rules

The section of a manifest specifying transition rules begins with the line "!START_TRANSITION_RULES" and ends with the line "!END_TRANSITION_RULES" (no quotation marks). A single transition rule is specified by a line of the form

(RATE) NAME1 -> NAME2

or

(RATE) NAME1 + NAME2 -> NAME3 + NAME4

where "RATE" is a (nonnegative) number designating the reaction rate of the transition rule, and NAME1, NAME2, NAME3, and NAME4 are alphanumeric species identifiers. Species identifiers may contain letters, numbers, commas, and underscores, but not other punctuation or whitespace. For example, a rule that species X1 and Y1 react to form species X2 and Y2 with rate 10 (in arbitrary units of 1/time) would be specified with:

(10) X1 + Y1 -> X2 + Y2

Rates may be placed anywhere in a transition rule statement except inside an identifier or the "->" token. For example, the following transition rules are all valid:

(10) X + Y -> A + B
X + Y ->(10) A + B
X + Y -> A + B (10)
Initial state
The section of a manifest specifying the initial state of the surface begins with the line "!START_INIT_STATE" and ends with the line "!END_INIT_STATE". An initial state is specified in a comma- or whitespace-separated format, with rows separated by newlines and columns are separated by either whitespace or commas (either is acceptable).

The contents of each cell in the grid are specified just as in the transition rule section. The name used in the initial state section must exactly match the name used in the transition rule section, or transition rules will not be applied to that position.

## Colormap

The section specifying a colormap begins with the line "!START_COLORMAP" and ends with the line "!END_COLORMAP". The colormap determines the colors that will be used to represent each species, and can be used to group multiple species to be displayed as one color. The colormap is optional; if no colormap is specified, then colors will be assigned automatically.

A single species' color can be specified with a line of the form

NAME: (R, G, B)

where "NAME" is the name of the species (of the same format used to specify species names in the transition rules and initial state sections) and "R", "G", and "B" are integers between 0 and 255 (inclusive) specifying the red, green, and blue components of the specified color. For instance, to make species "X1" display as black, one would write

X1: (0, 0, 0)

and to make species "Y1" display as vibrant red, one would write

Y1: (255, 0, 0)

Multiple species can also be grouped to be displayed with the same color and given a single display name using a statement of the form

{CLASS} NAME1, NAME2,...: (R, G, B)

where "CLASS" is the (optional) name under which the species should be displayed, "NAME1, NAME2,..." is a comma-separated list of species names of any length, and "R", "G", and "B" are as above. For instance, to display all of the species X1, X2, X3, and X4 as red and all of the species Y1, Y2, Y3, and Y4 as green, one would write

{X} X1, X2, X3, X4: (255, 0, 0)
{Y} Y1, Y2, Y3, Y4: (0, 255, 0)

## General Settings

Any line not in a transition rule, initial state, or colormap block is interpreted as a general setting. General settings take the form

SETTING = VALUE

with the obvious meanings. Useful settings to know are:

* speedup_factor: A nonnegative real number dictating the speed of simulation playback. Larger numbers mean faster playback.
* rng_seed: Integer specifying the random number seed used by the simulation. Set this value to make simulations reproducible.
* max_duration: A nonnegative number specifying the maximum length of simulation in arbitrary time units (the same arbitrary time units specified by transition rule reaction rates).
* node_display: Determines whether the state of each position on the grid (node) is overlaid, in text, on that node. Set to "text" to overlay text, or "color" to only show node color (default)
* pixels_per_node: Determines the size of a node, in pixels.
* draw_cell_borders: Iff True, black lines will be drawn around the edges of each cell.
* grid: Toggles between a square grid (default, "square") or a hex grid ("hex").

Acknowledgements
================

This package uses the pygbutton package by Al Sweigart (https://github.com/asweigart/pygbutton), as well as the random_color function made by adews (https://gist.github.com/adewes/5884820.). Thanks to both authors!
