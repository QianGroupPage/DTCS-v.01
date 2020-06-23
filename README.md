# lbl-crn
A chemical reaction network solver

## How to Install
if you are using Anaconda:

    conda install git
    conda install git+https://github.com/rithvikp/lbl-crn

if you are using normal Python:

    pip install git+https://github.com/rithvikp/lbl-crn

Once it's installed (it might take a while), run
`python -m lblcrn --version`
to check if everything installed.

## Walkthrough
This will walk you through the workflow of a predator-prey system.

    from lblcrn import *
    
First, you need to make a SpeciesManager to keep track of all your species.
Then, you can create some species.

    sm = SpeciesManager()
    prey = sm.sp('rabbit', Orbital('1s', 1.0))
    pred = sm.sp('fox', Orbital('1s', 2))
    
Then, you define your reaction system.

    rsys = RxnSystem(
        sm,

        Rxn(prey + pred, 2 * pred, 2),
        Rxn(prey, 2 * prey, 1),
        Rxn(pred, None, 1),

        Conc(pred, 1),
        Conc(prey, 1),
    )

Then, you solve the system. This simulates the system for 45 time units.

    time_series = simulate_crn(rsys, time_max=45, max_step=0.01)
    time_series.plot()
    
And then to plot the gaussian:

    time_series.xps.plot()

You can download this example [here](https://github.com/rithvikp/lbl-crn/blob/master/examples/predator_prey.ipynb)
(except not right now you can't because we're editing stuff).

## More Information
For more information, a quickstart guide, examples, and the api reference, 
please refer to the [documentation]()
