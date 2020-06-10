.. Examples

========
Examples
========

Basic Example
-------------

More Advanced Example
---------------------

Create a reaction system
^^^^^^^^^^^^^^^^^^^^^^^^

Reaction systems are created from various species which have specified orbitals,
names and schedules. These species are then combined into a reaction system. Below
is a sample predator-prey system that is modeled through this system.
<br/>
<br/>
TODO: Update api to follow this (move to optional parameters)
```python
sm = SpeciesManager()

x1 = sm.sp('x', {0:2})
x2 = sm.sp('y', {0:1})

rsys = RxnSystem(
    Rxn(x1 + x2, 2 * x2, 1.5),
    Rxn(x1, 2 * x1, 1),
    Rxn(x2, 1, 1),
    sm
)
```

Solve the reaction system
^^^^^^^^^^^^^^^^^^^^^^^^^

Simply specify the reaction system, the runtime, and other optional parameters for
differential equation solving to simulate the system.
```python
s = solve(rsys, time=40, max_step=1e-2)
```

#### Visualize and understand the system
Plotting species in the time domain can be done with a single command.
```python
s.basic_plot()
```
![Basic time domain plot](img/basictimedomain.png)
