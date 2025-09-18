# Digital Twin for Chemical Sciences (DTCS v.01) 

[article and citation] https://www.nature.com/articles/s43588-025-00857-y

[breifing] https://www.nature.com/articles/s43588-025-00859-w

[full article] https://chemrxiv.org/engage/chemrxiv/article-details/684a4f6d1a8f9bdab5b375a8

Directly visualizing chemical trajectories offers novel insights into catalysts, gas phase reactions, photo-induced dynamics, and quantum information processing. Identifying and tracking the exchange of matter to observe the creation and annihilation of chemical species is best achieved by closely coupling theory and experiment. We developed Digital Twin for Chemical Science (DTCS) v.01, a platform that mimics advanced characterization instruments, including those at Scientific User Facilities. DTCS v.01 addresses challenges in data acquisition, analysis, and model-driven interpretation via a physics-based, AI-accelerated approach. We validated this concept with ambient pressure X-ray Photoelectron Spectroscopy (APXPS) observations using a ubiquitous metal-water interfacial scenario, i.e., Ag/H2O, as a representative example. The inputs of DTCS v.01 are designed to mirror the experimental chemists' workflows, and the outputs can be directly compared to and are constantly updated from the experimental data. This integrated theoretical and experimental platform enhances user accessibility and facilitates the acquisition of standardized mechanistic insights.

## How to Install (new)

DTCS can be installed in Linux, MacOS, and Windows operating systems.
For Windows, we recommend using Windows Subsystem for Linux (WSL) that allows the use of Linux environment from within Windows. 

Make sure you are using Python 3.11 or above (try `python --version`)
In the current version, Python 3.11.7 is our default environment that has been tested most recently (24th April 2025)

Create a virtual environment. Run the following _in the demo directory_:
```sh
python -m venv --prompt "dtcs" .venv
source .venv/bin/activate
```

Now you should see (dtcs) in your prompt

Install the dependencies in a streamlined fashion. All dependencies and their versions are specified under the requirements.txt. For the purpose of this demo we want to make sure that you have exactly the environment this was tested on.

```sh
pip install -r ./requirements.txt
pip install "dtcs-2025.4.10.dev1+demo-py3-none-any.whl[demo]"
```

Finally, you need to add a kernel specification to your Jupyter notebook interface.

```sh
python -m ipykernel install --user --name dtcs-demo --display-name "DTCS Demo"
python -m jupyter kernelspec list
```

You should see `dtcs-demo` in that list of kernels now.
Then you can run Jupyter and open the demo notebook:

```sh
python -m jupyter notebook
```

Make sure you select `DTCS Demo` as the working kernel.

**Important Notes**

1. If you have multiple python installations, make sure you choose python 3.11 (or above) by replacing python with python3 or python3.11 resulting in the following commands.
```sh
python3 -m venv --prompt "dtcs" .venv
python3 -m ipykernel install --user --name dtcs-demo --display-name "DTCS Demo"
python3 -m jupyter kernelspec list
python3 -m jupyter notebook
```

2. If you have multiple pip installations, make sure you choose pip 3.11 (or above) by replacing pip with pip3 or pip3.11 resulting in the following commands.
```sh
pip3 install -r ./requirements.txt
pip3 install "dtcs-2025.4.10.dev1+demo-py3-none-any.whl[demo]"
```

3. If you use a macOS x86_64 (Intel), you may get an error while installing torch 2.6.0. In that case, replace torch==2.6.0 with torch==2.0.0 in the requirements.txt file and reinstall.


**If you want to use DTCS on NERSC JupyterLab** , it's highly recommended to create a seperate conda environment with
```
module load python/3.11
conda create --name dtcs-3.11 python=3.11
```

And then activate it with
```
conda activate dtcs-3.11
```
If this worked, you should see `(dtcs-3.11)` before your cursor in your shell. Then, you should install the package, with some additional dependencies.
```
pip install git+https://github.com/QianGroupPage/DTCS-v.01@demo-fixes
pip install opencv-python Pillow pygame pymatgen gpcam==7.2.5 fvgp==3.2.7
pip install numpy==1.26.2
conda install -c anaconda ipykernel
```

Finally, you need to add a kernel specification to your Jupyter notebook interface.

```
python -m ipykernel install --user --name dtcs-3.11 --display-name "DTCS 3.11"
jupyter kernelspec list
```

You should see `DTCS 3.11` in that list of kernels now, and additionally when you go to open a notebook.

## Walkthrough
This will walk you through the workflow of a simple Ag(111)-H2O system.
See reference here: https://pubs.acs.org/doi/full/10.1021/jacs.8b13672

    import sympy as sym
    from dtcs.spec.xps import XPSSpeciesManager
    from dtcs.spec.crn.bulk import CRNSpec, Rxn, RevRxn, Conc, ConcEq, ConcDiffEq, Schedule
    
First, you need to make a SpeciesManager to keep track of all your species.
Then, you can create some species.
The Species Manager lets us define species by specifying a name, as well as the 'signature' of the species we're defining, for example, the XPS binding energy locations.

    sm = XPSSpeciesManager()
    h2o_g = sm.make_species('H2O_g', 535.0, color='gray', latex='H_2O_g')
    o2_g = sm.make_species('O2_g', 535.0, color='gray', latex='O_{2g}')
    h2o = sm.make_species('H2O', 532.2, color='blue', latex='H_2O^*')
    oh = sm.make_species('OH', 0, color='red', latex='OH^*') #530.9
    o = sm.make_species('O', 0, color='aqua', latex='O^*') #530.0
    oh_h2o = sm.make_species('OH-H2O_{hb}', 0, color='black', latex='OH\!-\!H_2O^*') #531.6
    o_h2o = sm.make_species('O-H2O_{hb}', 0, color='black', latex='O\!-\!H_2O^*') #531.6
    h2o_multi = sm.make_species('multiH2O', 533.2, color='magenta', latex='H_2O_{multi}^*')

    oh_combined = sm.make_species('OH_combined', 530.9, color='red', latex='OH combined^*')  # We need to combine the OH in OH and OH-H2O
    o_combined = sm.make_species('O_combined', 530.0, color='aqua', latex='O combined^*')    # combine O in O and O-H2O
    h2o_hbond_combined = sm.make_species('H2O_{hb}_combined', 531.6, color='black', latex='H_2O_{hb} combined^*') # combine H2O_O and H2O_OH to form H2O_hb

    sm
    
Then, you define your reaction system.

    crn = CRNSpec(
        Rxn(o + h2o_g, o_h2o, k=3.915042),                  # H2O on O* Adsorption with HB (1)
        Rxn(oh + h2o_g, oh_h2o, k=1.664002),                # H2O on OH* Adsorption with HB (2)
        Rxn(o_h2o, oh + oh, k=6.220646),                    # O-H2O react to form 2 OHs (3)
        Rxn(oh + oh, o_h2o, k=0.160755),
        Rxn(oh_h2o, h2o + oh, k=0.299507),                  # OH-H2O Diffusion (4)
        Rxn(o_h2o, h2o + o, k=0.167130),                    # O-H2O Diffusion (5)
        Rxn(h2o, h2o_g, k=0.794455),                        # H2O Desorption, (6)
        Rxn(h2o_g, h2o, k=0.629363),                        # H2O Adsorption, (7)
        Rxn(oh_h2o, oh + h2o_g, k=0.300480),                # H2O on OH* Desorption with HB (8)
        Rxn(o_h2o, o + h2o_g, k=0.127713),                  # H2O on O* Desorption with HB (9)
        Rxn(oh_h2o + h2o_g, h2o_multi, k=1.267427),         # H2O on OH-H2O* Adsorption (10)
        Rxn(h2o_multi, oh_h2o + h2o_g, k=0.394500),         # H2O on OH-H2O* Desorption (11)
    
                                                                 
        ConcEq(oh_combined,oh + oh_h2o),                         
        ConcEq(o_combined,o + o_h2o),
        ConcEq(h2o_hbond_combined,o_h2o + oh_h2o),
    
        Conc(h2o_g, 1),                                     
        Conc(o, 0.25),
        sm,
        time=20,
    )


Then, you solve the system. This simulates the system for 20 time units as defined above.

    cts = crn.simulate()
    
And then to plot the APXPS spectra with selected surface species:

    cts.plot(species=[h2o, oh_combined, o_combined, h2o_hbond_combined, h2o_multi])

    plt.xlabel('Time')
    plt.ylabel('Concentration')


## More Information
DTCS v.01 is under Non-cost/Academic license. 

Unauthorized distributions are prohibited.

However, we welcome collaborations and charge no cost for academics. 

If you are interested in using DTCS v.01, please contact Dr.Jin Qian (jqian2@lbl.gov) at LBNL.

https://qiangrouppage.lbl.gov/

## Key developers to DTCS v.01 are:
Jin Qian, Sid Menon, Andrew Bogdan

We also acknowledge Asmita Jana, Ethan Crumlin, Rebecca Hamlyn, and Johannes Mahl for their insightful discussions

