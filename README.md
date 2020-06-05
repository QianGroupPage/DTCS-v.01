# lbl-crn
A chemical reaction network solver

## How to Install
The easy install is just
`pip install git+https://github.com/rithvikp/lbl-crn`

## Testing
To test if everything is working properly
1. Download this [example](https://github.com/rithvikp/lbl-crn/blob/master/examples/predator_prey.ipynb)
2. Open the file in Jupyter Notebook and run it - if you see all the graphs and
whatnot, the install worked.

## Structure
``` bash
lbl_crn/
├── common/ # code used by all submodules
│   ├── utils.py # miscellaneous functions
│   └── b.py
│   
├── surface_crn/
│   ├── examples/
│   │   ├── a.ipynb
│   │   └── b.py
│   ├── a.py
│   └── b.py
│   
├── homogenous_crn/
│   ├── examples/
│   │   ├── a.ipynb
│   │   └── b.py
│   └── b.py
│   
├── script.sh
├── otherscript.sh
└── requirements.txt
```

## TODO
  - Decide on repo structure
