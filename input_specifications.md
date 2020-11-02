## Input Specifications

The `xps_data_processing` module assumes that all files generated from an experiment are stored in a single
folder. In addition to machine-generated measurement files, the folder contains a manually 
entered digital notebook.

Any file whose name contains "digital" followed immediately by "notebook" is recognized as the digital notebook, for 
instance, "digital_notebook" or "Digital  Notebook 2020-08-20". Both `_` and space are recognized as filename 
separators.

All filenames for measurement files are divided into sections separated by `_`. The last section is 
the measurement number corresponding to a line describing the corresponding measurement
in the digital notebook. For instance, `Ir_foil_20200708_0003.txt`

### Digital Notebook

#### General Comments Blocks

The processor reads the first several blocks as general comment blocks until 
it finds any block in which a line starts with a number.


#### Condition Block

The first line contains temperature, pressure, comments, each separated by `,`, `.`, or `:`.

The following lines each describe a specific measurement under this condition. Each line must 
be in the following format: 

`[measurement number] +  [space] + [species names in the measurement separated by commas]`, and 
optionally followed by `[colon] + [additional comments separated by commas]`.

For example, `036 Survey, Ir 4f, O 1s, VB at 735 eV: `


#### Comment Blocks between Condition Blocks

To add a comment block between two condition line blocks, don't start any line 
with a number. The system will take such a block as a comment for an event between 
the previous condition line block and the next (TODO, currently additional comments
are skipped).






