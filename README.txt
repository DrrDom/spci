SPCI
-----
Tool for mining structure-property relationships from chemical datasets

Description
-----
Retrieves structure-property relationship from datasets in a chemically meaningful way.
Returns estimated contributions of fragments to the investigated property of compounds from a data set and can estimate contribution of different physicochemical factors as well.

Features
-----
1. Easy to use straightforward workflow with GUI.
2. Automatic model building and cross-validation.
3. Build models for imbalanced data set using multiple oversampling approach.
4. Model prediction.
5. Several fragmentation schemes to compute fragment contributions of:
  a) ring systems;
  b) Murcko scaffolds;
  c) common functional groups and rings;
  d) user-defined fragments;
  e) automatically generated fragments (based on SMARTS pattern matched broken bonds).

Visualization of results
-----
1. Built-in visualization.
2. rspci - R package for custom visualization (https://github.com/DrrDom/rspci)
3. Online tool for visualization, plot customization and figure downloading (link will be published soon)

Manual
-----
Short manual is included.

Publications
-----
Structural interpretation was published in the paper http://dx.doi.org/10.1002/minf.201300029
Integrated structural and physicochemical interpretation was published in the paper http://dx.doi.org/10.1021/acs.jcim.6b00371

Citation
-----
If you use this tool please cite this repository and the publication J. Chem. Inf. Model.  56, 8, 1455-1469 (http://dx.doi.org/10.1021/acs.jcim.6b00371)

Home page
-----
http://qsar4u.com/pages/sirms_qsar.php

License
-----
GPLv3

What'snew
-----
**1.0.0**
- RDKit is used as a backend instead of Indigo
- multiple undersampling was implemented
- changed default descriptors, that make this version incompatible with previous models and vice versa.
- updated sirms descriptors
- many small fixes and improvements

