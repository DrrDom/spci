# SPCI

Automatic tool for mining structure-property relationships from chemical data sets

#### Description

Retrieves structure-property relationship from data sets in a chemically meaningful way.  
Returns estimated contributions of fragments to the investigated property of compounds from a data set and can estimate contribution of different physicochemical factors as well.

#### Features

1. Easy to use straightforward workflow with GUI.
2. Automatic model building and cross-validation.
3. Build models for imbalanced data set using the multiple oversampling approach.
4. Prediction with built models.
5. Several fragmentation schemes to compute fragment contributions of:
  - common functional groups and rings;  
  - Murcko scaffolds;  
  - user-defined fragments;  
  - automatically generated fragments (based on SMARTS pattern matching broken bonds);  
  - per atom fragmentation.

#### Visualization and analysis of results

1. Built-in visualization.
2. rspci - R package for custom visualization and analysis (https://github.com/DrrDom/rspci)
3. Online tool for visualization, plot customization and figure downloading (http://158.194.101.252:3838/spci-vis/). Demo version is here (http://158.194.101.252:3838/spci-vis-demo/)
4. Per atom contributions can be visualized with RDKit similarity maps.
#### Manual

The short manual is included.

#### Citation

1.	Polishchuk, P. G.; Kuz'min, V. E.; Artemenko, A. G.; Muratov, E. N., Universal Approach for Structural Interpretation of Qsar/Qspr Models. Mol. Inf. 2013, 32, 843-853 - http://dx.doi.org/10.1002/minf.201300029 - structural interpretation.
2.	Polishchuk, P.; Tinkov, O.; Khristova, T.; Ognichenko, L.; Kosinskaya, A.; Varnek, A.; Kuz’min, V., Structural and Physico-Chemical Interpretation (SPCI) of QSAR Models and Its Comparison with Matched Molecular Pair Analysis. J. Chem. Inf. Model. 2016, 56, 1455-1469 - http://dx.doi.org/10.1021/acs.jcim.6b00371 - integrated structural and physicochemical interpretation.

#### Home page

http://qsar4u.com/pages/sirms_qsar.php

#### License

LGPLv3

#### What's new

1.0.0 (03.07.2018)
- RDKit is used as a backend instead of Indigo
- multiple undersampling was implemented
- changed default descriptors, that make this version incompatible with previous models and vice versa.
- updated sirms descriptors
- many small fixes and improvements

1.1.0 (07.02.2021)
- added support of RDKit descriptors
- added per atom fragmentation
- reorganized as a Python package
- console scripts have prefix spci_*

1.1.1 (23.03.2021)
- changed license to LGPLv3
- fixed arguments in scpi_descriptors
