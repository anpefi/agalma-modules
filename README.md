# agalma-modules
My own collection of agalma modules. Agalma (https://bitbucket.org/caseywdunn/agalma) is a set of analysis pipelines for phylogenomic analyses developed by the Dunn Lab at Brown University.

## Installation and use
These modules require [agalma](https://bitbucket.org/caseywdunn/agalma). Please refer to the agalma's documentation for details on the pipeline model used. To install them just clone this repository and copy/link the files to a directory in your PATH (maybe to ~/bin). Then you can call the module with the python interpreter.
```
git clone https://github.com/anpefi/agalma-modules.git
ln -s agalma-modules/agalma/* ~/bin/

#to run
python <module_name.py> <options>
```
Note that some modules could require some extra third-party software.


## Modules

## homologizeSRV.py
This module is a version of the original agalma's homologize.py with three major changes: blastn is used for nucleotide sequences instead of tblastx, a value for the inflation parameter of MCL can be passed, and BLAST hits are filtered by using SRVs (Bit-Score Ratio Values, the bit score relative to the best)as threshold instead of the raw 200 bitscore used in homologize.

## orthofinder.py
This module uses [OrthoFinder 0.4](http://www.stevekellylab.com/software/orthofinder) to get clusters of homologous genes (orthogroups in the terminology of OrthoFinder), that could be incorporated into the phylogeny pipelines of agalma. Requires orthofinder.py in a directory of the PATH with execution permissions (maybe you should add ""#!/bin/env python" at the beginning of orthofinder.py).
