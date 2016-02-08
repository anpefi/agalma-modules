# agalma-modules
This is my own collection of agalma modules (or pipelines using the agalma terminology). Agalma (https://bitbucket.org/caseywdunn/agalma) is a set of analysis pipelines for phylogenomic analyses developed by the Dunn Lab at Brown University. These modules were created to satisfy our needs, although they could be useful for agalma's general users.

## Installation and use
These modules require that [agalma](https://bitbucket.org/caseywdunn/agalma) is correctly installed. Please refer to the agalma's documentation for details on the pipeline model used. To install the modules just clone this repository (it requires git) and copy/link the files to a directory in your PATH (maybe ~/bin). Then you can call the module with the python interpreter.
```
git clone https://github.com/anpefi/agalma-modules.git
ln -s agalma-modules/agalma/* ~/bin/

# Running command
python <module_name.py> <options>
```
*Note that some modules could require some extra third-party software.*


## Modules

## homologizeSRV.py
This module is a version of the original homologize in agalma with three major changes: blastn is used for nucleotide sequences instead of tblastx, a value for the inflation parameter of MCL can be passed, and BLAST hits are filtered by using SRVs (Bit-Score Ratio Values, the bit score relative to the best)as threshold instead of the raw 200 bitscore used in homologize.
This module should be used as an alternative to the original homologize module in agalma, at the beginning of the phylogenomic pipeline.

## orthofinder.py
This module uses [OrthoFinder 0.4](http://www.stevekellylab.com/software/orthofinder) to get clusters of homologous genes (orthogroups in the terminology of OrthoFinder), that could be incorporated into the phylogeny pipelines of agalma. **Requires that the main program, orthofinder.py be in a directory of the PATH with execution permissions** (maybe you should add "#!/bin/env python" at the beginning of orthofinder.py).
This module should be used as an alternative to the original homologize module in agalma, at the beginning of the phylogenomic pipeline.
