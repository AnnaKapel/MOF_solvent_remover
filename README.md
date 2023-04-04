# MOF Solvent Remover

This program was designed to clean CIF files of Metal-Organic Frameworks (MOFs) from free solvents, counterions and bound solvents keeping into account the charge of removed fragments.

## Installation
This code requires CSD Python API ([https://www.ccdc.cam.ac.uk/solutions/csd-core/components/csd-python-api/](https://www.ccdc.cam.ac.uk/solutions/csd-core/components/csd-python-api/)). You can install it from the source.

Other dependencies:

 - pandas
 - mendeleev

## File requirements
Program was designed to work with CIF files.

 - **Filename in format REFCODE.cif or REFCODE_xxx.cif** (the script accesses CSD to determine presence of terminal oxygens). With invalid refcode all terminal oxygens will be removed.
 - **P1 symmetry** (You can convert your MOFs using Materials Studio or get them in P1 from CSD using this code https://github.com/uowoolab/CSD-cleaner - 100% compatible)
 - **MOF is considered a periodic structure**

## Solvent removal options - flags

 - `--data_path` Provide the path to the folder with input MOFs
 - `--export_path` Provide path to the output folder
 - `--n_processes` Specify number of processes for parallelization
 - `-v`, `--verbose` Pass to mute command line output
 - `--keep_bound` Keeps the bound solvent in the MOF
 - `--keep_oxo` Keeps all the terminal oxygens
## Outputs
**MOFs_removed_solvent folde**r contains all the edited structures that had solvent removed

If --keep_bound - **Free_solvent_removal_results.csv** - description of columns:

 - **CIF** - Filename
 - **Solvent** - YES or NO if no solvent identified
 - **Free solvent** - list of atoms grouped by molecules in free solvent
 - Number of free solvent molecules
 - **Counterions** - list of atoms grouped by molecules in counterions (charged fragments)
 - **Number of counterions**
 - **Charge_removed** - total charge removed from structure, can be positive, negative or 0
 - **Total atoms** - total atoms detected as solvent
 - **Atoms removed** - total atoms removed from CIF by the parcer
 - **atoms_match_flag** - TRUE if the cif parcer left out any atoms -> Total atoms and Atoms removed don't match
 - **Metal counterion flag** - TRUE if identified counterion that contains metals -> higher chance of wrong charge assignment
 - **huge_counterion_flag** - TRUE if there are structures with counterions with a lof of metals and high mass that the code fails to assign charges to 

If all solvent should be removed - **Solvent_removal_results.csv** - additional columns:

 - **Bound solvent** -  list of atoms grouped by molecules in bound solvent
 - **Number of bound molecules**
 - **flag_double** - flags if there are double bonds near the binding cite
 - **flag_aromatic** - flags if there are aromatic solvents removed
 - **Entry terminal oxo** - the script checks the CSD entry to determine if there are terminal oxygens in the structure. True if terminal oxo is present in the entry, FAILED REFCODE if the refcode in the filename wasn't found in the CSD
 - **OH_removed** - the code removes OH bound to the metals because statistically they are most often water with missing hydrogens.
If it does so, it flags the molecule so you can check If it is actually OH or water with missing H
 - **oxo_OH** - identifies terminal oxo atoms that likely can be OH or water on the atoms that can have those issues (U and Zr)