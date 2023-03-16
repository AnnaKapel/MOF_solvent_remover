#!/usr/bin/env python3
import argparse
import glob
import os
import time
from multiprocessing import Pool
import ccdc
from ccdc import io
from ccdc.io import CrystalWriter

from modules.assign_charges import *
from modules.cif_manipulation import *
from modules.outputs import *
from modules.solvent_analysis import *

parser = argparse.ArgumentParser()
parser.add_argument(
    "--files_path", default=os.getcwd(), help="specify the folder where MOFs are stored"
)
parser.add_argument(
    "--export_path", default=os.getcwd(), help="specify the output folder"
)
parser.add_argument(
    "--n_processes", default=4, help="specify number of parallel processes"
)
parser.add_argument(
    "-v", "--verbose", action="store_true", help="if set, turns off outputs"
)
args = parser.parse_args()

cwd = args.files_path
output_dir = args.export_path

def worker(file):
    try:
        start_time_MOF = time.time()
        print(f"Analyzing {file}...")
        input_cif = file
        # read in the .cif, extract the underlying molecule,
        # identify the unique sites, metal sites, binding sites

        cif = readentry(input_cif)
        mol = cif.molecule
        asymmol = cif.asymmetric_unit_molecule

        uniquesites = get_unique_sites(mol, asymmol)
        metalsites = get_metal_sites(uniquesites)
        binding_sites = get_binding_sites(
            metalsites, uniquesites
        )  # outputs list of atoms bound to metals

        # Now get the localized oxidation state contribution of each atom
        # need delocalized bond contributions
        dVBO = delocalisedLBO(mol)
        # then need ring bond contributions - NEW TO dev_194
        rVBO = ringVBOs(mol)
        # finally combine delocal/aromatic bond conrtibutions with localized bonding
        AON = iVBS_Oxidation_Contrib(uniquesites, rVBO, dVBO)
        # Previous only assigns an oxidation contribution to unique images of atoms,
        # also need to assign these values to redundant sites:
        rAON = redundantAON(AON, mol)

        # the atoms in rAON calculated and atoms in molecule_no_metals are different
        # redoing the rAON dictionary to have atom labels instead of atoms
        rAON_atomlabels = get_rAON_atomlabels(rAON)

        # creates a copy of the initial molecule
        # cleans the molecule from free solvents
        (
            molecule_work,
            free_solvents,
            counterions,
            solvent_stats_batch1,
        ) = remove_free_solvent(mol, rAON_atomlabels)

        # identifying the oxo molecules
        oxo_mols, solvent_stats_batch2 = get_oxo(uniquesites, file)

        # removes metals from the molecule
        molecule_no_metals = remove_metals(molecule_work)

        # creates the initial list of possible solvents
        solvent_mols, solvent_stats_batch3 = define_solvents(
            molecule_no_metals, rAON_atomlabels
        )

        # rewriting binding sites as labels
        binding_sites_labels = [site.label for site in binding_sites]

        # check the possible solvents
        solvent_mols_checked, solvent_stats_batch4 = check_solvent(
            solvent_mols, binding_sites_labels, mol, uniquesites
        )

        # combine all the lists of items to remove
        solvents_to_remove = get_solvents_to_remove(
            solvent_mols_checked, free_solvents, counterions, oxo_mols
        )

        # checking if there is any solvent
        solvent_present_flag = False
        if len(solvents_to_remove) > 0:
            solvent_present_flag = True

        # getting the coordinates of the atoms for removal
        solvent_coordinates = get_coordinates(mol, solvents_to_remove)
        
        ######## test ##########
        # mol_no_solvent = remove_solvent(mol, solvents_to_remove)
        # with CrystalWriter('test_cif.cif') as filehandle:
        #     filehandle.write(mol)

        # putting together dictionary of otput stats
        output_solvent_stats = {
            **solvent_stats_batch1,
            **solvent_stats_batch2,
            **solvent_stats_batch3,
            **solvent_stats_batch4,
        }

        total_solv_atoms = len(solvents_to_remove)

        print("--- %s seconds ---" % (time.time() - start_time_MOF))

        # removing solvent from cif file if solvent is present
        if solvent_present_flag:
            removed_atoms = output_cif(output_dir, file, solvent_coordinates, cwd)

        # output a csv with solvent removal stats
        output_row = output_csv(
            output_solvent_stats, solvent_present_flag, total_solv_atoms, file, removed_atoms
        )

        if not args.verbose:
            command_line_output(output_solvent_stats, solvent_present_flag)

        print("--- %s seconds ---" % (time.time() - start_time_MOF))

        return output_row
    
    except Exception as e:
        print(f'ERROR >>>>> {file}')


if __name__ == "__main__":
    start_time = time.time()

    export_df_path = os.path.join(args.export_path, "Solvent_removal_results.csv")
    failed_cifs_path = os.path.join(args.export_path, "Filed_cifs.txt")

    files = [os.path.join(cwd, file) for file in os.listdir(cwd) if file.endswith(".cif")]

    print(files)
    print(f"### {len(files)} .cif file(s) detected in {cwd} ###")

    pool = Pool(processes = args.n_processes)

    for res in pool.imap(worker, files):
        if res is not None:
            if os.path.exists(export_df_path):
                with open(export_df_path, 'a',  newline='') as file_obj:
                    writerObj = csv.writer(file_obj)
                    writerObj.writerow(res)        
            else:
                with open(export_df_path, 'w',  newline='') as fileObj:
                    writerObj = csv.writer(fileObj)
                    writerObj.writerow(['CIF', 'Solvent', 'Bound solvent', 'Number of bound mols', 
                                        'Free solvent', 'Number of free solvent molecules',
                                        'Counterions', 'Number of counterions', 'Terminal oxo',
                                        'Number of terminal oxo', 'Charge_removed', 'Total atoms',
                                        'Atoms removed', 'atoms_match_flag', 'flag_double',
                                        'flag_aromatic', 'Metal counterion flag', 'terminal_oxo_flag',
                                        'Entry terminal oxo', 'huge_counterion_flag', 'OH_removed',
                                        'oxo_OH'])
                    writerObj.writerow(res)

    print("--- %s seconds ---" % (time.time() - start_time))
