import os
import shutil
import csv

from .cif_manipulation import remove_solvents_from_file


def command_line_output(solvent_stats, solvent_present_flag):
    """Printing information to the terminal"""

    if solvent_present_flag:
        if solvent_stats.get("free_solvent_flag"):
            free_solvent = solvent_stats.get("free_solvents_output")
            print(f"Identified free solvent:\n{free_solvent}")
            print("-----")
            print(f"Number of molecules: {len(free_solvent)}")
            print("-----")
        if solvent_stats.get("counterions_flag"):
            counterions = solvent_stats.get("counterions_output")
            print(f"Identified counterions:\n{counterions}")
            print("-----")
            print(f"Number of molecules: {len(counterions)}")
            print("-----")
        if solvent_stats.get("bound_solvent_flag"):
            bound_solvent = solvent_stats.get("solvents_for_output")
            print(f"Identified bound solvent molecules:\n{bound_solvent}")
            print("-----")
            print(f"Number of molecules: {len(bound_solvent)}")
            print("-----")
        print("*********************")
    else:
        print("No solvent or counterions identified")
        print("*********************")


def output_csv(
    solvent_stats, solvent_present_flag, total_solv_atoms, file, removed_atoms
):
    file = os.path.basename(file)

    # forming output row depending on presence of solvent
    if solvent_present_flag:
        # getting items from storage dict
        solvents_for_output = solvent_stats.get("solvents_for_output")
        free_solvents_output = solvent_stats.get("free_solvents_output")
        counterions_output = solvent_stats.get("counterions_output")
        terminal_oxo = solvent_stats.get("oxo_mols")

        atom_count_match = True
        if total_solv_atoms != removed_atoms:
            atom_count_match = False

        # forming output row to append to the csv
        output_row = [
            file,
            "YES",
            solvents_for_output,
            len(solvents_for_output),
            free_solvents_output,
            len(free_solvents_output),
            counterions_output,
            len(counterions_output),
            terminal_oxo,
            len(terminal_oxo),
            solvent_stats.get("charge_removed"),
            total_solv_atoms,
            removed_atoms,
            atom_count_match,
            solvent_stats.get("flag_double"),
            solvent_stats.get("flag_aromatic"),
            solvent_stats.get("metal_counterion_flag"),
            solvent_stats.get("terminal_oxo_flag"),
            solvent_stats.get("entry_oxo"),
            solvent_stats.get("huge_counterion_flag"),
            solvent_stats.get("OH_removed"),
            solvent_stats.get("oxo_OH"),
        ]

    else:
        # getting the output items from the storage dictionary
        terminal_oxo_flag = solvent_stats.get("terminal_oxo_flag")
        entry_oxo = solvent_stats.get("entry_oxo")
        oxo_OH = solvent_stats.get("oxo_OH")

        # forming output row to append to the csv
        output_row = [
            file,
            "NO",
            ".",
            ".",
            ".",
            ".",
            ".",
            ".",
            ".",
            ".",
            ".",
            ".",
            ".",
            ".",
            ".",
            ".",
            ".",
            terminal_oxo_flag,
            entry_oxo,
            ".",
            ".",
            oxo_OH,
        ]

    return output_row

def output_cif(path, file, solvent_coordinates, cwd):
    output_directory = os.path.join(path, "MOFs_removed_solvent")
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    new_cif_file = shutil.copy2(os.path.join(cwd, file), output_directory)
    with open(new_cif_file, "r+") as cif_file:
        lines = cif_file.readlines()
        cif_file.seek(0)
        cif_file.truncate()
        file_content, removed_atoms = remove_solvents_from_file(lines, solvent_coordinates)
        for line in file_content:
            cif_file.write(line)
        cif_file.close()

    return removed_atoms


def remove_solvent(mol, solvents_to_remove):
    for atom in mol.atoms:
        if atom.label in solvents_to_remove:
            mol.remove_atom(atom)
    return mol