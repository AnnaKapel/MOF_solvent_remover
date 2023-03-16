from ccdc.io import EntryReader
import os


def get_rAON_atomlabels(rAON_old):
    """Creating a dictionary of rAON, where each value corresponds to atom
    label instead of an atom object"""
    rAON_atomlabels = {}
    for key, value in rAON_old.items():
        label = key.label
        rAON_atomlabels[label] = value
    return rAON_atomlabels


def remove_free_solvent(molecule, rAON_atomlabels):  # and counterions
    """Function cleans the molecule from free solvent.
    The solvent is considered free if the initial fragment of a molecule
    object doesn't contain a metal"""
    structure = molecule.copy()

    # setting the variables for collecting statistics
    free_solvents = []
    counterions = []
    free_solvents_output = []
    counterions_output = []
    charge_removed = 0
    free_solvent_flag = False
    counterions_flag = False
    metal_counterion_flag = "."
    huge_counterion_flag = "."

    # checking if the component of structure contains metals
    for component in structure.components:
        flag_metal = False
        for atom in component.atoms:
            if atom.is_metal is True:
                flag_metal = True
                break
        labels = []
        for atom in component.atoms:
            labels.append(atom.atomic_symbol)

        # checking for metal-containing counterions and free solvents
        if flag_metal is True:
            if component.is_polymeric is False:
                # removing the metal-containing counterions and writing down their charge
                metal_counterion_flag = "TRUE"
                counterion_unit = []
                charge_unit = get_component_charge(rAON_atomlabels, component)

                # checking if it is one of the irrationally charged huge counterions
                if charge_unit > 10:
                    huge_counterion_flag = True

                charge_removed += charge_unit
                for atom in component.atoms:
                    counterions.append(atom.label)
                    counterion_unit.append(atom.label)
                counterions_output.append(counterion_unit)
        else:
            free_output_unit = []
            counterion_unit = []
            charge_unit = get_component_charge(rAON_atomlabels, component)

            # checking if the molecule is neutral - free solvent
            if charge_unit == 0:
                for atom in component.atoms:
                    free_solvents.append(atom.label)
                    free_output_unit.append(atom.label)
                free_solvents_output.append(free_output_unit)
            else:
                # hard coding ignoring the units that often appear as disorder as
                # neutral units
                neutral_units = [
                    ["O"],
                    ["N"],
                    ["H"],
                    ["O", "O"],
                    ["N", "N"],
                    ["O", "O", "O"],
                    ["D"],
                ]
                flag_neutral = False

                if len(labels) <= 3:
                    for unit in neutral_units:
                        if labels == unit:
                            flag_neutral = True
                            break
                    if len(labels) == 2:
                        if "O" in labels and "H" in labels:
                            flag_neutral = True
                        elif "O" in labels and "C" in labels:
                            flag_neutral = True
                    elif len(labels) == 3:
                        labels.sort()
                        if labels == ["C", "O", "O"]:
                            flag_neutral = True

                # writing the fragment to the corresponding list
                if flag_neutral is False:
                    charge_removed += charge_unit
                    for atom in component.atoms:
                        counterions.append(atom.label)
                        counterion_unit.append(atom.label)
                    counterions_output.append(counterion_unit)
                else:
                    for atom in component.atoms:
                        free_solvents.append(atom.label)
                        free_output_unit.append(atom.label)
                    free_solvents_output.append(free_output_unit)

    # getting the information for the outputs and statistics

    if len(free_solvents) > 0:
        free_solvent_flag = True
    if len(counterions) > 0:
        counterions_flag = True

    # removing the free solvent from the molecule object for it not to be
    # in the way for the future work
    for atom in structure.atoms:
        if atom.label in free_solvents or atom.label in counterions:
            structure.remove_atom(atom)

    # setting dictionary for storing the statistics
    statistics_output = {
        "free_solvents_output": free_solvents_output,
        "counterions_output": counterions_output,
        "charge_removed": charge_removed,
        "metal_counterion_flag": metal_counterion_flag,
        "huge_counterion_flag": huge_counterion_flag,
        "free_solvent_flag": free_solvent_flag,
        "counterions_flag": counterions_flag,
    }

    return structure, free_solvents, counterions, statistics_output


def get_component_charge(rAON, component):
    """Takes out the charge of the individual atoms by searching for their
    labels in the rAON dictionary

    Returns an int or a float - atom charge"""
    charge = 0
    for atom in component.atoms:
        key = atom.label
        charge_unit = rAON.get(key)
        charge += charge_unit
    return charge

def get_refcode(file):
            """Gets the CSD entry refcode. Works with filenames type AAAAAA.cif
            and AAAAA_xxx.cif

            Retuns: refcode as a str"""

            ref_code = file.replace(".cif", "")
            if "\\" or "/" in ref_code:
                ref_code = os.path.basename(ref_code)
            if ref_code.isupper() is True and "_" not in ref_code:
                return ref_code
            if "_" in ref_code:
                pos = ref_code.find("_")
                ref_code = ref_code[:pos]
                return ref_code

def get_oxo(molecule, file):
    """Getting the terminal oxygens from the molecule that are supposed to
    be removed as water. Filters out the real oxo and flags if there are
    any terminal oxygens present.

    Returns: list of atom labels of oxygens to remove"""

    def check_entry(file):
        """Checks the name of the MOF in the CSD entry if it contains something
        related to -oxo-

        Retuns: True if the oxo is present in the entry
        False if not"""

        global entry_oxo
        entry_oxo = False

        # getting MOF information from CSD
        ref_code = get_refcode(file)
        csd_reader = EntryReader("CSD")
        mof = csd_reader.entry(ref_code)

        # checking if terminal oxygens have to stay
        name = mof.chemical_name
        if "oxo" or "oxa" in name:
            oxo_names = [
                "-oxo-",
                "-dioxo-",
                "-trioxo-",
                "-tetraoxo-",
                "-pentaoxo-",
                "-hexaoxo-",
                "-heptaoxo-",
                "-octaoxo-",
                "-nonaoxo-",
                "-decaoxo-",
                "-undecaoxo-",
                "-undecaoxa-",
                "-dodecaoxo-",
                "-bis(oxo-",
                "-icosaoxo-",
                "-tetracosaoxo-",
            ]
            for oxo_name in oxo_names:
                if oxo_name in name:
                    entry_oxo = True

        if entry_oxo is True:
            return True
        else:
            return False

    terminal_oxo_flag = False
    oxo_OH = False

    # list of the most probable metals to have oxo on them
    most_probable_oxo = ["W", "U", "Mo", "V", "Np", "Ti"]
    terminal_oxo = []
    oxo_present = check_entry(file)
    possible_oxo = []
    corresponding_metals = []
    # writes to list of possible oxo atoms and the metals that correspond to them
    # the list contains atom objects
    for atom in molecule:
        if atom.atomic_symbol == "O":
            if len(atom.neighbours) == 1:
                bond = atom.bonds
                if bond[0].bond_type == "Single":
                    probable_metal = atom.neighbours
                    if probable_metal[0].is_metal is True:
                        possible_oxo.append(atom)
                        corresponding_metals.append(probable_metal[0])

                        # some atoms are probable to have both oxo and OH or water
                        if probable_metal[0].atomic_symbol == ("Zr" or "U"):
                            oxo_OH = True

    """If the entry says that the oxos are present:
    The metals list is checked for the presence of most probable ones.
    If the most probable ones are present, the oxos on them are kept and all
    the other ones are removed.
    If there are no most probable metals, all the oxos are kept.

    If the entry says that the oxos are not present:
    All the oxos are removed"""

    # skipping further steps if no oxo
    if len(possible_oxo) != 0:
        # checking for presence of most probable metals
        flag_probable = False
        if oxo_present:
            for metal in corresponding_metals:
                if metal.atomic_symbol in most_probable_oxo:
                    flag_probable = True
                    break
            if flag_probable is True:
                for index, metal in enumerate(corresponding_metals):
                    if metal.atomic_symbol not in most_probable_oxo:
                        terminal_oxo.append(possible_oxo[index].label)
            else:
                for atom in possible_oxo:
                    terminal_oxo.append(atom.label)
        else:
            for atom in possible_oxo:
                terminal_oxo.append(atom.label)

        if len(terminal_oxo) > 0:
            terminal_oxo_flag = True

    statistics_output = {
        "entry_oxo": entry_oxo,
        "terminal_oxo_flag": terminal_oxo_flag,
        "oxo_OH": oxo_OH,
        "oxo_mols": terminal_oxo,
    }

    return terminal_oxo, statistics_output


def remove_metals(work_mol):
    """Removes the metals from the molecule object by cheking the atoms
    Output: molecule object with no metals"""
    for atom in work_mol.atoms:
        if atom.is_metal is True:
            work_mol.remove_atom(atom)
    return work_mol


def define_solvents(mol_no_metals, charges):
    """Selecting initial list of bound solvents. If the fragment is not
    charged it is considered a solvent.
    OH on the metals is usually water with missing H, so it's written is a
    bound solvent as well.

    Returns: a list of the atom labels that have to be removed"""

    OH_removed = "."

    mol_fragments = []
    mol_fragments_list = mol_no_metals.components

    # removing extra atoms selected by CSD API
    for component in mol_fragments_list:
        if len(component.atoms) > 1:
            mol_fragments.append(component)
    # collecting fragments with zero charge
    not_charged = []
    for fragment in mol_fragments:
        charge = 0
        for atom in fragment.atoms:
            key = atom.label
            charge_unit = charges.get(key)
            charge += charge_unit
        if charge == 0:
            not_charged.append(fragment)
        else:
            # checking for OH that is usually water
            if len(fragment.atoms) == 2:
                atoms = fragment.atoms
                atoms_labels = [atom.atomic_symbol for atom in atoms]
                if "O" in atoms_labels and "H" in atoms_labels:
                    OH_removed = True
                    not_charged.append(fragment)

    statistics_output = {"OH_removed": OH_removed}

    return not_charged, statistics_output


def check_solvent(solvents, binding_sites_labels, molecule, unique_sites):
    """This function checks the solvent molecule object if they satisfy
    the conditions:
    1. Not a carbonyl
    2. Bound to only one metal

    Also it composes solvent lists for the output - 2D array.

    Returns: a list of atom labels to remove.
    """

    solvents_final = []
    solvent_type_label = {}

    # setting variables for outputs to write to dict
    solvents_for_output = []
    flag_aromatic = "."
    flag_double = "."
    bound_solvent_flag = False

    # writing a list of bridging oxo to not remove bridging neutral_units
    bridging_oxo = []
    for atom in unique_sites:
        if atom.atomic_symbol == "O":
            neighbours = atom.neighbours
            n_count = 0
            for unit in neighbours:
                if unit.is_metal is True:
                    n_count += 1
            if n_count > 1:
                bridging_oxo.append(atom.label)

    for solvent in solvents:
        solvent_labels = []
        # checking if there are carbonyls present and if they are -
        # skipping the carbonyls
        if len(solvent.atoms) == 2:
            symbols = [atom.atomic_symbol for atom in solvent.atoms]
            if "C" in symbols and "O" in symbols:
                continue

        # making a dictionary atom object-label to use for flags
        atoms_list = solvent.atoms
        for i in range(len(atoms_list)):
            label = atoms_list[i].label
            solvent_labels.append(label)
            solvent_type_label[label] = atoms_list[i]
        # checking if the solvent is bound to only one metal
        count = 0
        for i in range(len(solvent_labels)):
            if solvent_labels[i] in binding_sites_labels:
                count += 1

        if count == 1:
            # checking for bridging stuff
            flag_bridging = False
            for atom in solvent_labels:
                if atom in bridging_oxo:
                    flag_bridging = True

            if flag_bridging is False:
                solvents_for_output.append(solvent_labels)

                # wriitng atom labels as list for solvent removal
                for label in solvent_labels:
                    solvents_final.append(label)

        # getting flags for the outputs
    for solvent in solvents_for_output:
        for atom in solvent:
            if atom in binding_sites_labels:
                atom_obj = solvent_type_label.get(atom)
                for bond in atom_obj.bonds:
                    if bond.bond_type == "Aromatic":
                        flag_aromatic = True
                        break
                    if bond.bond_type == "Double":
                        flag_double = True
                        break
        if flag_double == "True" and flag_aromatic == "True":
            break

    # checking presence of the bound solvent
    if len(solvents_final) != 0:
        bound_solvent_flag = True

    statistics_output = {
        "solvents_for_output": solvents_for_output,
        "flag_aromatic": flag_aromatic,
        "flag_double": flag_double,
        "bound_solvent_flag": bound_solvent_flag,
    }

    return solvents_final, statistics_output


def get_solvents_to_remove(solvent_mols, free_mols, counterion_mols, oxo_mols):
    """Puts together an overall list of solvents to remove from the text file.

    Retuns: list of atom labels"""

    final_solvents = []
    final_solvents.extend(solvent_mols)
    final_solvents.extend(free_mols)
    final_solvents.extend(counterion_mols)
    final_solvents.extend(oxo_mols)
    return final_solvents
