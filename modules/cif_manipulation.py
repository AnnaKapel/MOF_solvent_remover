from ccdc import io
from ccdc import crystal

def readentry (input_cif):
    """reads the input cif into a CCDC.Crystal type object,
    and converts to standard atom labelling convention
    (necessary for the rest of the code to function properly"""

    #read in the cif to a crystal object
    with io.CrystalReader(input_cif, format='cif') as readcif:
        cif = readcif[0]
        cif.assign_bonds()
    readcif.close()

    # to remove duplicate atoms, need the empirical formula
    formula = cif.formula
    elamnt = formula.split(' ')

    #now convert to standard labelling convention and identify
    #duplicate atoms to be removed
    with open(input_cif, 'r') as file:
        file.seek(0)
        newstring = str()
        lines = file.readlines()
        loop_pos = 0
        start = 0
        end = 0
        columncount = 0

        for i, line in enumerate(lines):
            #locate atom type and site label columns
            if 'loop_' in line:
                loop_pos = i
            if ("_atom" in line) and (not "_geom" in line) and (not "_aniso" in line):
                start = loop_pos + 1
                end = i + 1
        for i in range(start, end):
            if 'atom_site_type_symbol' in lines[i]:
                type_pos = columncount
            if 'atom_site_label' in lines[i]:
                label_pos = columncount
            columncount += 1
        # keep track of relative position of type and label
        counting = {}
        cutoff = {}
        # to_remove = []
        for i in range(end, len(lines)):
            if 'loop_' in lines[i]:
                break
            #lines with atom information will contain a ., so only look at these
            if '.' in lines[i]:
                # split lines by whitespace
                col = lines[i].split()
                #keep count of how many of each element type
                if not col[type_pos] in counting:
                    counting[col[type_pos]] = 1
                elif col[type_pos] in counting:
                    counting[col[type_pos]] +=1
                #new atom labels
                newlabel = f'{col[type_pos]}{counting[col[type_pos]]}'
                lines[i] = lines[i].replace(col[label_pos], newlabel)
                #cutoff repeated atoms
                if newlabel in elamnt:
                    cutoff[col[type_pos]] = counting[col[type_pos]]
                # if col[type_pos] in cutoff:
                #     if counting[col[type_pos]] > cutoff[col[type_pos]]:
                #         to_remove.append(lines[i])

        # print('To remove:', to_remove)
        #remove unnecessary atoms
        # for i in to_remove:
        #     lines.remove(i)
        #     print(i)

        #combine to new string
        for i in lines:
            newstring += i
        #read into new crystal object and assign bonds
        newcif = crystal.Crystal.from_string(newstring, format='cif')
        newcif.assign_bonds()
        file.close()
    return (newcif)

def get_coordinates(molecule, solvents):
    '''Gets the coordinate lists of the individual atoms for removal

    Returns: a list of lits of three atom coordinates as str with 5 decimal places'''
    atoms_coordinates = []
    for atom in molecule.atoms:
        if atom.label in solvents:
            coord_unit = []
            coord_x = coord_unit.append(atom.fractional_coordinates.x)
            coord_y = coord_unit.append(atom.fractional_coordinates.y)
            coord_z = coord_unit.append(atom.fractional_coordinates.z)
            #firmatting the coordinates to match the cif format
            for i in range(len(coord_unit)):
                if 0 < coord_unit[i] < 1:
                    coord_unit[i] = format(coord_unit[i], '.5f')
                elif coord_unit[i] == -0.0:
                    coord_unit[i] = '0.00000'
                elif coord_unit[i] < 0:
                    if coord_unit[i] >= -1:
                        coord_unit[i] += 1
                    elif coord_unit[i] < -1:
                        coord_unit[i] +=2
                    coord_unit[i] = format(coord_unit[i], '.5f')
                elif coord_unit[i] >= 1:
                    if coord_unit[i] >= 2:
                        coord_unit[i] -= 2
                    else:
                        coord_unit[i] -= 1
                    coord_unit[i] = format(coord_unit[i], '.5f')
            #writing the individual coordinate list to the total list of coordinates
            atoms_coordinates.append(coord_unit)
    return atoms_coordinates

def write_output_cif():
    pass

def remove_solvents_from_file(lines, coordinates):
    '''Goes through the cif file and removes the lines that correspond
    to the atoms that shoulb be removed. Works as a text parcer coparing strings

    Returns: a list of strings that are the lines that should be kept in the file'''
    atom_labels = []
    global atom_count
    atom_count = 0

    #checking for the pymatgen formatting - different cif files
    pymatgen_flag = False
    pymatgen = lines[0]
    content = pymatgen.split(' ')
    content_fixed = [x for x in content if x]

    #setting the positions of the needed elements in the file
    if content_fixed == ['#', 'generated', 'using', 'pymatgen\n']:
        pymatgen_flag = True
        label_position = 1
        len_content = 6
        x = 3
        y = 4
        z = 5
    else:
        label_position = 0
        len_content = 5
        x = 2
        y = 3
        z = 4

    #going through the list of coordinates and check each line for the presense
    #of the set of three coordinates
    for coord in coordinates:
        for line in lines:
            content = line.split(' ')
            content_fixed = [x for x in content if x]
            if content_fixed[0] == '#':
                continue
            check = False

            #ignoring the lines that clearly are not coordinates
            if len(content_fixed) < len_content:
                continue

            #making the coodinares strings fit a specific pattern for comparison
            for i in range(x, z + 1):
                if content_fixed[i][0] == '1':
                    content_fixed[i] = '0' + content_fixed[i][1:]
                elif content_fixed[i][0] == '2':
                    content_fixed[i] = '1' + content_fixed[i][1:]
                elif content_fixed[i] == '-0.00000':
                    content_fixed[i] = '0.00000'
                elif content_fixed[i][0] == '-':
                    num = float(content_fixed[i])
                    if num >= -1:
                        num += 1
                    elif num >= -2:
                        num += 2
                    content_fixed[i] = format(num, '.5f')
                if len(content_fixed[i]) != 7 and content_fixed[i][0].isnumeric():
                        num = float(content_fixed[i])
                        content_fixed[i] = format(num, '.5f')
                if content_fixed[i][0].isnumeric() and '\n' in content_fixed[i]:
                    content_fixed[i] = content_fixed[i].strip()
                    num = float(content_fixed[i])
                    content_fixed[i] = format(num, '.5f')

            #comparing the coordinates to the line
            if content_fixed[x][:-1] == str(coord[0][:-1]):
                if content_fixed[y][:-1] == str(coord[1][:-1]):
                    if content_fixed[z][:-1] == str(coord[2][:-1]):
                        check = True
            #if the coordinates match the line it gets removed
            #atom label is extracted from such line
            if check == True:
                label = content_fixed[label_position]
                atom_labels.append(label)
                lines.remove(line)
                atom_count += 1

    #going through the list of atom labels to remove all the lines that
    #correspond to these atom labels
    for atom in atom_labels:
        for line in lines:
            content = line.split(' ')
            content_fixed = [x for x in content if x]
            if atom in content_fixed:
                lines.remove(line)
    return lines