import os
import mendeleev
from ccdc import descriptors
from ccdc import molecule

def get_unique_sites(mole, asymmole):
    #blank list for unique sites
    uniquesites = []
    labels = []
    asymmcoords = []
    molecoords = []
    duplicates = []
    for atom in asymmole.atoms:
        asymmcoords.append(atom.coordinates) #writing coordinates of atoms in asymetric unit
    for atom in mole.atoms:
        if atom.coordinates in asymmcoords: #if they are coordinates of atom in asymetric unit
            if not atom.coordinates in molecoords: #and they are not in list of coordinates in whole MOF
                if not atom.label in labels: #and we didn't have this atom before
                    uniquesites.append(atom)
                    molecoords.append(atom.coordinates)
                    labels.append(atom.label)
                else:
                    duplicates.append(atom)
            else:
                duplicates.append(atom)
    if len(duplicates) >= 1:
        for datom in duplicates:
            for atom in uniquesites:
                if any([
                    (datom.coordinates == atom.coordinates),
                    (datom.label == atom.label)
                    ]):
                    if datom.atomic_symbol == atom.atomic_symbol:
                        if len(datom.neighbours) > len(atom.neighbours):
                            uniquesites.remove(atom)
                            uniquesites.append(datom)
                    elif not datom.label in labels:
                        uniquesites.append(datom)
    return uniquesites

def get_metal_sites(sites):
    metalsites = []
    for site in sites:
        if site.is_metal == True:
            metalsites.append(site)
    return metalsites

def get_binding_sites(metalsites, uniquesites):
    binding_sites = set()
    for metal in metalsites:
        for ligand in metal.neighbours:
            for site in uniquesites:
                if ligand.label == site.label:
                    binding_sites.add(site)
    return (binding_sites)

# def get_sphere_2(sphere_1, sites):
#     sphere2 = {}
#     for metal in sphere_1:
#         for ligand in sphere_1[metal]:
#             sphere = []
#             for environ in ligand.neighbours:
#                 for site in sites:
#                     if environ.label == site.label:
#                         sphere.append(site)
#             sphere2[ligand] = sphere
#     return sphere2

def ringVBOs(mole):
    ringVBO = {}
    unassigned = mole.atoms
    ringcopy = mole.copy()
    oncycle_atoms = []
    offcycle_atoms = []
    oncycle_labels = []
    offcycle_labels = []
    #remove all the metals, this 
    #prevents metal-containing rings (i.e. pores)
    #from interfering
    for atom in ringcopy.atoms:
        if atom.is_metal:
            ringcopy.remove_atom(atom)
    #collect all the cyclic atoms
    for atom in ringcopy.atoms:
            if atom.is_cyclic:
                if not atom in oncycle_atoms:
                    oncycle_atoms.append(atom)
                    oncycle_labels.append(atom.label)
    #we also need everything that the cyclic atoms are bound to
    for atom in oncycle_atoms:
        for neighbour in atom.neighbours:
            if not neighbour in oncycle_atoms:
                if not neighbour in offcycle_atoms:
                    offcycle_atoms.append(neighbour)
                    offcycle_labels.append(neighbour.label)
    cyclicsystem = (oncycle_atoms + offcycle_atoms)
    #remove every atom that isn't part of or directly bound to a cycle
    for atom in ringcopy.atoms:
        if not atom in cyclicsystem:
            ringcopy.remove_atom(atom)
    #find all non-cyclic bonds
    #bonds between cycles, break and cap with H
    for bond in ringcopy.bonds:
        if not bond.is_cyclic:
            #bonds between cycles
            if all((member.label in oncycle_labels for member in bond.atoms)):
                member1 = bond.atoms[0]
                member2 = bond.atoms[1]
                Hcap1 = molecule.Atom('H', coordinates = member1.coordinates)
                Hcap2 = molecule.Atom('H', coordinates = member2.coordinates)
                Hcap1_id = ringcopy.add_atom(Hcap1)
                Hcap2_id = ringcopy.add_atom(Hcap2)
                ringcopy.add_bond(bond.bond_type, Hcap1_id, member2)
                ringcopy.add_bond(bond.bond_type, Hcap2_id, member1)
                ringcopy.remove_bond(bond)
    #cap off-cycle atoms
    for offatom in offcycle_atoms:
        #get the VBO for each off-cycle atom
        #(VBO with respect to cyclic atoms)
        offVBO = 0
        #quick check for delocalized systems in the ring
        #if there are any, get the delocalised bond orders
        if any (bond.bond_type == 'Delocalised' for bond in offatom.bonds):
            offdVBO = delocalisedLBO(offcycle_atoms)
        for bond in offatom.bonds:
            #Each bond contributes to Ligand Bond Order according to its type
            if bond.bond_type == 'Single':
                offVBO +=1
            elif bond.bond_type == 'Double':
                offVBO +=2
            elif bond.bond_type == 'Triple':
                offVBO +=3
            elif bond.bond_type == 'Quadruple':
                offVBO +=4
            elif bond.bond_type == 'Delocalised':
                offVBO += offdVBO[offatom]
            elif bond.bond_type == 'Aromatic':
                offVBO += 0
                print ("impossible Aromatic bond")
        #cap with appropriate element for VBO
        if offVBO == 1:
            offatom.atomic_symbol = 'H'
        elif offVBO == 2:
            offatom.atomic_symbol = 'O'
        elif offVBO == 3:
            offatom.atomic_symbol = 'N'
        elif offVBO == 4:
            offatom.atomic_symbol = 'C'
        elif offVBO == 5:
            offatom.atomic_symbol = 'P'
        elif offVBO == 6:
            offatom.atomic_symbol = 'S'
        elif offVBO  > 6:
            print ("no, that's too many")
    #for each cyclic system, reassign bonds, kekulize, and get VBO
    #the bond and atom pruning we did above ensures that fused cycles
    #will be treated as a single system
    #while non-fused cycles that are connected via bonding are treated
    #as seperate systems
    for cyclesys in ringcopy.components:
        #reassign bonds and kekulize
        cyclesys.assign_bond_types()
        cyclesys.kekulize()
        #quick check for delocalized systems in the ring
        #if there are any, get the delocalised bond orders
        if any (bond.bond_type == 'Delocalised' for bond in cyclesys.bonds):
            rdVBO = delocalisedLBO(cyclesys)
        #assign VBO for each on-cycle atom
        for ratom in cyclesys.atoms:
            rVBO = 0
            if ratom.label in oncycle_labels:
                for rbond in ratom.bonds:
                #Each bond contributes to Ligand Bond Order according to its type
                    if rbond.bond_type == 'Single':
                        rVBO +=1
                    elif rbond.bond_type == 'Double':
                        rVBO +=2
                    elif rbond.bond_type == 'Triple':
                        rVBO +=3
                    elif rbond.bond_type == 'Quadruple':
                        rVBO +=4
                    elif rbond.bond_type == 'Delocalised':
                        rVBO += rdVBO[ratom]
                    elif rbond.bond_type == 'Aromatic':
                        rVBO += 0
                        print ("impossible Aromatic bond")
                #the VBOs are currently associated to atom objects
                #in molecule objects that we have modified
                #we need these to be associated to atom objects in
                #the parent (unmodified) molecule object
                for matom in unassigned:
                    if matom.label == ratom.label:
                        ringVBO[matom] = rVBO
                        unassigned.remove(matom)
    return(ringVBO)

def assign_VBS(atom, rVBO, dVBO):
    """This function will assign a Valence-Bond-Sum (VBS) to an atom. 
    Takes one CCDC atom object, list of binding sites, and the metal-free
    molecule object as inputs"""
    VBO = 0
    if atom.is_metal:
        return(0)
    if atom in rVBO:
        VBO = rVBO[atom]
    else:
        for bond in atom.bonds:
            if any(batom.is_metal
            for batom in bond.atoms):
                VBO += 0
            #Each bond contributes to Ligand Bond Order according to its type
            elif bond.bond_type == 'Single':
                VBO +=1
            elif bond.bond_type == 'Double':
                VBO +=2
            elif bond.bond_type == 'Triple':
                VBO +=3
            elif bond.bond_type == 'Quadruple':
                VBO +=4
            elif bond.bond_type == 'Delocalised':
                VBO += dVBO[atom]
            elif bond.bond_type == 'Aromatic':
                VBO += rVBO[atom]
    return(VBO)

def delocalisedLBO(molecule):
    """returns a dict of all atoms with
    delocalised bonds and their (delocalized-only) VBS"""

    def TerminusCounter(atomlist):
        """Counts the number of termini in the delocalized bond system.
        takes a list of all atoms and list of all bonds in the delocalized system as inputs."""
        NTerminus = 0
        for member in atomlist:
            connectivity = 0
            for bond in member.bonds:
                if bond.bond_type == 'Delocalised':
                    connectivity += 1
            if connectivity is 1:
                NTerminus += 1
        return(NTerminus)

    def delocal_crawl (atomlist):
        """returns a list of all atoms in a delocalised system,
        takes at least one delocalised bonding atom as input"""
        for delocatom in atomlist:
            for bond in delocatom.bonds:
                if bond.bond_type == 'Delocalised':
                    for member in bond.atoms:
                        if not member in atomlist:
                            atomlist.append(member)
                            return(delocal_crawl(atomlist))
        return(atomlist)

    delocal_dict = {}
    for atom in molecule.atoms:
        if all([
            (any(
                bond.bond_type == 'Delocalised'
                for bond in atom.bonds
                )),
                (not atom in delocal_dict)
            ]):
            delocal_dict[atom] = []
            delocal_system = delocal_crawl([atom])
            NTerminus = TerminusCounter(delocal_system)
            for datom in delocal_system:
                connectivity = 0
                delocLBO = 0
                for neighbour in datom.neighbours:
                    if neighbour in delocal_system:
                        connectivity += 1
                if connectivity == 1:
                    #terminus
                    delocLBO = (NTerminus+1)/NTerminus
                if connectivity > 1 :
                    #node
                    delocLBO = (connectivity+1)/connectivity
                delocal_dict[datom] = delocLBO
    return(delocal_dict)

# def iVBS_FormalCharge(atom):
#     """determines the formal charge on an atom
#     that is NOT part of an aromatic or delocalized
#     bond system"""
#     VBO = 0
#     if atom.is_metal:
#         return(VBO)
#     CN = 0
#     for neighbour in atom.neighbours:
#         if not neighbour.is_metal:
#             CN+=1
#     valence = valence_e(atom)
#     charge = 0
#     for bond in atom.bonds:
#         if any(batom.is_metal for batom in bond.atoms):
#             VBO += 0
#         #Each bond contributes to Ligand Bond Order according to its type
#         elif bond.bond_type == 'Single':
#             VBO +=1
#         elif bond.bond_type == 'Double':
#             VBO +=2
#         elif bond.bond_type == 'Triple':
#             VBO +=3
#         elif bond.bond_type == 'Quadruple':
#             VBO +=4
#     #need the unpaired electrons
#     unpaired_e = (4 - abs(4 - valence))
#     #expanded valences require special handling
#     if VBO <= (unpaired_e):
#         charge = VBO - unpaired_e
#     # Expanded (2e) valences:
#     elif (VBO > unpaired_e) and (VBO < valence):
#         diff = VBO - unpaired_e
#         if diff <= 2:
#             UPE = valence - unpaired_e - 2
#         elif diff <= 4:
#             UPE = valence - unpaired_e - 4
#         elif diff <= 6:
#             UPE = valence - unpaired_e - 6
#         elif diff <= 8:
#             UPE = valence - unpaired_e - 8
#         charge = valence - (VBO + UPE)
#     elif VBO >= (valence):
#         charge = valence - VBO
#     return(charge)

def get_CN(atom):
    """returns the coordination number of an atom"""
    CN = 0
    for neighbour in atom.neighbours:
        if not neighbour.is_metal:
            CN+=1
    return(CN)

def valence_e(elmnt):
    """returns number of valence electrons of an element"""
    atom = mendeleev.element(elmnt.atomic_symbol)
    if atom.block == 's':
        valence = atom.group_id
    elif atom.block == 'p':
        valence = (atom.group_id - 10)
    elif atom.block == 'd':
        valence = atom.group_id
    elif atom.block == 'f':
        if atom.atomic_number in range(56,72):
            valence = atom.atomic_number - 57 + 3
        elif atom.atomic_number in range(88,104):
            valence = atom.atomic_number - 89 + 3
        else:
            raise ValueError("valence_e() >> Unexpected f block element", atom)
    elif atom.group_id == 18:
        valence = 8 if atom.symbol != 'He' else 2
    else:
        raise ValueError("valence_e() >> Unexpected valence electrons", atom)
    return valence

def carbocation_check (atom):
    """geometry checker for carbocations/carbanions returns
    tetrahedral(anion) or trigonal(cation) depending on bond angles"""
    abc = []
    #get atom neighbours
    for neighbours in atom.neighbours:
        if not neighbours.is_metal:
            abc.append(neighbours)
    #get all three relevant bond angles
    angle1 = descriptors.MolecularDescriptors.atom_angle(abc[0], atom, abc[1])
    angle2 = descriptors.MolecularDescriptors.atom_angle(abc[0], atom, abc[2])
    angle3 = descriptors.MolecularDescriptors.atom_angle(abc[1], atom, abc[2])
    #average the angels
    AVGangle = abs(angle1 + angle2 + angle3)/3
    #take the difference between the averaged bond angles and
    #ideal trigonal planar/tetrahedral bond angles
    tet = abs(AVGangle - 109.5)
    trig = abs(AVGangle - 120)
    if tet < trig:
        return('tetrahedral')
    if trig < tet:
        return('trigonal')

def carbene_type(atom):
    """distinguishes between singlet and triplet carbenes,
    only input suspected carbene atoms (2-coordinate carbon II)"""
    #get alpha-atoms
    alpha = atom.neighbours
    alpha_type = []
    #get element symbols for alpha atoms
    for a in alpha:
        if not a.is_metal:
            alpha_type.append(a.atomic_symbol)
    # if any alpha atom is a heteroatom, return "singlet"
    # these are Fischer carbenes
    for a in alpha_type:
        if not any([(a == 'C'),
                    (a == 'H')]):
            return('singlet')
    # if the carbene C is in a heterocycle,
    # return "singlet"
    # there are Arduengo carbenes (NHCs, CAACs)
    if atom.is_cyclic == True:
        for ring in atom.rings:
            for species in ring.atoms:
                if not species.atomic_symbol == 'C':
                    return('singlet')
    # for all other carbenes, return "triplet"
    # these are Schrock carbenes
    return('triplet')

def iVBS_Oxidation_Contrib(unique_atoms, rVBO, dVBO):
    """determines the oxidation state contribution for all
    unique atoms in a MOF. Returns a dictionary of Atom:Oxidaiton_Contribution
    pairs. Takes the unique sites in the MOF (without metal), the MOF molecule
    object (without metal) and a list of metal-binding sites as inputs"""
    VBS = 0
    CN = 0
    valence = 0
    oxi_contrib = {}
    # for each unique atom
    for atom in unique_atoms:
        # assign valence-bond-sum
        VBS = assign_VBS(atom, rVBO, dVBO)
        #determine coordination number
        CN = get_CN(atom)
        #  determine number of valence electrons
        valence = valence_e(atom)
        # get number of unpaired electrons in the free element
        unpaired_e = (4 - abs(4 - valence))

        #  metals do not contribute:
        if  atom.is_metal:
            oxi_contrib[atom] = 0
        # Normal valences:
        elif VBS <= (unpaired_e):
            oxi_contrib[atom] = unpaired_e - VBS
        # Expanded (2e) valences:
        elif (VBS > unpaired_e) and (VBS < valence):
            diff = VBS - unpaired_e
            if diff <= 2:
                UPE = valence - unpaired_e - 2
            elif diff <= 4:
                UPE = valence - unpaired_e - 4
            elif diff <= 6:
                UPE = valence - unpaired_e - 6
            elif diff <= 8:
                UPE = valence - unpaired_e - 8
            oxi_contrib[atom] = VBS + UPE - valence
        elif VBS >= (valence):
            oxi_contrib[atom] = VBS - valence
        
        # need to check for 3-coordinate carbocations,
        # 3-coordinate carbanions, carbenes, and heavier
        # homologues (these are not immediately detectable)
        if any([
            (atom.atomic_symbol == 'C'),
            (atom.atomic_symbol == 'Si'),
            (atom.atomic_symbol == 'Ge'),
            (atom.atomic_symbol is 'Pb')]):
             if not atom in rVBO:
              # 3 coordinate and VBS 3 could be
              # carbanion or carbocation
              if VBS == 3 and CN == 3:
                  geom = carbocation_check(atom)
                  if geom == "trigonal":
                      oxi_contrib[atom] = -1
                  if geom == "tetrahedral":
                      oxi_contrib[atom] = 1
             # VBS 2 and 2 coordinate is carbene,
             # but singlet or triplet?
             if VBS == 2 and CN == 2:
               carbene = carbene_type(atom)
               if carbene == 'singlet':
                 oxi_contrib[atom] = 2
               if carbene == 'triplet':
                 oxi_contrib[atom] = 0

        # Nitro groups frequently have both N-O bonds assigned
        # as double bonds, giving incorrect VBS of 5
        # and oxidation contribution of -2
        #this block catches this and applies a fix
        if all([
            (atom.atomic_symbol == 'N'),
            (VBS == 5 and CN == 3),
        ]):
            N_sphere1 = atom.neighbours
            O_count = 0
            for neighbour in N_sphere1:
                if neighbour.atomic_symbol == 'O':
                    O_count += 1
            geom = carbocation_check(atom)
            if O_count == 2 and geom == "trigonal":
                oxi_contrib[atom] = 0        
    return (oxi_contrib)

def redundantAON(AON, molecule):
    redAON = {}
    for rsite1 in molecule.atoms:
        for usite1 in AON:
            redAON[usite1] = AON[usite1]
            if rsite1.label == usite1.label:
                redAON[rsite1] = AON[usite1]
    return(redAON)

# def KnownONs():
#     KONs = {}
#     with open(os.path.join(cwd,'KnownON.csv')) as ONs:
#         for ON in ONs.readlines():
#             ONlist = []
#             splitON = ON.split(',')
#             for split in splitON:
#                 split.replace(',','')
#             while('' in splitON):
#                 splitON.remove('')
#             while('\n' in splitON):
#                 splitON.remove('\n')
#             for i in range(1,len(splitON)):
#                 ONlist.append(splitON[i])
#             KONs[splitON[0]] = ONlist
#     ONs.close()
#     return(KONs)

# def IonizationEnergies():
#     KIEs = {}
#     IElist = []
#     with open(os.path.join(cwd,'Ionization_Energies.csv')) as IEs:
#         for IE in IEs.readlines():
#             splitIE = IE.split(',')
#             for split in splitIE:
#                 split.replace(',','')
#                 if(r'\n' in split):
#                     split.replace(r'\n', '')
#             while('' in splitIE):
#                 splitIE.remove('')
#             while(r'\n' in splitIE):
#                 splitIE.remove(r'\n')
#             IElist.append(splitIE)
#         for entry in IElist:
#             if entry[1] in KIEs:
#                 KIEs[(entry[1])].append(entry[2])
#             else:
#                 KIEs[entry[1]] = []
#                 KIEs[entry[1]].append(entry[2])
#     IEs.close()
#     return(KIEs)

# def HighestKnownONs():
#     HKONs = {}
#     with open(os.path.join(cwd,'KnownON.csv')) as ONs:
#         for ON in ONs.readlines():
#             highest = 0
#             splitON = ON.split(',')
#             for split in splitON:
#                 split.replace(',','')
#             while('' in splitON):
#                 splitON.remove('')
#             while('\n' in splitON):
#                 splitON.remove('\n')
#             for i in range(1,len(splitON)):
#                 if int(splitON[i]) >= highest:
#                     highest = int(splitON[i])
#             HKONs[splitON[0]] = int(highest)
#     ONs.close()
#     return(HKONs)

# def ONprobabilities():
#     ONP = {}
#     with open(os.path.join(cwd,'Oxidation_Probabilities.csv')) as ONPs:
#         for ON in ONPs.readlines():
#             ONPlist = []
#             splitONP = ON.split(',')
#             for split in splitONP:
#                 split.replace(',','')
#             while('' in splitONP):
#                 splitONP.remove('')
#             while('\n' in splitONP):
#                 splitONP.remove('\n')
#             for i in range(1,len(splitONP)):
#                 ONPlist.append(float(splitONP[i]))
#             ONP[splitONP[0]] = ONPlist
#     ONPs.close()
#     return(ONP)

# '''END of the MOSAEC functionality'''