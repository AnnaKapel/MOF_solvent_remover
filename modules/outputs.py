import pandas as pd
import os
import csv
import random
import shutil
from .cif_manipulation import remove_solvents_from_file

# output_solvent_stats = {'free_solvents_output': [['O13', 'H7', 'H9'], ['O14', 'H8', 'H10']], 
#                         'counterions_output': [], 
#                         'charge_removed': 0, 
#                         'metal_counterion_flag': '.', 
#                         'huge_counterion_flag': '.', 
#                         'free_solvent_flag': True, 
#                         'counterions_flag': False, 
#                         'entry_oxo': False, 
#                         'terminal_oxo_flag': False, 
#                         'oxo_OH': False, 
#                         'OH_removed': '.', 
#                         'solvents_for_output': [['O11', 'H3', 'H5'], ['O12', 'H4', 'H6']], 
#                         'flag_aromatic': '.', 
#                         'flag_double': '.', 
#                         'bound_solvent_flag': True,
#                         'oxo_mols': []}

# solvent_present_flag = True

def command_line_output(solvent_stats, solvent_present_flag):
    '''Printing information to the terminal'''

    if solvent_present_flag:
        if solvent_stats.get('free_solvent_flag'):
            free_solvent = solvent_stats.get('free_solvents_output')
            print(f'Identified free solvent:\n{free_solvent}')
            print('-----')
            print(f'Number of molecules: {len(free_solvent)}')
            print('-----')
        if solvent_stats.get('counterions_flag'):
            counterions = solvent_stats.get('counterions_output')
            print(f'Identified counterions:\n{counterions}')
            print('-----')
            print(f'Number of molecules: {len(counterions)}')
            print('-----')
        if solvent_stats.get('bound_solvent_flag'):
            bound_solvent = solvent_stats.get('solvents_for_output')
            print(f'Identified bound solvent molecules:\n{bound_solvent}')
            print('-----')
            print(f'Number of molecules: {len(bound_solvent)}')
            print('-----')   
        print('*********************')       
    else:
        print('No solvent or counterions identified')
        print('*********************')

def output_csv(solvent_stats, solvent_present_flag, total_solv_atoms, folder_path, file):
    #-----------------------
    #TO DO: atoms removed number from parcer
    #atoms match flag
    #-----------------------
    
    #path to the output csv file
    path_to_csv = os.path.join(folder_path,'Solvent_removal_results.csv')

    terminal_oxo_flag = 1
    entry_oxo = 1
    oxo_OH =1 
    #forming output row depending on presence of solvent
    if solvent_present_flag:
        #getting items from storage dict
        solvents_for_output = solvent_stats.get('solvents_for_output')
        free_solvents_output = solvent_stats.get('free_solvents_output')
        counterions_output = solvent_stats.get('counterions_output')
        terminal_oxo = solvent_stats.get('oxo_mols')

        #forming output row to append to the csv
        output_row = [file, 'YES', solvents_for_output, len(solvents_for_output),
                    free_solvents_output, len(free_solvents_output), 
                    counterions_output, len(counterions_output),
                    terminal_oxo, len(terminal_oxo), solvent_stats.get('charge_removed'),
                    total_solv_atoms, 'atoms removed', 'atoms_match_flag',
                    solvent_stats.get('flag_double'), solvent_stats.get('flag_aromatic'),
                    solvent_stats.get('metal_counterion_flag'), solvent_stats.get('terminal_oxo_flag'),
                    solvent_stats.get('entry_oxo'), solvent_stats.get('huge_counterion_flag'),
                    solvent_stats.get('OH_removed'), solvent_stats.get('oxo_OH')] 
    
    else:
        #getting the output items from the storage dictionary
        terminal_oxo_flag = solvent_stats.get('terminal_oxo_flag')
        entry_oxo = solvent_stats.get('entry_oxo')
        oxo_OH = solvent_stats.get('oxo_OH')

        #forming output row to append to the csv
        output_row = [file, 'NO', '.', '.', '.','.','.','.','.','.','.',
                    '.','.','.','.','.','.', terminal_oxo_flag, entry_oxo,
                    '.', '.', oxo_OH]
    
    return output_row
    
    #if the file is not present -> creates the file and adds header
    # if os.path.exists(path_to_csv):
    #     with open(path_to_csv, 'a',  newline='') as fileObj:
    #         writerObj = csv.writer(fileObj)
    #         writerObj.writerow(output_row)
        
    # else:
    #     with open(path_to_csv, 'w',  newline='') as fileObj:
    #         writerObj = csv.writer(fileObj)
    #         writerObj.writerow(['CIF', 'Solvent', 'Bound solvent', 'Number of bound mols', 
    #                             'Free solvent', 'Number of free solvent molecules',
    #                             'Counterions', 'Number of counterions', 'Terminal oxo',
    #                             'Number of terminal oxo', 'Charge_removed', 'Total atoms',
    #                             'Atoms removed', 'atoms_match_flag', 'flag_double',
    #                             'flag_aromatic', 'Metal counterion flag', 'terminal_oxo_flag',
    #                             'Entry terminal oxo', 'huge_counterion_flag', 'OH_removed',
    #                             'oxo_OH'])
    #         writerObj.writerow(output_row)

def output_cif(path, file, solvent_coordinates, cwd):
    output_directory = os.path.join(path, 'MOFs_removed_solvent')
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    
    new_cif_file = shutil.copy2(os.path.join(cwd, file) , output_directory)
    with open(new_cif_file, 'r+') as cif_file:
        lines = cif_file.readlines()
        cif_file.seek(0)
        cif_file.truncate()
        file_content = remove_solvents_from_file(lines, solvent_coordinates)
        for line in file_content:
            cif_file.write(line)
        cif_file.close()


