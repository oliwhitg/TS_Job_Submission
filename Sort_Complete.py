import numpy as np
import os
import time

current_dir = os.getcwd()

for minimisation_number in [10000, 100000,510000,600000]:

    if os.path.isfile(current_dir+"/{:}_data_export.dat".format(minimisation_number)) == True:
        with open('{:}_data_export.dat'.format(minimisation_number), 'r') as output:
    
            run_data = np.genfromtxt(output, dtype=str)
            print(run_data)
            time.sleep(1)
            previous_results = True
    else:
        previous_results = False
    
    
    
    
    with open('{:}_data_export.dat'.format(minimisation_number), 'a+') as output:
    
        cwd = current_dir+"/Minimisations/"
    
        for rings in os.listdir(cwd):
            cwd = current_dir+"/Minimisations/"
 
            print(rings)
    
            if os.path.isdir(cwd + rings):
                cwd = current_dir + "/Minimisations/" + rings
                for temp_folder in os.listdir(cwd):
                    cwd = current_dir+"/Minimisations/"+rings
                    print('    ', temp_folder)
                    print(current_dir)
                    print(cwd)
                    print(os.listdir(cwd+"/"+temp_folder))
                    if os.path.isdir(cwd+"/"+temp_folder):
                        cwd = current_dir+"/Minimisations/"+rings+"/"+temp_folder
    
                        T = int(temp_folder[1:])
                        T_len = int(len(temp_folder)+2)
    
                        for seed_folder in os.listdir(cwd):
    
                            print('    ', '    ', seed_folder)
    
                            if os.path.isdir(cwd+"/"+seed_folder):
                                cwd = current_dir+"/Minimisations/"+rings+"/"+temp_folder + "/" + seed_folder
    
                                seed_len = int(len(seed_folder) - 5)
                                seed = int(seed_folder[T_len:seed_len])
    
                                for protocol in os.listdir(cwd):
                                    cwd = current_dir+"/Minimisations/"+rings+"/"+temp_folder + "/" + seed_folder 
                                    print('    ', '    ', '    ', protocol)
    
                                    if os.path.isdir(cwd+"/"+protocol):
                                        cwd = current_dir+"/Minimisations/"+rings+"/" + \
                                            temp_folder + "/" + seed_folder + "/" + protocol
                                        len_protocol = len(protocol)+1
    
                                        if os.path.isdir(current_dir+"/Minimisations/"+rings+"/"+temp_folder + "/" + \
                                            seed_folder + "/" + protocol + "/" + \
                                            protocol + \
                                            '_{:}'.format(minimisation_number)):
                                            cwd = current_dir+"/Minimisations/"+rings+"/"+temp_folder + "/" + \
                                                seed_folder + "/" + protocol + "/" + \
                                                protocol + \
                                                '_{:}'.format(minimisation_number)
    
                                            for intercept in os.listdir(cwd):
    
                                                print('    ', '    ', '    ',
                                                      '    ', intercept)
    
                                                if os.path.isdir(current_dir+"/Minimisations/"+rings+"/"+temp_folder + "/" + seed_folder + "/" + protocol + "/" + protocol + '_{:}'.format(
                                                        minimisation_number) + "/" + intercept):
                                                    cwd = current_dir+"/Minimisations/"+rings+"/"+temp_folder + "/" + seed_folder + "/" + protocol + "/" + protocol + '_{:}'.format(
                                                        minimisation_number) + "/" + intercept
    
                                                    I = int(intercept[10:])
    
                                                    for areas in os.listdir(cwd):
    
                                                        print('    ', '    ', '    ',
                                                              '    ', '    ', areas)
    
                                                        if os.path.isdir(current_dir+"/Minimisations/"+rings+"/"+temp_folder + "/" + seed_folder + "/" + protocol + "/" + protocol + '_{:}'.format(
                                                                minimisation_number) + "/" + intercept + "/" + areas):
                                                            cwd = current_dir+"/Minimisations/"+rings+"/"+temp_folder + "/" + seed_folder + "/" + protocol + "/" + protocol + '_{:}'.format(
                                                                minimisation_number) + "/" + intercept + "/" + areas
                                                            area = int(areas[5:])
                                                            if previous_results == False:
        #                                                        print('No previous Results')
                                                                if 'engtot.out' in os.listdir(cwd):
                                                                    with open(cwd+'/engtot.out', 'r') as f:
                                                                        array = np.genfromtxt(
                                                                            f)
    
                                                                    E = array[-1, 2]
    
                                                                    output.write('{:<25}{:<25}{:<25}{:<25}{:<25}{:<25}'.format(
                                                                        int(rings), T, seed, protocol+'_{:}'.format(minimisation_number), I, area))
                                                                    output.write(
                                                                        '{:<10}\n'.format(E))
                                                                    print('++')
                                                                else:
                                                                    print('--')
    
                                                            else:
        #                                                        print('Found old run_data and using it')
                                                                to_write = [
                                                                    str(int(rings)), str(T), str(seed), str(protocol+'_{:}'.format(minimisation_number)), str(I), str(area)]
                                                                new_bool = True
        #                                                        print(run_data)
                                                                for i in range(run_data.shape[0]):
        #                                                            print(run_data[i,:], to_write)
                                                                    if to_write[0]==run_data[i,0] and to_write[1]==run_data[i,1] and to_write[2]==run_data[i,2] and to_write[2]==run_data[i,2] and to_write[3]==run_data[i,3] and to_write[4]==run_data[i,4] and to_write[5]==run_data[i,5]:
    
                                                                        new_bool = False
                                                                        print('==')
    
                                                                if new_bool == True:
                                                                        if 'engtot.out' in os.listdir(cwd):
                                                                            with open(cwd+'/engtot.out', 'r') as f:
                                                                                array = np.genfromtxt(
                                                                                    f)
                                                                            E = array[-1, 2]
                                                                            output.write('{:<25}{:<25}{:<25}{:<25}{:<25}{:<25}'.format(
                                                                                int(rings), T, seed, protocol+'_{:}'.format(minimisation_number), I, area))
                                                                            output.write(
                                                                                '{:<10}\n'.format(E))
                                                                            print('++')
                                                                        else:
                                                                            print('--')
