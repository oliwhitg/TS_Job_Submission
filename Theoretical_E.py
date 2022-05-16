import numpy as np
import os
import time

current_dir = os.getcwd()
print(current_dir)
current_dir = '/u/mw/jesu2929/Holiday_Work/TS_Test_Environment'
def pbc_r(v, d):
    for i in range(3):
        if v[i] > d[i] / 2:
            v[i] -= d[i]
        elif v[i] < -d[i] / 2:
            v[i] += d[i]

    return np.linalg.norm(v)

def harmlin(r, u0, r_cut):
    if r < r_cut:
        return (-u0/r_cut)*r + u0
    else:
        return 0

for minimisation_number in [10000]:

    if os.path.isfile(current_dir+"/{:}_theoretical_export.dat".format(minimisation_number)) == True:
        with open('{:}_theoretical_export.dat'.format(minimisation_number), 'r') as output:
    
            run_data = np.genfromtxt(output, dtype=str)
            print(run_data)
            time.sleep(1)
            previous_results = True
    else:
        previous_results = False
    
    
    
    
    with open('{:}_theoretical_export.dat'.format(minimisation_number), 'a+') as output:
    
        cwd = current_dir+"/Minimisations/"
    
        for rings in os.listdir(cwd):
            cwd = current_dir+"/Minimisations/"
 
            print(rings)
    
            if os.path.isdir(cwd + rings) and rings=='567':
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

                                                                if 'testout.rst' in os.listdir(cwd):
                                                                    print('Found Restart')
                                                                    with open(cwd + '/testout.rst', 'r') as f:
                                                                        crds = np.genfromtxt(f, skip_header=5, max_rows=4800)
                                                                    with open(cwd + '/harmpairs.dat', 'r') as f:
                                                                        harmpairs = np.genfromtxt(f,skip_header=1)
                                                                    with open(cwd + '/crys.crds', 'r') as f:
                                                                        d = np.genfromtxt(f,
                                                                                             skip_header=(5 + 4800 + 3))
                                                                    with open(cwd+'/harm_au.inpt', 'r') as f:
                                                                        for i in range(8):
                                                                            skip = f.readline()
                                                                        u0 = float(f.readline().strip().split()[0][:-2])
                                                                        r_cut = float(f.readline().strip().split()[0][:-2])

                                                                    E_harm = 0
                                                                    r_mx = 3.0387
                                                                    r_xx = 4.9624
                                                                    for i in range(harmpairs.shape[0]):
                                                                        atom1, atom2 = int(harmpairs[i,0]), int(harmpairs[i,1])
                                                                        if atom1 > 3200 or atom2 > 3200:
                                                                            r0 = r_mx

                                                                        else:
                                                                            r0 = r_xx
                                                                        v = np.subtract(crds[atom1-1,:], crds[atom2-1,:])
                                                                        r = pbc_r(v,d)
                                                                        E_harm += 0.5*1.001*(r-r0)**2
                                                                    ## xx
                                                                    E_xx = 0
                                                                    for atom1 in range(3199):
                                                                        for atom2 in range(atom1+1, 3200):
                                                                            v = np.subtract(crds[atom1, :],crds[atom2, :])
                                                                            r = pbc_r(v,d)
                                                                            E_xx += harmlin(r,u0, r_cut)

                                                                    ## mm
                                                                    E_mm = 0
                                                                    for atom1 in range(3200, 4799):
                                                                        for atom2 in range(atom1+1, 4800):
                                                                            v = np.subtract(crds[atom1, :],crds[atom2, :])
                                                                            r = pbc_r(v,d)
                                                                            E_mm += harmlin(r,u0, r_cut)


                                                                    output.write('{:<25}{:<25}{:<25}{:<25}{:<25}{:<25}'.format(
                                                                        int(rings), T, seed, protocol+'_{:}'.format(minimisation_number), I, area))
                                                                    output.write(
                                                                        '{:<10}{:<10}{:<10}\n'.format(E_harm,E_xx,E_mm))
                                                                    print('++')
                                                                else:
                                                                    print('--')
    
#                                                            else:
#        #                                                        print('Found old run_data and using it')
#                                                                to_write = [
#                                                                    str(int(rings)), str(T), str(seed), str(protocol+'_{:}'.format(minimisation_number)), str(I), str(area)]
#                                                                new_bool = True
#        #                                                        print(run_data)
#                                                                for i in range(run_data.shape[0]):
#        #                                                            print(run_data[i,:], to_write)
#                                                                    if to_write[0]==run_data[i,0] and to_write[1]==run_data[i,1] and to_write[2]==run_data[i,2] and to_write[2]==run_data[i,2] and to_write[3]==run_data[i,3] and to_write[4]==run_data[i,4] and to_write[5]==run_data[i,5]:
#
#                                                                        new_bool = False
#                                                                        print('==')
#
#                                                                if new_bool == True:
#                                                                        if 'engtot.out' in os.listdir(cwd):
#                                                                            with open(cwd+'/engtot.out', 'r') as f:
#                                                                                array = np.genfromtxt(
#                                                                                    f)
#                                                                            E = array[-1, 2]
#                                                                            output.write('{:<25}{:<25}{:<25}{:<25}{:<25}{:<25}'.format(
#                                                                                int(rings), T, seed, protocol+'_{:}'.format(minimisation_number), I, area))
#                                                                            output.write(
#                                                                                '{:<10}\n'.format(E))
#                                                                            print('++')
#                                                                        else:
#                                                                            print('--')
