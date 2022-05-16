import numpy as np
import os
import shutil
import glob
import subprocess
import time
import sys
import getopt
import re
from shutil import copytree
from shutil import copy

def import_aux():
    ##################################################################
    #import aux
    with open('run.aux', 'r') as f:
        skip            = f.readline()
        output          = str(f.readline().split()[0])
        skip            = f.readline()
        dry_run         = (f.readline().split()[0])
        skip            = f.readline()
        run_type        = str(f.readline().split()[0])
        run_string      = f.readline().split()[0]
        rescale         = f.readline().split()[0]
        run_steps       = int(f.readline().split()[0])
        skip            = f.readline()
        restart_from    = f.readline().split()[0]
        input_string    = f.readline().split()[0]
        skip            = f.readline()
        low_ring_size   = int(f.readline().split()[0])
        upper_ring_size = int(f.readline().split()[0])
        skip            = f.readline()
        MC_T            = int(1000*float(f.readline().split()[0]))
        seed_start      = int(f.readline().split()[0])
        seed_end        = int(f.readline().split()[0])
        skip            = f.readline()
        alpha_target    = float(f.readline().split()[0])
        alpha_lenience  = float(f.readline().split()[0])
        p6_target       = float(f.readline().split()[0])
        p6_lenience     = float(f.readline().split()[0])
        skip            = f.readline()
        upper_A         = float(f.readline().split()[0])
        lower_A         = float(f.readline().split()[0])
        step_A          = float(f.readline().split()[0])
        skip            = f.readline()
        repulsion       = str(f.readline().split()[0])
        repulsion_except= int(f.readline().split()[0])
        repulsion_min   = float(f.readline().split()[0])
        repulsion_max   = float(f.readline().split()[0])
        repulsion_step  = float(f.readline().split()[0])

        skip            = f.readline()
        repair          = str(f.readline().split()[0])

        catchup         = str(f.readline().split()[0])
        catchup_step    = str(f.readline().split()[0])

    ###################################################################
    #convert global variables to local variables
    options = {}

    print(dry_run)
    if dry_run in ['True', 't', '1']:
        dry_run = True
    else:
        dry_run = False

    if catchup in ['True', 't', '1']:
        catchup = True
    else:
        catchup = False

    if run_type=='initial':
        options["run_type"] =0
    elif run_type=='repeat':
        options["run_type"]=1
    elif run_type=='novo':
        options["run_type"]=2
    else:
        print('Unrecognised run type')
        return 1


    # Ouput Prefix
    if not os.path.isdir(output):
        os.mkdir(output)
    rundir = output

    if restart_from=='no':
        options['restart_from'] = 0
    elif restart_from=='harmlin' or restart_from=='harmonic':
        options['restart_from'] = 1
    elif restart_from=='static' or restart_from=='s':
        options['restart_from'] = 2
    elif restart_from=='dynamic':
        options['restart_from'] = 3

    #Ring Size String eg. 567
    ring_size_string = ""
    i = low_ring_size
    while i <= upper_ring_size:
        ring_size_string=ring_size_string+"{:}".format(i)
        i+=1

    # Update rundir to include Ring size string
    if not os.path.isdir(rundir+'/'+ring_size_string+'/'):
        os.mkdir(rundir+'/'+ring_size_string+'/')
    rundir = rundir+'/'+ring_size_string

    options['output_prefix'] = rundir

    if run_string=='harmonic':
        options["run"] = 0
        print('Harmonic potential')

    elif run_string=='static':
        options["run"] = 1
        print('Static Ions')
        print('Beware : running this without harmonic potential will slow down calculations')

    elif run_string=='dynamic':
        options["run"] = 2
        print('Beware : running this without harmonic potential will slow down calculations')
    else:
        print('unrecognised run type')

    #############################################################################################

    if input_string=='netmc':
        if os.path.isdir(options['output_prefix'] + '/T{:}/'.format(MC_T))==False:
            os.mkdir(options['output_prefix'] + '/T{:}/'.format(MC_T))
        print('MC T : ', MC_T)

        for i in range(seed_start, seed_end+1):
            print('         Seed : ', i)
            in_dir = os.getcwd()+'/Results/'+ring_size_string+'/T{:}/t{:}_s{:}_same'.format(MC_T, MC_T, i)
            in_dir_novo = os.getcwd()+'/Minimisations/'+ring_size_string+'/T{:}/t{:}_s{:}_same'.format(MC_T, MC_T, i)
            out_dir = options['output_prefix'] + '/T{:}/t{:}_s{:}_same'.format(MC_T, MC_T, i)
            if os.path.isdir(out_dir)==False:
                os.mkdir(out_dir)
            if catchup == True:
                catchup_run(out_dir)
            else:
                ##########################################################################################################################################################
                if repair==True:
                    repair_run(out_dir, run_steps, False)
                ##########################################################################################################################################################
                print('Run type : ', run_type)
                if options["run_type"] == 0:
                    print('Initial run will start from netmc results, and run a harmonic energy minimisation')
                    initial_run(in_dir, options['restart_from'], out_dir, run_string, repulsion_min, repulsion_max, repulsion_step, lower_A, upper_A, step_A, run_steps, dry_run)
                ##########################################################################################################################################################
                elif options["run_type"] ==1:
                    print('Repeat run will take the results from a previous run, and run additional steps restarting from the last run')
                    print('This result will start from results \$POTENTIAL_n and run to \$POTENTIAL_n+1')
                    print('We will be running {:} potential'.format(run_string))
                    repeat_run(in_dir, options['restart_from'], out_dir, run_string, repulsion_min, repulsion_max, repulsion_step, lower_A, upper_A, step_A, run_steps, dry_run)
                ##########################################################################################################################################################
                elif options["run_type"] ==2:
                    print('Novo run takes a previously minimised structure and feeds it into a new potential')
                    print('We will be running a {:} potential on a structure already minimised by {:}'.format(run_string, restart_from))
                    novo_run(in_dir_novo, options['restart_from'], out_dir, run_string, repulsion_min, repulsion_max, repulsion_step, lower_A, upper_A, step_A, run_steps, dry_run)
                ##########################################################################################################################################################



def initial_run(in_dir, restart_from, out_dir, run_string, repulsion_min, repulsion_max, repulsion_step, lower_A, upper_A, step_A, run_steps, dry_run):
    if run_string == 'harmonic':

        if os.path.isdir(out_dir + '/' + run_string) == False:
            os.mkdir(out_dir + '/' + run_string)

            print('Entirely new run')

        if os.path.isdir(out_dir + '/' + run_string) == True:
            print('Already got a harmonic, will fill in the gaps')
            #shutil.rmtree(out_dir + '/' + run_string)
        os.mkdir(out_dir + '/' + run_string + '/' + run_string + '_' + str(run_steps))
        current_running_dir = out_dir + '/' + run_string + '/' + run_string + '_' + str(run_steps)
        ###################################################

        for intercept_100 in range(int(repulsion_min * 100), int(100 * (repulsion_max + 2 * repulsion_step) + 1),
                                   int(100 * repulsion_step)):

            if int(100 * (repulsion_max + repulsion_step)) <intercept_100 <= int(100 * (repulsion_max + 2* repulsion_step)):
                intercept = 0
            else:
                intercept = intercept_100 / 100.0

            print('                                 Intercept : ', intercept)

#            if os.path.isdir(out_dir + '/' + run_string + '/' + run_string + '_0' + '/intercept_{:}'.format(
#                    int(intercept*100))) == False:

            os.mkdir(out_dir + '/' + run_string + '/' + run_string + '_' + str(run_steps) + '/intercept_{:}'.format(int(intercept*100)))
            current_running_dir = out_dir + '/' + run_string + '/' + run_string + '_' + str(run_steps) + '/intercept_{:}'.format(int(intercept*100))
            for area_1000 in range(int(1000 * lower_A), int(1000 * upper_A + 1), int(1000 * step_A)):
                area = area_1000 / 1000.0


                print('                                             Area : ', area_1000)

                if os.path.isdir(out_dir + '/' + run_string + '/' + run_string + '_' + str(run_steps) + '/intercept_{:}'.format(
                        int(intercept*100)) + '/area_{:}/'.format(int(area_1000))) == False:

                    os.mkdir(out_dir + '/' + run_string + '/' + run_string + '_' + str(run_steps) + '/intercept_{:}'.format(
                        int(intercept*100)) + '/area_{:}/'.format(int(area_1000)))
                    current_running_dir = out_dir + '/' + run_string + '/' + run_string + '_' + str(run_steps) + '/intercept_{:}'.format(
                        int(intercept*100)) + '/area_{:}/'.format(int(area_1000))

                    ## unconvinced by indir
                    start_harmlin(True, in_dir, current_running_dir, intercept, area, 0, run_steps,
                                  dry_run)

                else:
                    print('Already done area {:}'.format(area_1000))

#            else:
#                print('Already done intercept {:}'.format(intercept*100))


    ###############################################
    else:
        print('Trying to run high level potential without harmonics first')
        print('Expected to blow up, please run harmonics first')
        print('set harm steps to 1 if you need to do this')
        return 1

    return

def catchup_run(out_dir):
    dry_run=False
    for run_string in os.listdir(out_dir):
        print('Catching up {:}'.format(run_string))
        run_string_length = len(run_string) + 1
        dict={}
        for runs in os.listdir(out_dir + '/' + run_string):
            if '{:}_'.format(run_string) in runs:
                restart_run_count = int(runs[run_string_length:])
                print(restart_run_count)
                dict['{:}'.format(restart_run_count)] = {}

        print('We have {:} total runs'.format(len(dict.keys())))
        print('[',[float(i) for i in dict.keys()],']')
        for runs in dict.keys():
            print('Run {:}:'.format(runs))
            if 1>0:
                restart_run_count = int(runs)

                for intercepts in os.listdir(
                        out_dir + '/' + run_string + '/' + run_string + '_{:}'.format(restart_run_count)):

                    if 'intercept_' in intercepts:
                        intercept = int(intercepts[10:])
                        print('    Intercept : {:}'.format(intercept))
                        dict['{:}'.format(runs)]['{:}'.format(int(intercept))] = []
                        for area_1000 in os.listdir(out_dir + '/' + run_string + '/' + run_string + '_{:}'.format(
                                restart_run_count) + '/intercept_{:}'.format(int(intercept))):
                            if 'area_' in area_1000:
                                area = int(area_1000[5:])

                                if 'testout.rst' in os.listdir(
                                        out_dir + '/' + run_string + '/' + run_string + '_{:}'.format(
                                                restart_run_count) + '/intercept_{:}'.format(
                                                int(intercept)) + '/area_{:}'.format(int(area))):
                                    dict['{:}'.format(runs)]['{:}'.format(int(intercept))].append(area)
                        print('        Areas : ', dict['{:}'.format(runs)]['{:}'.format(int(intercept))])
        ##################################
        ## correct intercepts
        ## Sorting through runs in reverse order of size

        sorted_runs = [int(i) for i in dict.keys()]
        sorted_runs.sort(reverse=True)
        ref_run = max(sorted_runs)

        modded_dict = dict.copy()

        for runs_counter in range(1,len(sorted_runs)):

            runs = sorted_runs[runs_counter]
            steps_to_run = ref_run - runs

            run_intercepts = modded_dict['{:}'.format(runs)].keys()
            ref_intercepts = modded_dict['{:}'.format(ref_run)].keys()

            ############################################################################################################
            ## this only works for catching up intercepts, areas to come
            to_run_dict = {}
            to_run_intercepts = list(set(run_intercepts) - set(ref_intercepts))

            for intercepts in run_intercepts:
                if intercepts in to_run_intercepts:
                    to_run_dict['{:}'.format(intercepts)] = dict['{:}'.format(ref_run)]['{:}'.format(intercepts)]

                    ########################################################################################################

                    ## modded version of repeat run

                    os.mkdir(out_dir + '/' + run_string + '/' + run_string + '_{:}'.format(ref_run) + '/' + intercepts)
                    #current_running_dir = out_dir + '/' + run_string + '/' + run_string + '_{:}'.format(ref_run) + '/' + intercepts)

                    modded_dict['{:}'.format(ref_run)]['{:}'.format(intercepts)] = []

                    for area in dict['{:}'.format(ref_run)]['{:}'.format(intercepts)]:

                        os.mkdir(out_dir + '/' + run_string + '/' + run_string + '_{:}'.format(ref_run) + '/' + intercepts + '/' + area)
                        current_running_dir = out_dir + '/' + run_string + '/' + run_string + '_{:}'.format(ref_run) + '/' + intercepts + '/' + area

                        for files in os.listdir(out_dir + '/' + run_string + '/' + run_string + '_{:}'.format(runs) + '/' + intercepts + '/' + area):
                            copy( out_dir + '/' + run_string + '/' + run_string + '_{:}'.format(runs) + '/' + intercepts + '/' + area + '/' + files, current_running_dir)

                        intercept_100 = float(intercepts[10:])
                        area_1000 = float(area[5:])

                        if run_string == 'harmonic':

                            start_harmlin(False, 'blank', current_running_dir, float(int(intercept_100) / 100), float(int(area_1000) / 1000), ref_run, steps_to_run, dry_run)
                        elif run_string == 'static':
                            start_TS_static(False, 'blank', current_running_dir, float(int(intercept_100) / 100), float(int(area_1000) / 1000), ref_run, steps_to_run, dry_run)


                        modded_dict['{:}'.format(ref_run)]['{:}'.format(intercepts)].append(area)
                else:
                    run_intercepts_areas = modded_dict['{:}'.format(runs)]['{:}'.format(intercepts)]
                    ref_intercepts_areas = modded_dict['{:}'.format(ref_run)]['{:}'.format(intercepts)]

                    to_run_intercepts_areas = list(set(run_intercepts_areas) - set(ref_intercepts_areas))

                    for area in to_run_intercepts_areas:
                        os.mkdir(out_dir + '/' + run_string + '/' + run_string + '_{:}'.format(int(ref_run)) + '/intercept_{:}'.format(intercepts) + '/area_{:}'.format(area))
                        current_running_dir = out_dir + '/' + run_string + '/' + run_string + '_{:}'.format(int(ref_run)) + '/intercept_{:}'.format(intercepts) + '/area_{:}/'.format(area)

                        for files in os.listdir(out_dir + '/' + run_string + '/' + run_string + '_{:}'.format(int(runs)) + '/intercept_{:}'.format(intercepts) + '/area_{:}/'.format(area)):
                            #print(out_dir + '/' + run_string + '/' + run_string + '_{:}'.format(int(runs)) + '/intercept_{:}'.format(intercepts) + '/area_{:}/'.format(area)+files)
                            copy(out_dir + '/' + run_string + '/' + run_string + '_{:}'.format(int(runs)) + '/intercept_{:}'.format(intercepts) + '/area_{:}/'.format(area)+files, current_running_dir)

                        intercept_100 = intercepts
                        area_1000 = area

                        if run_string == 'harmonic':

                            start_harmlin(False, 'blank', current_running_dir, float(int(intercept_100) / 100),float(int(area_1000) / 1000), ref_run, steps_to_run, dry_run)
                        elif run_string == 'static':
                            start_TS_static(False, 'blank', current_running_dir, float(int(intercept_100) / 100),float(int(area_1000) / 1000), ref_run, steps_to_run, dry_run)

                        modded_dict['{:}'.format(ref_run)]['{:}'.format(intercepts)].append(area)


        ################################################################################################################


            ## add runs done to list of intercepts, areas in ref_run



    return 0

def repair_run(out_dir, run_steps, dry_run):
    if os.path.isdir(out_dir) == False:
        print('Nothing run yet, nothing to repair')
        return 0
    print('Repair only corrects the most recent run')
    run_bool = False
    for run_string in os.listdir(out_dir):

        if run_bool == False:           #This prevents simultaneously running different repair jobs

            print('Repairing {:} ...'.format(run_string))
            run_string_length = len(run_string) + 1
            list_runs = []
            for folders in os.listdir(out_dir + '/' + run_string):
                if '{:}_'.format(run_string) in folders:
                    list_runs.append(int(folders[run_string_length:]))
            restart_run_count = int(max(list_runs))



            for intercepts in os.listdir(out_dir + '/' + run_string + '/' + run_string + '_{:}'.format(restart_run_count)):
                if 'intercept_' in intercepts:
                    intercept_100 = int(intercepts[10:])
                    for areas in os.listdir(out_dir + '/' + run_string + '/' + run_string + '_{:}'.format(restart_run_count)+'/intercept_{:}'.format(int(intercept_100))):
                        if 'area_' in areas:
                            area = int(areas[5:])
                            if 'testout.rst' not in os.listdir(out_dir + '/' + run_string + '/' + run_string + '_{:}'.format(restart_run_count)+'/intercept_{:}'.format(int(intercept_100))+'/area_{:}'.format(int(area))):
                                print('Incomplete run')
                                if 'crys.crds' not in os.listdir(out_dir + '/' + run_string + '/' + run_string + '_{:}'.format(restart_run_count)+'/intercept_{:}'.format(int(intercept_100))+'/area_{:}'.format(int(area))):
                                    print('No starting run, something is fucked, imma delete this area')
                                    shutil.rmtree(out_dir + '/' + run_string + '/' + run_string + '_{:}'.format(restart_run_count)+'/intercept_{:}'.format(int(intercept_100))+'/area_{:}'.format(int(area)))
                                    break
                                else:
                                    current_running_dir = out_dir + '/' + run_string + '/' + run_string + '_{:}'.format(restart_run_count)+'/intercept_{:}'.format(int(intercept_100))+'/area_{:}/'.format(int(area))
                                    print('Calculations not yet run')
                                    run_bool = True
                                    if run_string == 'harmonic':
                                        start_harmlin(False, 'blank', current_running_dir, float(intercept_100 / 100),
                                                      float(area / 1000), restart_run_count+run_steps,
                                                      run_steps, dry_run)
                                    elif run_string == 'static':
                                        start_TS_static(False, 'blank', current_running_dir, float(intercept_100 / 100),
                                                        float(area / 1000), restart_run_count+run_steps,
                                                        run_steps, dry_run)
    return 0



def repeat_run(in_dir, restart_from, out_dir, run_string, repulsion_min, repulsion_max, repulsion_step, lower_A, upper_A, step_A, run_steps, dry_run):
    if os.path.isdir(out_dir + '/' + run_string) == False:
        print('Don\'t know where to restart this from as requested run ({:}) does not exist'.format(out_dir+'/'+run_string))
        return 1
    else:

        run_string_length = len(run_string) + 1
        list_runs = []
        for folders in os.listdir(out_dir + '/' + run_string):
            if '{:}_'.format(run_string) in folders:
                list_runs.append(int(folders[run_string_length:]))
        restart_run_count = int(max(list_runs))
        #            copytree(out_dir+'/'+run_string+'/'+run_string+'_{:}'.format(restart_run_count), out_dir+'/'+run_string+'/'+run_string+'_{:}'.format(restart_run_count+1))


        ##############################################################################

        #should probably use the intercept and seed from previous runs, not make new ones
        dict_allowed_runs = {}
        for intercepts in os.listdir(out_dir + '/' + run_string + '/' + run_string + '_{:}'.format(restart_run_count)):
            if 'intercept_' in intercepts:
                intercept = int(intercepts[10:])
                dict_allowed_runs['{:}'.format(int(intercept))] = []
                for area_1000 in os.listdir(out_dir + '/' + run_string + '/' + run_string + '_{:}'.format(restart_run_count)+'/intercept_{:}'.format(int(intercept))):
                    if 'area_' in area_1000:
                         area = int(area_1000[5:])
                         if 'testout.rst' in os.listdir(out_dir + '/' + run_string + '/' + run_string + '_{:}'.format(restart_run_count)+'/intercept_{:}'.format(int(intercept))+'/area_{:}'.format(int(area))):
                             dict_allowed_runs['{:}'.format(int(intercept))].append(area)

        ##############################################################################

        os.mkdir(out_dir + '/' + run_string + '/' + run_string + '_{:}'.format(restart_run_count+run_steps))
        current_running_dir = out_dir + '/' + run_string + '/' + run_string + '_{:}'.format(restart_run_count+run_steps)
        
        ##############################################################################

        for intercept_100 in dict_allowed_runs.keys():
            
            os.mkdir(out_dir + '/' + run_string + '/' + run_string + '_{:}'.format(
                restart_run_count+run_steps) + '/intercept_{:}'.format(int(intercept_100)))
            current_running_dir = out_dir + '/' + run_string + '/' + run_string + '_{:}'.format(
                restart_run_count+run_steps) + '/intercept_{:}'.format(int(intercept_100))

            for area_1000 in dict_allowed_runs['{:}'.format(int(intercept_100))]:
                os.mkdir(out_dir + '/' + run_string + '/' + run_string + '_{:}'.format(
                    restart_run_count+run_steps) + '/intercept_{:}'.format(int(intercept_100)) + '/area_{:}'.format(
                    int(area_1000)))
                current_running_dir = out_dir + '/' + run_string + '/' + run_string + '_{:}'.format(
                    restart_run_count+run_steps) + '/intercept_{:}'.format(int(intercept_100)) + '/area_{:}/'.format(
                    int(area_1000))


                for files in os.listdir(
                        out_dir + '/' + run_string + '/' + run_string + '_{:}/intercept_{:}/area_{:}/'.format(
                                restart_run_count , intercept_100, area_1000)):
                    copy(out_dir + '/' + run_string + '/' + run_string + '_{:}/intercept_{:}/area_{:}/'.format(
                        restart_run_count, intercept_100, area_1000) + files, current_running_dir)
                if run_string == 'harmonic':

                    start_harmlin(False, 'blank', current_running_dir, float(int(intercept_100)/100), float(int(area_1000)/1000), restart_run_count+run_steps,
                                  run_steps, dry_run)
                elif run_string == 'static':
                    start_TS_static(False, 'blank', current_running_dir, float(int(intercept_100)/100), float(int(area_1000)/1000), restart_run_count+run_steps,
                                    run_steps, dry_run)



                

#        ##############################################################################        
#        for intercept_100 in range(int(repulsion_min * 100), int(100 * (repulsion_max + repulsion_step) + 1),
#                                   int(100 * repulsion_step)):
#            if intercept_100 == int(100 * (repulsion_max + repulsion_step)):
#                intercept = 0
#                intercept_100 = 0
#            else:
#                intercept = intercept_100 / 100
#
#            os.mkdir(out_dir + '/' + run_string + '/' + run_string + '_{:}'.format(
#                restart_run_count) + '/intercept_{:}'.format(int(intercept_100)))
#            current_running_dir = out_dir + '/' + run_string + '/' + run_string + '_{:}'.format(
#                restart_run_count) + '/intercept_{:}'.format(int(intercept_100))
#            for area_1000 in range(int(1000 * lower_A), int(1000 * upper_A + 1), int(1000 * step_A)):
#                area = area_1000 / 1000
#                os.mkdir(out_dir + '/' + run_string + '/' + run_string + '_{:}'.format(
#                    restart_run_count) + '/intercept_{:}'.format(int(intercept_100)) + '/area_{:}'.format(
#                    int(area_1000)))
#                current_running_dir = out_dir + '/' + run_string + '/' + run_string + '_{:}'.format(
#                    restart_run_count) + '/intercept_{:}'.format(int(intercept_100)) + '/area_{:}/'.format(
#                    int(area_1000))
#                if os.path.isdir(out_dir + '/' + run_string + '/' + run_string + '_{:}/intercept_{:}/area_{:}/'.format(
#                        restart_run_count - 1, intercept_100, area_1000)) == True:
#                    for files in os.listdir(
#                            out_dir + '/' + run_string + '/' + run_string + '_{:}/intercept_{:}/area_{:}/'.format(
#                                    restart_run_count - 1, intercept_100, area_1000)):
#                        copy(out_dir + '/' + run_string + '/' + run_string + '_{:}/intercept_{:}/area_{:}/'.format(
#                            restart_run_count - 1, intercept_100, area_1000) + files, current_running_dir)
#                    if run_string == 'harmonic':
#                        start_harmlin(False, 'blank', current_running_dir, intercept, area, restart_run_count,
#                                      run_steps, dry_run)
#                    elif run_string == 'static':
#                        start_TS_static(False, 'blank', current_running_dir, intercept, area, restart_run_count,
#                                        run_steps, dry_run)
#                else:
#                    print('Failed to restart as from dir ({:}) doesn\'t exist'.format(out_dir + '/' + run_string + '/' + run_string + '_{:}/intercept_{:}/area_{:}/'.format(
#                        restart_run_count - 1, intercept_100, area_1000)))


    return 0

def novo_run(in_dir, restart_from, out_dir, run_string, repulsion_min, repulsion_max, repulsion_step, lower_A, upper_A, step_A, run_steps, dry_run):
    ### yeah this is clearly fucked
    print('Novo run')
    if restart_from==0:
        print('Error in calling novo here, must restart from a valid configuration')
    elif restart_from==1:
        restart_from_string = "harmonic"
    elif restart_from==2:
        restart_from_string = "static"
    elif restart_from==3:
        restart_from_string = "dynamic"
    else:
        print('FUCK')    

    ##################################################

    run_string_length = len(restart_from_string) + 1
    list_runs = []

    for folders in os.listdir(in_dir + '/' + restart_from_string):
        if '{:}_'.format(restart_from_string) in folders:
            list_runs.append(int(folders[run_string_length:]))
    print(list_runs)
    restart_run_count = int(max(list_runs))

    ##############################################################################

    # should probably use the intercept and seed from previous runs, not make new ones
    dict_allowed_runs = {}
    for intercepts in os.listdir(out_dir + '/' + restart_from_string + '/' + restart_from_string + '_{:}'.format(restart_run_count)):
        if 'intercept_' in intercepts:
            intercept = int(intercepts[10:])
            dict_allowed_runs['{:}'.format(int(intercept))] = []
            for area_1000 in os.listdir(out_dir + '/' + restart_from_string + '/' + restart_from_string + '_{:}'.format(
                    restart_run_count) + '/intercept_{:}'.format(int(intercept))):
                if 'area_' in area_1000:
                    area = int(area_1000[5:])
                    if 'testout.rst' in os.listdir(out_dir + '/' + restart_from_string + '/' + restart_from_string + '_{:}'.format(
                    restart_run_count) + '/intercept_{:}'.format(int(intercept)) + '/area_{:}'.format(
                            int(area))):
                        dict_allowed_runs['{:}'.format(int(intercept))].append(area)

    ##############################################################################

    if os.path.isdir(out_dir + '/' + run_string) == True:
        print('Error : run string already exists')
        print('Delete previous run and rerun?')
        break_input = input()
        if break_input == 'y' or break_input == 'yes':
            print('Ok, deleting ...')
            shutil.rmtree(out_dir + '/' + run_string)
        else:
            print('Breaking')
            return 1

    os.mkdir(out_dir + '/' + run_string)
    os.mkdir(out_dir + '/' + run_string + '/' + run_string + '_' + str(run_steps))
    current_running_dir = out_dir + '/' + run_string + '/' + run_string + '_' + str(run_steps)

    ###################

    for intercept_100 in dict_allowed_runs.keys():

        os.mkdir(out_dir + '/' + run_string + '/' + run_string + '_' +
            str(run_steps) + '/intercept_{:}'.format(int(intercept_100)))
        current_running_dir = out_dir + '/' + run_string + '/' + run_string + '_' + \
                              str(run_steps) + '/intercept_{:}'.format(int(intercept_100))

        for area_1000 in dict_allowed_runs['{:}'.format(int(intercept_100))]:
            os.mkdir(out_dir + '/' + run_string + '/' + run_string + '_' + str(run_steps) + '/intercept_{:}'.format(int(intercept_100)) + '/area_{:}'.format(
                int(area_1000)))
            current_running_dir = out_dir + '/' + run_string + '/' + run_string + '_' + str(run_steps) + '/intercept_{:}'.format(int(intercept_100)) + '/area_{:}/'.format(
                int(area_1000))

            for files in os.listdir(
                    out_dir + '/' + restart_from_string + '/' + restart_from_string + '_{:}/intercept_{:}/area_{:}/'.format(
                        restart_run_count, intercept_100, area_1000)):
                copy(out_dir + '/' + restart_from_string + '/' + restart_from_string + '_{:}/intercept_{:}/area_{:}/'.format(
                    restart_run_count, intercept_100, area_1000) + files, current_running_dir)

            if run_string == 'static':
                start_TS_static(True, in_dir, current_running_dir, float(int(intercept_100)/100), area, restart_run_count,
                                run_steps, dry_run)
            else:
                start_TS_dynamic(True, in_dir, current_running_dir, intercept, float(int(area_1000)/1000), 0, run_steps,
                                 dry_run)


    return 0






def start_harmlin(new, input_folder, running_output_dir, intercept, area, restart_run_count, run_steps, dry_run):

    if new==False:
        for filename in glob.glob(running_output_dir+'fort*'):
            os.remove(filename)
        for filename in glob.glob(running_output_dir+'*.out*'):
            os.remove(filename)
        for filename in glob.glob(running_output_dir + 'header*'):
            os.remove(filename)

        if os.path.isfile(running_output_dir + 'testout.rst')==True:
            shutil.move(running_output_dir + 'testout.rst', running_output_dir + 'crys.crds{:}'.format(restart_run_count))
            copy(running_output_dir + 'filenames.inpt', running_output_dir + 'filenames_original.inpt')
            with open(running_output_dir + 'filenames.inpt', 'w') as f:
                with open(running_output_dir + 'filenames_original.inpt', 'r') as g:
                    for i in range(5):
                        f.write(g.readline())
                    f.write('crys.crds{:}\n'.format(restart_run_count))
                    skip = g.readline()
                    for i in range(42 - 6):
                        f.write(g.readline())
            shutil.move(running_output_dir + 'harm_au.inpt', running_output_dir + 'old_harm_au.inpt')
            with open(running_output_dir + 'old_harm_au.inpt', 'r') as f:
                with open(running_output_dir + 'harm_au.inpt', 'w') as g:
                    for i in range(8):
                        g.write(f.readline())
                    old_repulsion = f.readline().strip().split()
                    val = float(old_repulsion[0][:-2])/2
                    print(old_repulsion, '    ', val)
                    g.write('{:}d0     u0 (11)\n'.format(val))
                    for i in range(3):
                        g.write(f.readline())
    
                    g.write('{:}d0     u0 (22)\n'.format(val))
                    skip = f.readline()
                    g.write(f.readline())
        else:
            print('Previous run at {:} failed!'.format(running_output_dir))

    # copy important harmlin files
    # code must be compiled on dirac first
    else:
        restart_run_count=run_steps
        for files in os.listdir(input_folder):
            copy(input_folder + '/' + files, running_output_dir)
        copy(os.getcwd() + '/harmlin/inputs/cf.inpt', running_output_dir)
        copy(os.getcwd() + '/harmlin/inputs/lightcf_licl.inpt', running_output_dir)
        copy(os.getcwd() + '/harmlin/inputs/filenames.inpt', running_output_dir)
        copy(os.getcwd() + '/harmlin/inputs/crystal_cell.inpt', running_output_dir)
        #
        copy(os.getcwd() + '/harmlin/inputs/harm_au.inpt', running_output_dir)
        if intercept==0:
            modify_harmlin(running_output_dir, False, float(intercept))
        else:
            modify_harmlin(running_output_dir, True, float(intercept))

        #
        copy(os.getcwd() + '/harmlin/inputs/runtime_licl.inpt', running_output_dir)
        #modify_runtime_harmlin(running_output_dir, run_steps, restart)
        copy(os.getcwd() + '/harmlin/inputs/combined_tersoff_3d_lj_D_zmove6000_harm_nodrsav_lin_2021.x', running_output_dir)

        scale_old_crds_harmpairs(running_output_dir, area)
    modify_runtime_harmlin(running_output_dir, run_steps)
    with open(running_output_dir + '/sub_script', 'w') as f:
        f.write('#!/bin/csh\n#$ -j y\n#$ -o .\n#$ -l s_rt=700:00:00\n#$ -l h_rt=700:00:00\n')
        f.write('#$ -N harmonic{:}\n\n'.format(restart_run_count))
        f.write('set fromdir = jesu2929@coulson.chem.ox.ac.uk:{:}\n'.format(running_output_dir))
        f.write('set todir   = jesu2929@coulson.chem.ox.ac.uk:{:}\n'.format(running_output_dir))

        f.write('setenv WORK /tmp/$USER/$JOB_ID\n\n')

        f.write('mkdir -p $WORK\n')
        f.write('scp -r $fromdir\* $WORK\n')
        
        f.write('cd $WORK\n')
        f.write('./combined_tersoff_3d_lj_D_zmove6000_harm_nodrsav_lin_2021.x\n')
        f.write('rm disp.out\n')
        f.write('scp *.dat *.log $fromdir\n')
        f.write('scp * $todir\n')
        f.write('rm -Rf $WORK\n')
    if dry_run == False:
        sub_string = str(running_output_dir + 'sub_script')

        subprocess.call(['scp', sub_string, 'jesu2929@dirac.chem.ox.ac.uk:~/From_Coulson/.'])
        subprocess.call(['ssh', 'jesu2929@dirac.chem.ox.ac.uk', 'qsub', 'From_Coulson/sub_script']) 


def start_TS_static(new, input_folder, running_output_dir, intercept, area, restart_run_count, run_steps, dry_run):
    for filename in glob.glob(running_output_dir + 'fort*'):
        os.remove(filename)
    for filename in glob.glob(running_output_dir + '*.out*'):
        os.remove(filename)
    for filename in glob.glob(running_output_dir + 'header*'):
        os.remove(filename)


    if new==True:
        #restart_run_count = 0
        ## remove harmlin ones
        for files in os.listdir(running_output_dir):
            if files in os.listdir(os.getcwd()+'/harmlin/inputs'):
                os.remove(running_output_dir + '/' + files)
        for files in os.listdir(os.getcwd()+'/TS_Static/test2'):
            copy(os.getcwd()+'/TS_Static/test2/'+files, running_output_dir)

    print('Restart Run Count : ', restart_run_count)
    shutil.move(running_output_dir + 'testout.rst', running_output_dir + 'crys.crds{:}'.format(restart_run_count))
    copy(running_output_dir + 'filenames.inpt', running_output_dir + 'filenames_original.inpt')

    #if new==False:
    with open(running_output_dir + 'filenames.inpt', 'w') as f:
        with open(running_output_dir + 'filenames_original.inpt', 'r') as g:
            for i in range(5):
                f.write(g.readline())
            f.write('crys.crds{:}\n'.format(restart_run_count))
            skip = g.readline()
            for i in range(42 - 6):
                f.write(g.readline())
    modify_runtime_static(running_output_dir, run_steps)
    with open(running_output_dir + '/sub_script', 'w') as f:
        f.write('#!/bin/csh\n#$ -j y\n#$ -o .\n#$ -l s_rt=700:00:00\n#$ -l h_rt=700:00:00\n')
        f.write('#$ -N Static{:}\n\n'.format(restart_run_count))
        f.write('set fromdir = jesu2929@coulson.chem.ox.ac.uk:{:}\n'.format(running_output_dir))
        f.write('set todir   = jesu2929@coulson.chem.ox.ac.uk:{:}\n'.format(running_output_dir))

        f.write('setenv WORK /tmp/$USER/$JOB_ID\n\n')

        f.write('mkdir -p $WORK\n')
        f.write('scp -r $fromdir\* $WORK\n')

        f.write('cd $WORK\n')
        f.write('./dip_10000.x\n')
        f.write('scp *.dat *.log $fromdir\n')
        f.write('scp * $todir\n')
        f.write('rm -Rf $WORK\n')
    if dry_run == False:
        sub_string = str(running_output_dir + 'sub_script')

        subprocess.call(['scp', sub_string, 'jesu2929@dirac.chem.ox.ac.uk:~/From_Coulson/.'])
        subprocess.call(['ssh', 'jesu2929@dirac.chem.ox.ac.uk', 'qsub', 'From_Coulson/sub_script'])





def modify_harmlin(dir, rep, fuck):
    copy(dir+'harm_au.inpt', dir+'harm_au_original.inpt')
    with open(dir+'harm_au_original.inpt', 'r') as input_file:
        with open(dir+'harm_au.inpt', 'w') as output_file:
            f = input_file.readline()
            output_file.write(f)

            for i in range(6):
                output_file.write(input_file.readline())

            if rep==True:
                output_file.write('\n0.5d0     u0 (11)  \n'+
                                  '{:.4f}d0    cut (11) \n'.format(float(fuck*4.9624))+
                                    '0.0d0     u0 (12)  \n'+
                                    '0.0d0     cut (12) \n'+
                                    '0.5d0     u0 (22)  \n'+
                                    '{:.4f}d0     cut (22) \n'.format(float(fuck*4.9624)))
            else:
                output_file.write('\n0.0d0     u0 (11)  \n'+
                                  '{:.4f}d0    cut (11) \n'.format(float(fuck*4.9624))+
                                    '0.0d0     u0 (12)  \n'+
                                    '0.0d0     cut (12) \n'+
                                    '0.0d0     u0 (22)  \n'+
                                    '{:.4f}d0     cut (22) \n'.format(float(fuck*4.9624)))
    return

def modify_runtime_harmlin(dir,steps):
    copy(dir+'runtime_licl.inpt', dir+'runtime_licl_original.inpt')
    with open(dir+'runtime_licl_original.inpt', 'r') as input_file:
        with open(dir+'runtime_licl.inpt', 'w') as output_file:
            output_file.write('{:<15}Number Steps in the run\n'.format(steps))
            skip = input_file.readline()
            for i in range(17):
                output_file.write(input_file.readline())
            output_file.write(input_file.readline())
            for i in range(53-19):
                output_file.write(input_file.readline())
    return

def modify_runtime_static(dir, steps):
    #copy(dir+'runtime_licl.inpt', dir+'runtime_licl_original.inpt')
    with open(dir+'runtime_licl_original.inpt', 'r') as input_file:
        with open(dir+'runtime.inpt', 'w') as output_file:
            output_file.write('{:<15}Number Steps in the run\n'.format(steps))
            skip = input_file.readline()
            for i in range(17):
                output_file.write(input_file.readline())
            output_file.write(input_file.readline())
            for i in range(53-19):
                output_file.write(input_file.readline())
    return

def scale_old_crds_harmpairs(dir, area):
    shutil.move(dir+'crys.crds', dir+'unscaled_crys.crds')
    with open(dir+'unscaled_crys.crds', 'r') as filein:
        crds = np.genfromtxt(filein)
    with open(dir+'dimensions.dat','r') as filein:
        dim = np.genfromtxt(filein)
    dim = np.array([dim[0],dim[1],200])
    crds = np.multiply(crds, 1.0/0.52917721090380)
    dim = np.multiply(dim, 1.0/0.52917721090380)
    crds = np.multiply(crds, np.sqrt(area))
    dim = np.multiply(dim, np.sqrt(area))
    print('Dimensions : ', dim) 
    n_atoms = crds.shape[0]
    nSi = int(n_atoms/3)
    nO  = int(n_atoms-nSi)
    with open(dir+'crys.crds', 'w') as fileout:
        fileout.write('T\nF\nF\nF\nF\n')
        for i in range(nSi, n_atoms):
            fileout.write('{:<26}{:<26}{:<26}\n'.format(crds[i,0], crds[i,1], crds[i,2]))
        for i in range(0,nSi):
            fileout.write('{:<26}{:<26}{:<26}\n'.format(crds[i,0], crds[i,1], crds[i,2]))
        fileout.write('{:<26}{:<26}{:<26}\n'.format(1.0,0.0,0.0))
        fileout.write('{:<26}{:<26}{:<26}\n'.format(0.0,1.0,0.0))
        fileout.write('{:<26}{:<26}{:<26}\n'.format(0.0,0.0,1.0))
        fileout.write('{:<26}\n'.format(dim[0]))
        fileout.write('{:<26}\n'.format(dim[1]))
        fileout.write('{:<26}\n'.format(dim[2]))
    
    
    shutil.move(dir+'harmpairs.dat', dir+'unscaled_harmpairs.dat')
    with open(dir+'unscaled_harmpairs.dat', 'r') as filein:
        harmpairs = np.genfromtxt(filein, skip_header=1)
    no_pairs = harmpairs.shape[0]
    with open(dir+'harmpairs.dat', 'w') as fileout:
        fileout.write('{:}\n'.format(no_pairs))
        for i in range(no_pairs):
            for j in range(2):
                if harmpairs[i,j]<=nSi:
                    harmpairs[i,j] += nO
                else:
                    harmpairs[i,j] -= nSi
            fileout.write('{:<26}{:<26}\n'.format(int(harmpairs[i,0]), int(harmpairs[i,1])))
    
    return

def scale_crds_harmpairs(dir, area):
    shutil.move(dir+'/crys.crds', dir+'/unscaled_crys.crds')
    with open(dir+'/unscaled_crys.crds', 'r') as filein:
        crds = np.genfromtxt(filein, skip_header=5, skip_footer=6)
        dim = np.genfromtxt(filein, skip_header=8+crds.shape[0])
    crds = np.multiply(crds, 1.0/0.52917721090380)
    dim = np.multiply(dim, 1.0/0.52917721090380)
    crds = np.multiply(crds, np.sqrt(area))
    dim = np.multiply(dim, np.sqrt(area))

    n_atoms = crds.shape[0]
    nSi = int(n_atoms/3)
    nO  = int(n_atoms-nSi)
    with open(dir+'/crys.crds', 'w') as fileout:
        fileout.write('T\nF\nF\nF\nF\n')
        for i in range(nSi, n_atoms):
            fileout.write('{:<26}{:<26}{:<26}\n'.format(crds[i,0], crds[i,1], crds[i,2]))
        for i in range(0,nSi):
            fileout.write('{:<26}{:<26}{:<26}\n'.format(crds[i,0], crds[i,1], crds[i,2]))
        fileout.write('{:<26}{:<26}{:<26}\n'.format(1.0,0.0,0.0))
        fileout.write('{:<26}{:<26}{:<26}\n'.format(0.0,1.0,0.0))
        fileout.write('{:<26}{:<26}{:<26}\n'.format(0.0,0.0,1.0))
        fileout.write('{:<26}\n'.format(dim[0]))
        fileout.write('{:<26}\n'.format(dim[1]))
        fileout.write('{:<26}\n'.format(dim[2]))


    shutil.move(dir+'/harmpairs.dat', dir+'/unscaled_hermpairs.dat')
    with open(dir+'/unscaled_harmpairs.dat', 'r') as filein:
        harmpairs = np.genfromtxt(filein, skip_header=1)
    no_pairs = harmpairs.shape[0]
    with open(dir+'/harmpairs.dat', 'w') as fileout:
        fileout.write('{:}\n'.format(no_pairs))
        for i in range(no_pairs):
            for j in range(2):
                if harmpairs[i,j]<=nSi:
                    harmpairs[i,j] += nO
                else:
                    harmpairs[i,j] -= nSi
            fileout.write('{:<26}{:<26}\n'.format(int(harmpairs[i,0]), int(harmpairs[i,1])))

    return

import_aux()


#    elif input_string=='alpha':
#
#        options['input_location'] = ring_size_string+'/'
#        options['input_sublocations'] = []
#        with open('p6_alpha_details.out') as f:
#            array = np.genfromtxt(f)
#        p6_allowed      = np.argwhere((array[:,2]>=(p6_target-0.5*p6_lenience)) & (array[:,2]<=(p6_target+0.5*p6_lenience)))
#        alpha_allowed   = np.argwhere((array[:,3]>=(alpha_target-0.5*alpha_lenience)) & (array[:,3]<=(alpha_target+0.5*alpha_lenience)))
#
#
#        options['output_location'] = options['output_prefix']+'/alpha_{:}/p6_{:}/'.format(alpha_target, p6_target)
#        if not os.path.isdir(options['output_prefix']+'/alpha_{:}'.format(alpha_target)):
#            os.mkdir(options['output_prefix']+'/alpha_{:}'.format(alpha_target))
#        if not os.path.isdir(options['output_prefix']+'/alpha_{:}/p6_{:}/'.format(alpha_target,p6_target)):
#            os.mkdir(options['output_prefix']+'/alpha_{:}/p6_{:}/'.format(alpha_target,p6_target))
#
#        for vals in p6_allowed:
#            if vals in alpha_allowed:
#                options['file_sublocations'].append('T{:}/t{:}_s{:}_same'.format(array[vals,0], array[vals,0], array[vals,1]))
#
#                if not os.path.isdir(options['output_prefix'] + '/alpha_{:}/p6_{:}/t{:}_s{:}/'.format(alpha_target, p6_target, array[vals,0], array[vals,1])):
#                    os.mkdir(options['output_prefix'] + '/alpha_{:}/p6_{:}/t{:}_s{:}/'.format(alpha_target, p6_target, array[vals,0], array[vals,1]))
#
#                for i in range(1+int((upper_A-lower_A)/step_A)):
#                    area = (lower_A+i*step_A)
#                    if not os.path.isdir(options['output_prefix'] + '/alpha_{:}/p6_{:}/t{:}_s{:}/{:}'.format(alpha_target, p6_target, array[vals,0], array[vals,1], int(1000*area))):
#                        os.mkdir(options['output_prefix'] + '/alpha_{:}/p6_{:}/t{:}_s{:}/{:}'.format(alpha_target, p6_target, array[vals,0], array[vals,1], int(1000*area)))
#
#
#                    shutil.copytree(options['input_location'] + 'T{:}/t{:}_s{:}_same'.format(array[vals,0], array[vals,0], array[vals,1]), options['output_prefix'] + '/alpha_{:}/p6_{:}/t{:}_s{:}/{:}'.format(alpha_target, p6_target, array[vals,0], array[vals,1], int(1000*area)))
#                    scale_crds(options['output_prefix'] + '/alpha_{:}/p6_{:}/t{:}_s{:}/{:}'.format(alpha_target, p6_target, array[vals,0], array[vals,1], int(1000*area)))



#print('Number of arguments:', len(sys.argv), 'arguments.')
#print('Argument List:', str(sys.argv))

def make_extract_restart(working_dir, area_value):
    with open(working_dir+'/area_{:}/extract_restart.sh'.format(area_value), 'w') as filename:
        filename.write('#!/bin/bash\n')
        filename.write('\n')

        filename.write('if [ -d "output" ]\n')
        filename.write('then\n')
        filename.write('    rm -r output\n')
        filename.write('fi\n')

        filename.write('if [ -d "output_restart" ]\n')
        filename.write('then\n')
        filename.write('    rm -r output_restart\n')
        filename.write('fi\n')

        filename.write('\n')
        filename.write('mkdir output_restart\n')
        filename.write('\n')
        filename.write('\n')
        filename.write('mv fort* output_restart/\n')
        filename.write('mv *.out* output_restart/\n')
        filename.write('mv header.hdr output_restart/\n')
        filename.write('#mv testout.rst output_restart/\n')

def make_extract(working_dir, area_value):
    with open(working_dir+'/area_{:}/extract.sh'.format(area_value), 'w') as filename:
        filename.write('#!/bin/bash\n')
        filename.write('\n')

        filename.write('if [ -d "output" ]\n')
        filename.write('then\n')
        filename.write('    rm -r output\n')
        filename.write('fi\n')

        filename.write('\n')
        filename.write('mkdir output_restart\n')
        filename.write('\n')
        filename.write('\n')
        filename.write('mv fort* output_restart/\n')
        filename.write('mv *.out* output_restart/\n')
        filename.write('mv header.hdr output_restart/\n')
        filename.write('mv testout.rst output_restart/\n')


def make_subscript(working_dir, area_value):
    with open(working_dir+'/area_{:}/sub_script'.format(area_value), 'w') as filename:
        filename.write('!/bin/csh \n')
        filename.write('#$ -j y\n')
        filename.write('#$ -o . \n')
        filename.write('#$ -l s_rt=1600:00:00\n')
        filename.write('#$ -l h_rt=1600:00:00\n')
#        filename.write('#$ -pe smp 16')
        filename.write('#$ -N area_{:}\n'.format(area_value))
        filename.write('')
        filename.write('set fromdir = jesu2929@coulson.chem.ox.ac.uk:{:}\n'.format(working_dir+'/area_{:}/'.format(area_value)))
        filename.write('set todir =   jesu2929@coulson.chem.ox.ac.uk:{:}\n'.format(working_dir+'/area_{:}/'.format(area_value)))
        filename.write('')
        filename.write('setenv WORK /tmp/$USER/code_cycle/job_data/$JOB_ID\n')
        filename.write('')
        filename.write('mkdir -p $WORK\n')
        filename.write('scp -r $fromdir\* $WORK\n')

        filename.write('cd $WORK\n')
        filename.write('{:}\n'.format(working_dir+'/area_{:}/dip_10000.x'.format(area_value)))
        filename.write('scp * $todir\n')
        filename.write('rm -Rf $WORK\n')



def scale_crds_harmpairs(dir, area):
    shutil.move(dir+'/crys.crds', dir+'/unscaled_crys.crds')
    with open(dir+'/unscaled_cryrs.crds', 'r') as filein:
        crds = np.genfromtxt(filein, skip_header=5, skip_footer=6)
        dim = np.genfromtxt(filein, skip_header=8+crds.shape[0])
    crds = np.multiply(crds, 1.0/0.52917721090380)
    dim = np.multiply(dim, 1.0/0.52917721090380)
    crds = np.multiply(crds, np.sqrt(area))
    dim = np.multiply(dim, np.sqrt(area))

    n_atoms = crds.shape[0]
    nSi = int(n_atoms/3)
    nO  = int(n_atoms-nSi)
    with open(dir+'/crys.crds', 'w') as fileout:
        fileout.write('T\nF\nF\nF\nF\n')
        for i in range(nSi, n_atoms):
            fileout.write('{:<26}{:<26}{:<26}\n'.format(crds[i,0], crds[i,1], crds[i,2]))
        for i in range(0,nSi):
            fileout.write('{:<26}{:<26}{:<26}\n'.format(crds[i,0], crds[i,1], crds[i,2]))
        fileout.write('{:<26}{:<26}{:<26}\n'.format(1.0,0.0,0.0))
        fileout.write('{:<26}{:<26}{:<26}\n'.format(0.0,1.0,0.0))
        fileout.write('{:<26}{:<26}{:<26}\n'.format(0.0,0.0,1.0))
        fileout.write('{:<26}'.format(dim[0]))
        fileout.write('{:<26}'.format(dim[1]))
        fileout.write('{:<26}'.format(dim[2]))


    shutil.move(dir+'/harmpairs.dat', dir+'/unscaled_hermpairs.dat')
    with open(dir+'/unscaled_harmpairs.dat', 'r') as filein:
        harmpairs = np.genfromtxt(filein, skip_header=1)
    no_pairs = harmpairs.shape[0]
    with open(dir+'/harmpairs.dat', 'w') as fileout:
        fileout.write('{:}\n'.format(no_pairs))
        for i in range(no_pairs):
            for j in range(2):
                if harmpairs[i,j]<=nSi:
                    harmpairs[i,j] += nO
                else:
                    harmpairs[i,j] -= nSi
            fileout.write('{:<26}{:<26}\n'.format(int(harmpairs[i,0]), int(harmpairs[i,1])))

    return


def main(options):
    if (options["restart"]==False):
        folders = glob.glob('area_*', recursive=True)
        for folder in folders:
            try:
                shutil.rmtree(folder)
            except OSError as e:
                print("Error : %s" % (e.strerror))


    working_dir = os.getcwd()


#    with open('test_opt_pe.dat', 'r') as f:
#        array = np.genfromtxt(f)

    with open("crys.crds", "r") as f:
        array = np.genfromtxt(f, skip_header=5, skip_footer=6)



    for i in range(array.shape[0]):
        area_i = array[i,0]
        area_label = str(int(area_i*1000))
        print(area_label)

        pbx = array[i,1]
        pby = array[i,2]
        pbz = array[i,3]
        print(pbx,'  ', pby, '  ', pbz)

        factor = 1.0/0.52917721090380
        print(pbx*factor,'  ', pby*factor, '  ', pbz*factor)
        if (options["restart"]==True):
            os.chdir("area_{}".format(area_label))
            make_extract_restart(working_dir, area_label)
            subprocess.call(['./extract_restart.sh'])

            if (options["server"]==False):
                subprocess.call(['nq', './dip_11000.x'])
            else:
                make_subscript(working_dir, area_label)
                subprocess.call(['scp', 'sub_script', 'jesu2929@dirac.chem.ox.ac.uk:~/From_Coulson/.'])
                print('Copied')
                time.sleep(1000)
                subprocess.call(['ssh', 'dirac', 'qsub', 'From_Coulson/sub_script'])                                                    #
            os.chdir('../')

        if (options["restart"]==False):
            os.mkdir("area_{}".format(area_label))

            shutil.copy('test_area_{}_sample_0_crys.crds'.format(i), 'area_{}/'.format(area_label))
            with open('test_area_{}_sample_0_crys.crds'.format(i), 'r') as f:
                crds = np.genfromtxt(f)

            print(max(crds[:,0]))
            print(max(crds[:,1]))
            print('\n')
            print(pbx*factor)
            print(pby*factor)

            with open("area_{}/crys.crds".format(area_label), 'w') as f:
                f.write('T\nF\nF\nF\nF\n')
                for i in range(crds.shape[0]):
                    f.write('{:<20}{:<20}{:<20}\n'.format(crds[i,0], crds[i,1], crds[i,2]))

            for files in ['cf.inpt','a.x','a.f','dip_10000.x','crystal_cell.inpt','dip_11000.x','tube_19_0_15.crds','ts_sio2.inpt','runtime.inpt','q','lj.inpt','lightcf_licl.inpt','filenames.inpt','extract.sh','dynmat.inpt','test_opt_pe.dat','split_to_areas.py']:
                shutil.copy(files, 'area_{}/'.format(area_label))

            with open('area_{}/crystal_cell.inpt'.format(area_label), 'w') as g:


                g.write('   1.0d0  0.0d0  0.0d0\n')
                g.write('   0.0d0  1.0d0  0.0d0     Cell matrix.\n')
                g.write('   0.0d0  0.0d0  1.0d0\n')
                g.write('{:<25} L_x\n'.format(pbx*factor))
                g.write('{:<25} L_y\n'.format(pby*factor))
                g.write('{:<25} L_z\n'.format(200.0))
                g.write('           20       Number of unit cells in x direction (starting from a crystal only)\n')
                g.write('           12       Number of unit cells in y direction (starting from a crystal only)\n')
                g.write('           1        Number of unit cells in z direction (starting from a crystal only)\n')
                g.write('           2        Number of type I particles in unit cell.\n')
                g.write('x_tet.mat\n')
                g.write('           2        Number of type II particles in unit cell.\n')
                g.write('m_tet.mat\n')
                g.write('   4.6d0            Unit cell length - x direction.\n')
                g.write('   8.0d0            Unit cell length - x direction.\n')
                g.write('   100.0d0          Unit cell length - x direction.\n')
            with open('area_{:}/runtime.inpt'.format(area_label), 'w') as filename:
                print(run_steps)
                time.sleep(1000)
                filename.write('{:}         Number of steps in the run.\n'.format(run_steps))
                filename.write('0.0001d0    Translational temperature.\n')
                filename.write('151.315     Relative molecular mass (gmol^-1)\n')
                filename.write('800         Number of ionic molecular units.\n')
                filename.write('2           Number of ionic species.\n')
                filename.write('1600        Number of species of type 1.\n')
                filename.write('-1.38257d0  Permanent charge on type 1 ions.\n')
                filename.write('16.0d0      Atomic mass of type 1.\n')
                filename.write('.true.      Polarizable ?\n')
                filename.write('800         Number of species of type 2.\n')
                filename.write('2.76514d0   Permanent charge on type 2 ions.\n')
                filename.write('28.08d0     Atomic mass of type 2.\n')
                filename.write('.false.     Polarizable ?\n')
                filename.write('25.0d0      Timestep (a.u.).\n')
                filename.write('rim         Type of run (rim,dippim,quadpim)\n')
                filename.write('ts          Type of run (epp,cim,aim)\n')
                filename.write('.false.     Steepest descent minimisation ? (if .true. require time step + annealing steps + T/F for anneal radii only).\n')
                filename.write('.false.     Conjugate gradient minimisation?\n')
                filename.write('.true.      Restart ? (if .true. then next line is \'inpt\' or \'rst\' if cell lengths are to be read from the \'crystal_cell\' or restart files).\n')
                filename.write('inpt\n')
                filename.write('.false.     Set up velocities?\n')
                filename.write('.false.     Velocity rescale prior to main run ?\n')
                filename.write('.false.     Random displacement of ions.\n')
                filename.write('.true.      Move ions?\n')
                filename.write('.false.     Do a dynamical matrix calculation?\n')
                filename.write('.true.      Relax input structure? (.true. will \'kill\' the velocities at a local k.e. maximum - steepest descent minimisation).\n')
                filename.write('1           Number of steps inbetween periodic output (energies).\n')
                filename.write('100         Number of steps inbetween periodic output (velocities etc).\n')
                filename.write('1           Number of steps inbetween periodic output (frictions etc).\n')
                filename.write('1           Number of steps inbetween periodic output (pressure etc).\n')
                filename.write('10          Number of steps inbetween rdf call in main loop.\n')
                filename.write('10          Number of steps inbetween S(k) call in main loop\n')
                filename.write('1           Number of ions to monitor.\n')
                filename.write('1           Ion number to monitor. (1)\n')
                filename.write('5.60d0      eta = <x>/boxlen.\n')
                filename.write('45.0d0      rcut (au).\n')
                filename.write('1.0d-7      convergence parameter.\n')
                filename.write('.false.     Nose-Hoover thermostat? (if true then enter a relaxation time)\n')
                filename.write('.false.     Periodic rescale of temperature?\n')
                filename.write('.false.     Periodic re-annealing?\n')
                filename.write('.false.     Isotropic barostat?\n')
                filename.write('.false.     Anisotropic barostat?\n')
                filename.write('.false.     Orthorhombic cell?\n')
                filename.write('.false.     T ramping ?\n')
                filename.write('.false.     p ramping?\n')
                filename.write('.false.     EFG calc. ?\n')
                filename.write('.false.     Box ramping?\n')
                filename.write('.false.     Box ramping?\n')





            os.chdir('area_{}'.format(area_label))
            subprocess.call(['nq', './dip_11000.x'])
            os.chdir('../')
            #    os.chdir('area_{}'.format(area_label))
            #    subprocess.call(['./dip_11000.x'])
            #    time.sleep(1000)
            #    os.chdir('../')
#options = get_options()
#main(options)


