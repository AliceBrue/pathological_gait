import numpy as np


def check_std(par_folder, optim_par_file, init_par_file, check_offsets=True):
    """
    Compute, sort and write in a txt file the difference between parameters of two par files
    INPUTS: - par_folder: name of the folder containing par files
            - optim_par_file: optimisation par file name
            - init_par_file: init par file name
            - check_offsets: if False, offset parameters are not taken into account, default=True
    OUTPUTS: - txt file listing sorted differences
    """

    #: Extract values from par files
    optim_par = open(par_folder+"/"+optim_par_file+".par", "r")
    init_par = open(par_folder+"/"+init_par_file+".par", "r")
    optim_lines = optim_par.readlines()
    init_lines = init_par.readlines()

    #: List init parameters (if different optim parameters)
    init_params = []
    for i in range(len(init_lines)):
        init_params += [init_lines[i].split()[0]]

    #: List differences between parameters
    param = []
    diff = []
    init_std = []
    j = 0
    for i in range(len(optim_lines)):
        optim_param = optim_lines[i].split()[0]
        init_param = init_lines[j].split()[0]
        if optim_param == init_param:
            param += [optim_param]
            diff += [abs(float(optim_lines[i].split()[1])-float(init_lines[j].split()[1]))]
            init_std += [float(init_lines[j].split()[3])]
            j += 1
        else:
            #: If substracted optim paramaters
            if optim_param in init_params[j:]:
                ind = init_params[j:].index(optim_param)
                j += 1
                for k in range(ind):
                    j += 1

    #: Sort differences
    n_param = len(param)
    diff = np.array(diff)
    sort_ind = np.flip(np.argsort(diff))

    #: Write results file
    diff_file = open(par_folder+"/diff_"+optim_par_file+"_"+init_par_file+".txt", "w+")
    diff_file.write("param diff init_std %init_std \n")
    diff_file.write("mean " + str(np.sum(diff)/n_param) + " " + str(np.sum(init_std)/n_param) + " " +
                    str((np.sum(diff)/n_param)*100/(np.sum(init_std)/n_param)) + "\n")
    for i in range(n_param):
        if not check_offsets:
            if param[i].split(".")[-1] != "offset":
                diff_file.write(param[sort_ind[i]] + " " + str(diff[sort_ind[i]]) + " " + str(init_std[sort_ind[i]]) +
                                " " + str(diff[sort_ind[i]]*100/init_std[sort_ind[i]]) + "\n")
        else:
            diff_file.write(param[sort_ind[i]] + " " + str(diff[sort_ind[i]]) + " " + str(init_std[sort_ind[i]]) + " " +
                            str(diff[sort_ind[i]]*100/init_std[sort_ind[i]]) + "\n")
