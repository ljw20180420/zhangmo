# Program to perform TAD analysis:
# - Calculate the DI, HMM states and topological domains coordinates.
# - Plot the observed DI and true DI (Hidden Markov Model).

# Usage: python2.7 HiCtool_TAD_analysis.py [-h] [options]
# Options:
#  -h, --help                show this help message and exit
#  --action                  Action to perform: full_tad_analysis, plot_chromosome_DI. 
#  -i INPUT_FILE             Input contact matrix file if action is "full_tad_analysis" or DI values if action is "plot_chromosome_DI".
#  -c CHROMSIZES_PATH        Path to the folder chromSizes with trailing slash at the end.
#  -s SPECIES                Species. It has to be one of those present under the chromSizes path.
#  --isGlobal                Insert 1 if the input matrix is a global matrix, 0 otherwise.  
#  --tab_sep                 Insert 1 if the input matrix is in a tab separated format, 0 if it is in compressed format.
#  --chr                     If action is "full_tad_analysis": chromosome or list of chromosomes between square brackets to select specific maps for the analysis. If action is "plot_chromosome_DI" insert a single chromosome to plot the DI values.
#  --data_type               Data type to label your data, example: observed, normalized, etc.
#  --full_chromosome         Insert 1 to plot DI and HMM states for the entire chromosome, 0 otherwise.
#  --coord                   List of two integers with start and end coordinates to plot the DI values and HMM values.
#  --input_file_hmm          Input HMM states file if action is "plot_chromosome_DI" to plot also the HMM states.
#  --plot_legend             If action is "plot_chromosome_DI", insert 1 to plot the legend, 0 otherwise.
#  --plot_grid               If action is "plot_chromosome_DI", insert 1 to plot the grid, 0 otherwise.

from optparse import OptionParser
import numpy as np
import os
import os.path
from os import path

parameters = {'action': None,
              'input_file': None,
              'chromSizes_path': None,
              'isGlobal': None,
              'tab_sep': None,
              'chr': None,
              'species': None,
              'data_type': None,
              'full_chromosome': None,
              'coord': None,
              'input_file_hmm': None,
              'plot_legend': None,
              'plot_grid': None
              }


def save_list(a_list, output_file):
    """
    Save a list in a txt file.
    Arguments:
        a_list (obj): name of the list to save.
        output_file (str): output file name in txt format.
    Output: 
        txt file containing the saved list.       
    """
    with open (output_file,'w') as fout:
        n = len(a_list)
        for i in range(n):
            fout.write('%s\n' %a_list[i])

def save_matrix_rectangular(a_matrix, output_file):
    """
    Save an inter-chromosomal contact matrix in the HiCtool compressed format to txt file.
    1) Data are reshaped to form a vector.
    2) All the consecutive zeros are replaced with a "0" followed by the
    number of times zeros are repeated consecutively.
    3) Data are saved to a txt file.
    Arguments:
        a_matrix (numpy matrix): input contact matrix to be saved
        output_file (str): output file name in txt format
    Output:
        txt file containing the formatted data
    """
    import numpy as np
    n_row = np.shape(a_matrix)[0]
    n_col = np.shape(a_matrix)[1]
    vect = np.reshape(a_matrix,[1,n_row*n_col]).tolist()[0]
    with open (output_file,'w') as fout:
        k = len(vect)
        i = 0
        count = 0
        flag = False # flag to set if the end of the vector has been reached
        while i < k and flag == False:
            if vect[i] == 0:
                count+=1
                if (i+count == k):
                    w_out = str(0) + str(count)
                    fout.write('%s\n' %w_out)
                    flag = True
                    break
                while vect[i+count] == 0 and flag == False:
                    count+=1
                    if (i+count == k):
                        w_out = str(0) + str(count)
                        fout.write('%s\n' %w_out)
                        flag = True
                        break
                if flag == False:
                    w_out = str(0) + str(count)
                    fout.write('%s\n' %w_out)
                    i+=count
                    count = 0
            else:
                fout.write('%s\n' %vect[i])
                i+=1 

def save_matrix(a_matrix, output_file):
    """
    Save an intra-chromosomal contact matrix in the HiCtool compressed format to txt file.
    1) The upper-triangular part of the matrix is selected (including the
    diagonal).
    2) Data are reshaped to form a vector.
    3) All the consecutive zeros are replaced with a "0" followed by the
    number of times zeros are repeated consecutively.
    4) Data are saved to a txt file.
    Arguments:
        a_matrix (numpy matrix): input contact matrix to be saved
        output_file (str): output file name in txt format
    Output:
        txt file containing the formatted data
    """
    import numpy as np
    n = len(a_matrix)
    iu = np.triu_indices(n)
    vect = a_matrix[iu].tolist()
    with open (output_file,'w') as fout:
        k = len(vect)
        i = 0
        count = 0
        flag = False # flag to set if the end of the vector has been reached
        while i < k and flag == False:
            if vect[i] == 0:
                count+=1
                if (i+count == k):
                    w_out = str(0) + str(count)
                    fout.write('%s\n' %w_out)
                    flag = True
                    break
                while vect[i+count] == 0 and flag == False:
                    count+=1
                    if (i+count == k):
                        w_out = str(0) + str(count)
                        fout.write('%s\n' %w_out)
                        flag = True
                        break
                if flag == False:
                    w_out = str(0) + str(count)
                    fout.write('%s\n' %w_out)
                    i+=count
                    count = 0
            else:
                fout.write('%s\n' %vect[i])
                i+=1     

def load_matrix(input_file):
    """
    Load an HiCtool compressed square (and symmetric) contact matrix from a txt file and parse it.
    Arguments:
        input_file (str): input file name in txt format (generated by the function 
        "save_matrix").
    Return: numpy array containing the parsed values stored in the input txt file to build a contact matrix.      
    """
    import numpy as np    
    
    print("Loading " + input_file + "...")
    with open (input_file,'r') as infile:
        matrix_vect = []        
        for i in infile:
            if i[0] == "0" and i[1] != ".":
                for k in range(int(i[1:-1])):
                    matrix_vect.append(0)
            else:
                j = i[:-1]            
                matrix_vect.append(float(j))
  
    k = len(matrix_vect)
    matrix_size = int((-1+np.sqrt(1+8*k))/2)
    
    iu = np.triu_indices(matrix_size)
    output_matrix_1 = np.zeros((matrix_size,matrix_size)) # upper triangular plus the diagonal
    output_matrix_1[iu] = matrix_vect
    
    diag_matrix = np.diag(np.diag(output_matrix_1)) # diagonal
    output_matrix_2 = np.transpose(output_matrix_1) # lower triangular plus the diagonal
    output_matrix = output_matrix_1 + output_matrix_2 - diag_matrix
    print("Done!")
    return output_matrix

    
def save_matrix_tab(input_matrix, output_filename):
    """
    Save a contact matrix in a txt file in a tab separated format. Columns are
    separated by tabs, rows are in different lines.
    Arguments:
        input_matrix (numpy matrix): input contact matrix to be saved
        output_filename (str): output file name in txt format
    Output:
        txt file containing the tab separated data
    """
    with open (output_filename, 'w') as f:
            for i in range(len(input_matrix)):
                row = [str(j) for j in input_matrix[i]]
                f.write('\t'.join(row) + '\n')

def load_matrix_tab(input_file):
    """
    Load a contact matrix saved in a tab separated format using the function
    "save_matrix_tab".
    Arguments:
        input_file (str): input contact matrix to be loaded.
    Return: numpy array containing the parsed values stored in the input tab separated txt file to build a contact matrix.
    """
    import numpy as np
    
    print("Loading " + input_file + "...")
    with open (input_file, 'r') as infile:
        lines = infile.readlines()
        temp = []
        for line in lines:
            row = [float(i) for i in line.strip().split('\t')]
            temp.append(row)
            
        output_matrix = np.array(temp)
    print("Done!")
    return output_matrix


    
def load_DI_values(input_file):
    """
    Load a DI txt file generated with "calculate_chromosome_DI".
    Arguments:
        input_file (str): input file name in txt format.
    Return: 
        List of the DI values.        
    """
    import numpy as np
    
    fp = open(input_file,'r+')
    lines = fp.read().split('\n')
    lines = lines[:-1]
    
    return np.nan_to_num(np.array([float(line) for line in lines])).tolist()


def extract_single_map(input_global_matrix,
                       tab_sep,
                       chr_row,
                       chr_col,
                       species='hg38',
                       bin_size=1000000,
                       data_type='observed',
                       save_output=True,
                       save_tab=False):
    """
    Extract a single contact matrix for a pair of chromosomes from the global matrix (all-by-all chromosomes).
    Arguments:
        input_global_matrix (object | str): global contact matrix. This can be passed either as
        an object of the workspace or a string of the filename saved to file.
        tab_sep (bool): if "input_global_matrix" is passed with a filename, then this boolean 
        tells if the global matrix was saved in tab separated format (True) or not (False).
        chr_row (str): chromosome in the rows of the output contact matrix.
        chr_col (str): chromosome in the columns of the output contact matrix. If chr_col is 
        equal to chr_row then the intra-chromosomal map is extracted.
        species (str): species label in string format.
        bin_size (int): bin size in bp of the contact matrix.
        data_type (str): which kind of data type you are extracting: "observed" or "normalized".
        save_output (bool): if True, save the contact matrix in HiCtool compressed txt file.
        save_tab (bool): if True, save the contact matrix in tab separated format.
    Return: 
        Contact matrix in numpy array format.
    Outputs:
        Txt file with the contact matrix in HiCtool compressed format if "save_output=True".
        Txt file with the contact matrix in tab separated format if "save_tab=True".
    """            
    chromosomes = open(parameters['chromSizes_path'] + parameters['species'] + '.chrom.sizes', 'r')
    chromosomes_list = []
    chr_dim = []
    d_chr_dim = {}
    while True:
        try:
            line2list = next(chromosomes).split('\n')[0].split('\t')
            chromosomes_list.append(line2list[0])
            chr_dim.append(int(line2list[1])/bin_size)
            d_chr_dim[line2list[0]] = int(line2list[1])/bin_size
        except StopIteration:
            break
    
    d_chr_dim_inc = {}
    k=1
    for i in chromosomes_list:
        d_chr_dim_inc[i] = sum(chr_dim[:k])
        k+=1
    
    if isinstance(input_global_matrix,str):
        if tab_sep == False:
            full_matrix = load_matrix(input_global_matrix)
        else:
            full_matrix = load_matrix_tab(input_global_matrix)
    else:
        full_matrix = input_global_matrix
    
        d_chr_dim_inc = {}
    k=1
    for i in chromosomes_list:
        d_chr_dim_inc[i] = sum(chr_dim[:k])
        k+=1
    
    if isinstance(input_global_matrix,str):
        if tab_sep == False:
            full_matrix = load_matrix(input_global_matrix)
        else:
            full_matrix = load_matrix_tab(input_global_matrix)
    else:
        full_matrix = input_global_matrix
    
    if chr_row == '1':
        row_start = 0
    else:
        row_start = d_chr_dim_inc[chromosomes_list[chromosomes_list.index(chr_row)-1]]
    row_end = row_start + d_chr_dim[chr_row]
    
    if chr_col == '1':
        col_start = 0
    else:
        col_start = d_chr_dim_inc[chromosomes_list[chromosomes_list.index(chr_col)-1]]
    col_end = col_start + d_chr_dim[chr_col]
    
    output_matrix = full_matrix[row_start:row_end,col_start:col_end]
    
    if chr_row == chr_col:
        if bin_size >= 1000000:
            bin_size_str = str(bin_size/1000000)
            my_filename = 'HiCtool_' 'chr' + chr_row + '_' + bin_size_str + 'mb_' + data_type + '.txt'
        elif bin_size < 1000000:
            bin_size_str = str(bin_size/1000)
            my_filename = 'HiCtool_' 'chr' + chr_row + '_' + bin_size_str + 'kb_' + data_type + '.txt'
        if save_output == True:
            save_matrix(output_matrix, my_filename)
    else:
        if bin_size >= 1000000:
            bin_size_str = str(bin_size/1000000)
            my_filename = 'HiCtool_' 'chr' + chr_row + '_chr' + chr_col + '_' + bin_size_str + 'mb_' + data_type + '.txt'
        elif bin_size < 1000000:
            bin_size_str = str(bin_size/1000)
            my_filename = 'HiCtool_' 'chr' + chr_row + '_chr' + chr_col + '_' + bin_size_str + 'kb_' + data_type + '.txt'
        if save_output == True:
            save_matrix_rectangular(output_matrix, my_filename)
    
    if save_tab == True:
        save_matrix_tab(output_matrix, my_filename.split('.')[0] + '_tab.txt')
    
    return output_matrix


def calculate_chromosome_DI(input_contact_matrix, 
                            a_chr,
                            isGlobal,
                            tab_sep=False,
                            data_type='normalized',
                            species='hg38',
                            save_file=True,
                            bin_size=40000,
                            output_path="tad_analysis"):
    """
    Function to calculate the DI values for a chromosome and save them 
    in a txt file.
    Arguments:
        input_contact_matrix (str | obj): normalized intra-chromosomal contact matrix at a bin size of 40kb passed as a filename (str)
        or an object. Either a single contact matrix or a global contact matrix can be passed (see following arguments).
        a_chr (str): chromosome number (example for chromosome 1: '1').
        isGlobal (bool): set True if your input matrix is a global matrix (all-by-all chromosomes).
        tab_sep (bool): set True if your input matrix is in a tab separated format. If the matrix is passed as an
        object, this parameter is not taken into consideration.
        species (str): species label in string format.
        save_file (bool): if True, saves the DI values to txt file.
    Returns: List with the DI values.
    Output: Txt file with the DI values if "save_file=True".
    """
    
    import copy
    
    if isGlobal == False:
        if isinstance(input_contact_matrix, str):
            if tab_sep == False:
                contact_matrix = load_matrix(input_contact_matrix)
            else:
                contact_matrix = load_matrix_tab(input_contact_matrix)
        else:
            contact_matrix = copy.deepcopy(input_contact_matrix)
    else:
        contact_matrix = extract_single_map(input_global_matrix=input_contact_matrix,
                                            tab_sep=tab_sep,
                                            chr_row=a_chr,
                                            chr_col=a_chr,
                                            species=species,
                                            bin_size=bin_size,
                                            data_type=data_type,
                                            save_output=False,
                                            save_tab=False)
        
    print("Calculating DI values...")
    n = contact_matrix.shape[0]

    # Calculation of the DI
    DI = [] # list of the DI for each bin
    len_var = 2000000//bin_size # range of upstream or downstream bins to calculate DI

    for locus in range(n): # 'locus' refers to a bin
        if locus < len_var:
            A = sum(contact_matrix[locus][:locus])
            B = sum(contact_matrix[locus][locus+1:locus+len_var+1])
        elif locus >= n-len_var:
            A = sum(contact_matrix[locus][locus-len_var:locus])
            B = sum(contact_matrix[locus][locus+1:])
        else:
            A = sum(contact_matrix[locus][locus-len_var:locus])
            B = sum(contact_matrix[locus][locus+1:locus+len_var+1])
    
        E = (A+B)/2 # expected number of reads
        
        if A==0 and B==0:
            di = 0
            DI.append(di)
        else:
            try:
                di = ((B-A)/(abs(B-A)))*((((A-E)**2)/E)+(((B-E)**2)/E))
            except ZeroDivisionError:
                di = 0
            DI.append(di)
    
    if save_file == True:
        save_list(DI, os.path.join(output_path, 'HiCtool_chr' + a_chr + '_DI.txt'))
    
    print("Done!")
    return DI
    
    
def calculate_chromosome_hmm_states(input_file_DI,
                                    a_chr,
                                    save_file=True,
                                    zero_threshold=0.4,
                                    output_path="tad_analysis"):
    """
    Function to calculate the HMM states (true DI values) for a chromosome and save 
    them in a txt file. It takes DI values as input.
    Arguments:
        input_file_DI (str | obj): txt file of the DI values generated with the function "calculate_chromosome_DI" or
        object with the DI values returned by "calculate_chromosome_DI".
        a_chr (str): chromosome number (example for chromosome 1: '1').
        save_file (bool): if True, saves the hmm states to txt file.
    Returns: List with the HMM states.
    Output: Txt file with the DI values if "save_file=True".
    """
    import numpy as np
    import hmmlearn.hmm as hmm
    
    print("Calculating true DI values...")
    if isinstance(input_file_DI,str):
        A = load_DI_values(input_file_DI)
    else:
        A = input_file_DI   
    
    # Guessed Transition Matrix
    TRANS_GUESS = np.array([[0.4, 0.3, 0.3],
                             [0.3, 0.4, 0.3],
                             [0.3, 0.3, 0.4]])
    
    # Guessed Emission Matrix   
    EMISS_GUESS = np.array([[0.4, 0.3, 0.3],
                             [0.3, 0.4, 0.3],
                             [0.3, 0.3, 0.4]])
    
    # Observed emissions                  
    emissions = []
    
    for i in range(0,len(A)):
        if A[i] >= zero_threshold:
            emissions.append(1)
        elif A[i] <= -zero_threshold:
            emissions.append(2)
        else:
            emissions.append(0)

    # Hidden Markov Model with discrete emissions
    # model = hmm.MultinomialHMM(n_components=3, init_params="")
    model = hmm.CategoricalHMM(n_components=3, n_features=3, init_params="") # the MultinomialHMM of old version is the CategoricalHMM of new version
    model.transmat_ = TRANS_GUESS
    model.emissionprob_ = EMISS_GUESS
    
    input_observations = np.array([emissions]).T
    model.fit(input_observations) # estimate model parameters
    
    # Find most likely state sequence corresponding to input_onservations using the Viterbi Algorithm
    logprob, likelystates_array = model.decode(input_observations, algorithm="viterbi")
    
    likelystates = likelystates_array.tolist()
    
    if save_file == True:
        save_list(likelystates, os.path.join(output_path, "HiCtool_chr" + a_chr + "_hmm_states.txt"))
    
    print("Done!")
    return likelystates


def load_hmm_states(input_file):
    """
    Load an HMM txt file generated with "calculate_chromosome_hmm_states".
    Arguments:
        input_file (str): input file name in txt format.
    Returns: 
        List of the HMM states.    
    """
    fp = open(input_file,'r+')
    lines = fp.read().split('\n')
    lines = lines[:-1]
    
    return [int(line) for line in lines]


def plot_chromosome_DI(input_file_DI, 
                       a_chr,
                       full_chromosome,
                       start_pos=0, 
                       end_pos=0,
                       input_file_hmm=None,
                       species='hg38',
                       plot_legend=True,
                       plot_grid=True,
                       bin_size=40000,
                       output_path="tad_analysis"):
    """
    Function to plot the DI and true DI values for a chromosome.
    Arguments:
        input_file_DI (str | obj): txt file of the DI values generated with the function "calculate_chromosome_DI" or
        object with the DI values returned by "calculate_chromosome_DI".
        a_chr (str): chromosome number (example for chromosome 1: '1').
        full_chromosome (bool): if True, plot the full chromosome "a_chr". In this case "start_pos" and "end_pos" parameters are not considered.
        start_pos (int): start coordinate for the plot in bp.
        end_pos (int): end coordinate for the plot in bp.
        input_file_hmm (str | obj): txt file of the true DI values generated with the function "calculate_chromosome_hmm_states" or
        object with the true DI values returned by "calculate_chromosome_hmm_states".
        species (str): species name (hg38, mm10, etc.).
        plot_legend (bool): if True, plot the legend.
        plot_grid (bool): if True, plot the grid.
    Output:
        Plot saved to pdf file.
    """    
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.use('Agg')

    chromosomes = open(parameters['chromSizes_path'] + parameters['species'] + '.chrom.sizes', 'r')
    chromosomes_list = []
    chr_dim = []
    d_chr_length = {}
    d_chr_dim = {}
    while True:
        try:
            line2list = next(chromosomes).split('\n')[0].split('\t')
            chromosomes_list.append(line2list[0])
            chr_dim.append(int(line2list[1])/bin_size)
            d_chr_length[line2list[0]] = int(line2list[1])
            d_chr_dim[line2list[0]] = int(line2list[1])/bin_size
        except StopIteration:
            break
    
    if full_chromosome == True:
        start_pos = 0
        end_pos = int(round(d_chr_dim[a_chr]))*bin_size
        start_index = 0
        end_index = int(round(d_chr_dim[a_chr])) + 1
    else:
        if end_pos == 0:
            print("ERROR: insert start and end coordinates")
            return
        start_index = int(round(start_pos/bin_size))
        end_index = int(round((end_pos)/bin_size))
    
        if end_pos > int(round(d_chr_dim[a_chr]))*bin_size and end_pos <= d_chr_length[a_chr]:
            end_pos = int(round(d_chr_dim[a_chr]))*bin_size
        elif end_pos > d_chr_length[a_chr]:
            print("ERROR: end coordinate exceeds chromosome dimension")
            return

    if isinstance(input_file_DI,str):
        DI = load_DI_values(input_file_DI)
    else:
        DI = input_file_DI   
        
    DI_part = DI[start_index:end_index]

    x = np.arange(start_pos,end_pos,bin_size)
    width = bin_size/1.5
    pos_DI = np.array(DI_part)
    neg_DI = np.array(DI_part)
    pos_DI[pos_DI <= 0] = np.nan
    neg_DI[neg_DI > 0] = np.nan
    
    if input_file_hmm == None:
        print("Plotting DI values...")

        plt.close("all")
        plt.bar(x, pos_DI, width, color="r", label="Positive DI")
        plt.bar(x, neg_DI, width, color="g", label="Negative DI")
        plt.xlim([x[0]-bin_size*8,x[-1]+bin_size*8])
        plt.ylim([min(DI_part)-25,max(DI_part)+25])
        plt.title("Directionality Index " + species + " [Chr " + a_chr +": " + str(start_pos) + "-" + str(end_pos) + "]")
        plt.xlabel("Base coordinates")
        plt.ylabel("Directionality Index (DI) values")
        plt.grid(plot_grid)
        if plot_legend == True:
            plt.legend(prop={'size': 8})
        plt.savefig(os.path.join(output_path, "HiCtool_chr" + a_chr + "_DI.pdf"), format = 'pdf')
        print("Done!")
    
    else:
        print("Plotting DI and true DI values...")     
        
        if isinstance(input_file_hmm,str):
            likelystates = load_hmm_states(input_file_hmm)
        else:
            likelystates = input_file_hmm
    
        DI_true = []
        for i in range(0,len(likelystates)):
            if likelystates[i] == 1:
                DI_true.append(min(DI_part)-12)
            elif likelystates[i] == 2:
                DI_true.append(min(DI_part)-15)
            else:
                DI_true.append(0)
        
        DI_true_part = DI_true[start_index:end_index]
        
        # Plot
        pos_DI_true = np.array(DI_true_part)
        neg_DI_true = np.array(DI_true_part)
        pos_DI_true[pos_DI_true != min(DI_part)-12] = np.nan
        neg_DI_true[neg_DI_true != min(DI_part)-15] = np.nan
        
        plt.close("all")
        plt.bar(x, pos_DI, width, color="r", label="Positive DI", linewidth = 0.1)
        plt.bar(x, neg_DI, width, color="g", label="Negative DI", linewidth = 0.1)
        plt.plot(x, pos_DI_true, marker=">", color="r", label="Positive true DI")
        plt.plot(x, neg_DI_true, marker="<", color="g", label="Negative true DI")
        plt.xlim([x[0]-bin_size*8,x[-1]+bin_size*8])
        plt.ylim([min(DI_part)-25,max(DI_part)+25])
        plt.title("Directionality Index " + species + " [Chr " + a_chr +": " + str(start_pos) + "-" + str(end_pos) + "]")
        plt.xlabel("Base coordinates")
        plt.ylabel("Directionality Index (DI) values")
        plt.grid(plot_grid)
        if plot_legend == True:
            plt.legend(prop={'size': 8})
        plt.savefig(os.path.join(output_path, "HiCtool_chr" + a_chr + "_DI_HMM.pdf"), format = 'pdf')
        print("Done!")

def save_topological_domains(a_matrix, output_file):
    """
    Function to save the topological domains coordinates to text file.
    Each topological domain coordinates (start and end) occupy one row and are
    tab separated.
    Arguments:
        a_matrix (numpy matrix): file to be saved with topological domains coordinates.
        output_file (str): output file name in txt format.
    Output:
        Tab separated txt file with topological domain start and end coordinates.
    """
    def compile_row_string(a_row):
        return str(a_row).strip(']').strip('[').lstrip().replace(' ','\t')
    with open(output_file, 'w') as f:
        for row in a_matrix:
            f.write(compile_row_string(row)+'\n')


def load_topological_domains(input_file):
    """
    Function to load the topological domains coordinates from txt file.
    Arguments:
        input_file (str): input file name generated with "calculate_topological_domains" in txt format.
    Return:
        List of lists with topological domain coordinates.
    """
    import csv
    print("Loading topological domain coordinates...")
    with open(input_file, 'r') as f:
        reader = csv.reader(f, dialect='excel', delimiter='\t')
        topological_domains = []
        for row in reader:
            row_int = [int(x) for x in row]
            topological_domains.append(row_int)
        print("Done!")
        return topological_domains
    

def calculate_chromosome_topological_domains(input_file_hmm,
                                             a_chr,
                                             bin_size=40000,
                                             output_path="tad_analysis"):
    """
    Function to calculate the topological domains coordinates of a chromosome. It takes the
    HMM states as input. Topological domains are stored in each line with tab separated start and end coordinates.
    Arguments:
        input_file_hmm (str | obj): txt file of the HMM states generated with the function 
        "calculate_chromosome_hmm_states" or object with the true DI values returned by "calculate_chromosome_hmm_states".
        a_chr (str): chromosome number (example for chromosome 1: '1').
    Returns: List of lists with topological domain coordinates.
    Output: Tab separated txt file with the topological domain coordinates.
    """
    import numpy as np    
    
    print("Calculating topological domain coordinates...")
    if isinstance(input_file_hmm,str):
        likelystates = load_hmm_states(input_file_hmm)
    else:
        likelystates = input_file_hmm
       
    # Start coordinates of the domains
    p = []
    for i in range(1,len(likelystates)):
        if (likelystates[i] == 1 and likelystates[i-1] == 2) or (likelystates[i] == 1 and likelystates[i-1] == 0):
            p.append(i * bin_size)
    
    # End coordinates of the domains
    n = []
    for i in range(1,len(likelystates)-1):
        if (likelystates[i] == 2 and likelystates[i+1] == 1) or (likelystates[i] == 2 and likelystates[i+1] == 0):
            n.append(i * bin_size)
    
    if len(p) == 0 or len(n) == 0:
        print("WARNING! No topological domains can be detected in chromosome " + a_chr)
        return
    
    p1 = 0
    n1 = 0
    p2 = 1
    n2 = 1
    
    # Step 1: checking if the first negative values are greater than the first positive value.
    while n[n1] < p[p1]:
        n1 = n1 + 1
        n2 = n2 + 1
    
    # Now we have removed all the first negative values before the first positive one.
    topological_domains = []
    while p1 < len(p)-1 and n1 < len(n)-1:
        # Step 2: checking if there are two consecutive positive values.
        while n[n1] > p[p2] and p2 < len(p)-1:
            p2 = p2 + 1
        # Now we have removed the possible gaps between consecutive positive states.
    
        # Step 3: checking if there are two consecutive negative values.
        while n[n2] < p[p2] and n2 < len(n)-1:
            n1 = n1 + 1
            n2 = n2 + 1
        # Now we have removed the possible gaps between consecutive negative states.
    
        # Step 4: identification of the Topological Domain.
        topological_domains.append([p[p1],n[n1]])
        p1 = p2
        n1 = n2
        p2 = p1 + 1
        n2 = n1 + 1
    
    save_topological_domains(np.matrix(topological_domains), os.path.join(output_path, "HiCtool_chr" + a_chr + "_topological_domains.txt"))
    print("Done!")
    return topological_domains


def full_tad_analysis(input_contact_matrix,
                      a_chr,
                      isGlobal,
                      tab_sep,
                      species='hg38',
                      data_type='normalized',
                      save_di=True,
                      save_hmm=True,
                      bin_size = 40000,
                      zero_threshold=0.4,
                      output_path="tad_analysis"):
    """
    Compute DI values, HMM states and topological domain coordinates for a chromosome.
    Arguments:
        input_contact_matrix (str | obj): normalized intra-chromosomal contact matrix at a bin size of 40kb passed as a filename (str)
        or an object. Either a single contact matrix or a global contact matrix can be passed (see following parameters).
        a_chr: chromosome number (example for chromosome 1: '1').
        isGlobal (bool): set True if your input matrix is a global matrix (all-by-all chromosomes).
        tab_sep (bool): set True if your input matrix is in a tab separated format. If the matrix is passed as an
        object, this parameter is not taken into consideration.
        species (str): 'hg38' or 'mm10' or any other species label in string format.
        data_type (str): data type, "observed" or "normalized".
        save_di (bool): if True, save the DI values to txt file.
        save_hmm (bool): if True, save the HMM states to txt file.
    Returns: List of lists with topological domain coordinates.
    Output: 
        Txt file containing topological domain coordinates.
        Txt file with the DI values if "save_di=True".
        Txt file with the HMM states if "save_hmm=True".
        Single chromosome contact matrix in compressed format if the input matrix is a global matrix.
    """        
    import copy
    
    if isGlobal == False:
        if isinstance(input_contact_matrix, str):
            if tab_sep == False:
                contact_matrix = load_matrix(input_contact_matrix)
            else:
                contact_matrix = load_matrix_tab(input_contact_matrix)
        else:
            contact_matrix = copy.deepcopy(input_contact_matrix)
    else:
        contact_matrix = extract_single_map(input_global_matrix=input_contact_matrix,
                                            tab_sep=tab_sep,
                                            chr_row=a_chr,
                                            chr_col=a_chr,
                                            species=species,
                                            bin_size=bin_size,
                                            data_type=data_type,
                                            save_output=False,
                                            save_tab=False)

    # DI VALUES
    DI = calculate_chromosome_DI(input_contact_matrix=contact_matrix,
                                 a_chr=a_chr,
                                 isGlobal=False,
                                 tab_sep=False,
                                 data_type=data_type,
                                 species=species,
                                 save_file=save_di,
                                 bin_size=bin_size,
                                 output_path=output_path)
    
    # HMM STATES
    HMM = calculate_chromosome_hmm_states(input_file_DI=DI,
                                          a_chr=a_chr,
                                          save_file=save_hmm,
                                          zero_threshold=zero_threshold,
                                          output_path=output_path)
    
    # TOPOLOGICAL DOMAIN COORDINATES
    tad = calculate_chromosome_topological_domains(input_file_hmm=HMM,
                                                   a_chr=a_chr,
                                                   bin_size=bin_size,
                                                   output_path=output_path)
    return tad



def hictools_call_tads(action='full_tad_analysis', contact_matrix_or_DI=None, chromSizes_path=None, isGlobal=False, tab_sep=True, chrs=None, species=None, data_type='normalized', full_chromosome=True, coord=None, input_file_hmm=None, plot_legend=True, plot_grid=True, bin_size=40000, zero_threshold=0.4, output_path="tad_analysis"):
    
    if species + ".chrom.sizes" not in os.listdir(chromSizes_path):
        available_species = ', '.join([x.split('.')[0] for x in  os.listdir(chromSizes_path)])
        raise Exception('Wrong species inserted! Check the species spelling or insert an available species: ' + available_species + '. If your species is not listed, please contact Riccardo Calandrelli at <rcalandrelli@eng.ucsd.edu>.')
    
    if not path.exists(output_path):
        os.mkdir(output_path)
    
    if action == 'full_tad_analysis':
        if isGlobal == None:
            raise Exception('-h for help or insert 1 if the contact matrix is global (all-by-all chromosomes), 0 otherwise!')
        else:
            pass
        if tab_sep == None:
            raise Exception('-h for help or insert 1 if the contact matrix is in tab separated format, 0 otherwise!')
        else:
            pass
        if data_type == None:
            raise Exception('-h for help or insert a custom label for the data type (observed, normalized, etc.)!')
        else:
            pass
        
        if isGlobal == False:
            if len(chrs) > 1:
                raise Exception('To perform the analysis on multiple chromosomes you must insert a global all-by-all chromosomes matrix.')
            else:
                pass
        
        for c in chrs:
            print("Performing TAD analysis on chr" + c + " ...")
            full_tad_analysis(contact_matrix_or_DI,
                              c,
                              isGlobal,
                              tab_sep,
                              species=species,
                              data_type=data_type,
                              save_di=True,
                              save_hmm=True,
                              bin_size=bin_size,
                              zero_threshold=zero_threshold,
                              output_path=output_path)
            print("Done!")
    

    elif action == 'plot_chromosome_DI':
        
        if full_chromosome == None:
            raise Exception('-h for help or insert 1 if you wish to plot the DI for the entire chromosome, 0 otherwise!')
        else:
            pass
        
        if  len(chrs) > 1:
            raise Exception("Only a single chromosome is accepted if action is plot_chromosome_DI!")
        
        if full_chromosome == False:
            start_pos = coord[0]
            end_pos = coord[1]
        else:
            start_pos = 0
            end_pos = 0
        
        plot_chromosome_DI(contact_matrix_or_DI, 
                           chrs[0],
                           full_chromosome,
                           start_pos, 
                           end_pos,
                           input_file_hmm,
                           species,
                           plot_legend,
                           plot_grid,
                           bin_size=bin_size,
                           output_path="tad_analysis")
