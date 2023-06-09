#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sun May 22 14:37:37 2016

@author: wxt

"""
import argparse, sys, logging, logging.handlers, numpy, domaincaller

currentVersion = domaincaller.__version__


def getargs(commands=None):
    ## Construct an ArgumentParser object for command-line arguments
    parser = argparse.ArgumentParser(usage = '%(prog)s <--uri cool -O output> [options]',
                                     description = '''A python implementation of original DI-based
                                     domain caller proposed by Dixon et al. (2012)''',
                                     formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('-v', '--version', action = 'version',
                        version = ' '.join(['%(prog)s', currentVersion]),
                        help = 'Print version number and exit')

    # Output
    parser.add_argument('--uri',
                        help = 'Cool URI.')
    parser.add_argument('-O', '--output')
    parser.add_argument('-D', '--DI-output')
    parser.add_argument('--window-size', default=2000000, type=int,
                        help='''Window size used for directionality index (DI) calculation.''')
    parser.add_argument('-W', '--weight-col', default='weight',
                        help = '''Name of the column in .cool to be used to construct the
                        normalized matrix. Specify "-W RAW" if you want to run with the raw matrix.''')
    parser.add_argument('--exclude', nargs = '*', default = ['chrY','chrM'],
                        help = '''List of chromosomes to exclude.''')
    parser.add_argument('-p', '--cpu-core', type = int, default = 1,
                        help = 'Number of processes to launch.')
    parser.add_argument('--probs', type = float, default = 0.99,
                        help = 'median probability threshold for consecutive uniform state.')
    parser.add_argument('--minsize', type = int, default = 3,
                        help = 'minimal endurable consecutive bin number with the same state.')
    
    parser.add_argument('--logFile', default = 'domaincaller.log', help = '''Logging file name.''')
    
    ## Parse the command-line arguments
    if not commands:
        commands = sys.argv[1:]
    if not commands:
        commands.append('-h')
    args = parser.parse_args(commands)
    
    return args, commands



def call_tads(args):
    ## Root Logger Configuration
    logger = logging.getLogger()
    logger.setLevel(10)
    console = logging.StreamHandler()
    filehandler = logging.handlers.RotatingFileHandler(args.logFile,
                                                        maxBytes = 200000,
                                                        backupCount = 5)
    # Set level for Handlers
    console.setLevel('INFO')
    filehandler.setLevel('DEBUG')
    # Customizing Formatter
    formatter = logging.Formatter(fmt = '%(name)-25s %(levelname)-7s @ %(asctime)s: %(message)s',
                                    datefmt = '%m/%d/%y %H:%M:%S')
    ## Unified Formatter
    console.setFormatter(formatter)
    filehandler.setFormatter(formatter)
    # Add Handlers
    logger.addHandler(console)
    logger.addHandler(filehandler)
    
    ## Logging for argument setting
    arglist = ['# ARGUMENT LIST:',
                '# Output TAD file = {0}'.format(args.output),
                '# Output DI file = {0}'.format(args.DI_output),
                '# Cool URI = {0}'.format(args.uri),
                '# Window Size = {0}'.format(args.window_size),
                '# Column for matrix balancing = {0}'.format(args.weight_col),
                '# Excluded Chromosomes = {0}'.format(args.exclude),
                '# Number of processes used = {0}'.format(args.cpu_core),
                '# Log file name = {0}'.format(args.logFile)
                ]
    
    argtxt = '\n'.join(arglist)
    logger.info('\n' + argtxt)

    try:
        if not args.weight_col in ['weight']:
            correct = False
        else:
            correct = args.weight_col

        G = domaincaller.genomeLev.Genome(args.uri, balance_type=correct, window=args.window_size, exclude=args.exclude)
        logger.info('Fit Hidden Markov Models with up to 20 mixtures ...')
        #maxM = 3
        maxM = 20
        aic = numpy.zeros(maxM)
        candidates = []
        for i in range(1, maxM+1):
            logger.info('Mixture Number: {0}'.format(i))

            K = 0 # Number of Parameters for AIC Calculation
            # Gaussian Mixture Model Parameters
            K += (3 * (3 * i))
            # HMM parameters
            K += 9 # transition matrix
            #K += 3 # start probabilities

            model = G.oriHMMParams(numdists=i)
            model.fit(G.training_data, algorithm='baum-welch', max_iterations=500,
                        stop_threshold=1e-5, n_jobs=args.cpu_core, verbose=False) # max_iterations is decrease from 10000 to 500 by ljw because in matlab version, it is 500
            candidates.append(model)
            #print("Edges: {0}".format(model.get_params()['edges']))
            
            logL = sum(model.log_probability(seq) for seq in G.training_data)
            aic[i-1] = -2 * logL + 2 * K

            logger.info('Log likelihood: {0}, AIC value: {1}, Parameter Number: {2}'.format(logL, aic[i-1], K))
        
        # Apply AIC
        order = numpy.int(numpy.floor(numpy.log10(abs(aic.min())))) - 1
        div = numpy.power(10, order)
        
        # Relative probability
        for i in range(maxM):
            p_aic = numpy.exp((aic.min() - aic[i]) / (div * 2))
            if p_aic >= 0.9:
                idx = i
                break
        model = candidates[idx]
        logger.info('HMM with {0} mixtures achieved the best performance'.format(idx+1))

        logger.info('Inferring TADs ...')
        DIout = open(args.DI_output, 'w')
        out = open(args.output, 'w')
        for c in G.chroms:
            logger.info('{0} ...'.format(c))
            tdata = domaincaller.chromLev.extract_matrix(G.hic, c, correct)
            work = domaincaller.chromLev.Chrom(c, G.hic.binsize, tdata)
            work.callDomains(model, window=args.window_size, probs=args.probs, minsize=args.minsize)
            for d in work.domains:
                out.write('{0}\t{1}\t{2}\n'.format(c, d[0], d[1]))
            for i, v in enumerate(work.DIs):
                start = i * G.hic.binsize
                end = min(start + G.hic.binsize, G.hic.chromsizes[c])
                DIout.write('{0}\t{1}\t{2}\t{3:.4g}\n'.format(c, start, end, v))

        DIout.close()
        out.close()
        logger.info('Done!')
    except:
        raise
        # traceback.print_exc(file = open(args.logFile, 'a'))
        # sys.exit(1)


def run():
    # Parse Arguments
    args, commands = getargs()
    # Improve the performance if you don't want to run it
    if commands[0] not in ['-h', '-v', '--help', '--version']:
        call_tads(args)
        
if __name__ == '__main__':
    run()