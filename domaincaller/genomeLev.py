import cooler, logging
from domaincaller.chromLev import Chrom, extract_matrix
import numpy as np
from pomegranate.gmm import GeneralMixtureModel
from pomegranate.distributions import Normal
from pomegranate.hmm import DenseHMM

log = logging.getLogger(__name__)

### define the induced class of NormalDistribution which has a nonzero minimal standard deviation

class Genome(object):
    
    def __init__(self, uri, balance_type='weight', window=2000000, exclude=[]):
        
        self.hic = cooler.Cooler(uri)
        res = self.hic.binsize
        self.training_data = []
        log.debug('Calculating DIs for each chromosome ...')
        self.chroms = []
        for c in self.hic.chromnames:
            if c in exclude:
                continue
            log.debug('Chrom {0} ...'.format(c))
            self.chroms.append(c)
            tdata = extract_matrix(self.hic, c, balance_type)
            work = Chrom(c, res, tdata)
            work.calDI(window=window)
            self.training_data.append(work.DIs)

    def oriHMMParams(self, numdists=3):
        """
        Set initial parameters for the Hidden Markov Model (HMM).
        
        """
        # GMM emissions
        # 3 Hidden States:
        # 0--downstream (right boundary), 1--no bias, 2--upstream (left boundary)        

        model = GeneralMixtureModel([Normal(min_cov=1e-6)]*3*numdists, max_iter=500)
        model.fit(np.concatenate(self.training_data).reshape(-1,1).astype(np.float32))

        breakpoint()

        mus = np.array([distri.means[0] for distri in model.distributions])
        idx = np.argsort(mus)
        weight = model.weights[idx].reshape(3, numdists)
        starts = np.sum(weight, 1)
        weight = weight / starts.reshape(3,1)
        starts = starts / np.sum(starts)
        idx = idx.reshape(3, numdists)
        dists = []
        for i in range(3):
            if (numdists==1):
                dists.append(model.distributions[idx[i]][0])
            else:
                dists.append(GeneralMixtureModel(model.distributions[idx[i,:]], weight[i,:]))

        # transition matrix
        A = [[0.34, 0.33, 0.33],
            [0.33, 0.34, 0.33],
            [0.33, 0.33, 0.34]]

        hmm = DenseHMM.from_matrix(A, dists, starts, state_names=['0', '1', '2'], name='mixture{0}'.format(numdists))
        
        return hmm
        

