import requests, gzip, os, bioframe, numpy, scipy

### download something
file='https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/genes/hg19.ncbiRefSeq.gtf.gz'
r = requests.get(file)
with open(os.path.basename(file), 'wb') as f:
    f.write(r.content)

# decompress version
r = requests.get('https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/genes/hg19.ensGene.gtf.gz')
with open("./hg19.ensGene.gtf", 'wb') as f:
    f.write(gzip.decompress(r.content))

### digest genome
genome = bioframe.load_fasta('/home/ljw/hg19_with_random/hg19.fa')
emzyme = 'MboI'
bins = bioframe.digest(genome, emzyme)
bins.to_csv(path_or_buf=f'{emzyme}.bed', sep='\t', header=True, index=False)

### nested repeated measure
def nested_repeated_measure(df, colgroups):
    arep = scipy.sparse.csc_array(numpy.concatenate([numpy.eye(df.shape[0])]*df.shape[1], axis=0))
    agru = []
    for colgroup in colgroups:
        agru.append(numpy.zeros(df.shape))
        idx = [list(df.columns.values).index(col) for col in colgroup]
        agru[-1][:,idx] = 1
        agru[-1] = numpy.reshape(agru[-1],[-1,1])
    agru = scipy.sparse.csc_array(numpy.concatenate(agru, axis=1))
    respo = numpy.reshape(df.values,[-1,1])
    r0 = scipy.sparse.linalg.lsqr(arep, respo)[3]
    aall = scipy.sparse.hstack([arep,agru])
    r1 = scipy.sparse.linalg.lsqr(aall, respo)[3]
    rank0 = arep.shape[1]
    rank1 = aall.shape[1]-1
    dfn = rank1-rank0
    dfd = len(respo)-rank1
    return 1-scipy.stats.f.cdf((r0-r1)*dfd/r1/dfn, dfn, dfd)

### chi-square two multinomial distributions are the same
def chi_square_two_multinomial_same(arr1, arr2):
    s1, s2 = numpy.sum(arr1), numpy.sum(arr2)
    d0 = (arr1 + arr2) / (s1 + s2)
    E1, E2= s1 * d0, s2 * d0
    chi2 = numpy.sum((arr1 - E1)**2 / E1) + numpy.sum((arr2 - E2)**2 / E2)
    # return 1 - scipy.stats.chi2.logsf(chi2, len(arr1)-1)
    return chi2