import cooler, ljwlib.hic_module, matplotlib, seaborn, numpy

file = 'hic_matrices_MboI/WT-rep1.MboI.cool'
clrMboI = cooler.Cooler(file)
clr = ljwlib.hic_module.load_rename_add_normVec_cov_tot('hic_matrices/WT-rep1.allValidPairs.mcool', 25_000)
regionx = ('chr1',0, 10_000)
regiony = ('chr1',18_000_000, 28_550_000)
binsx = clrMboI.bins().fetch(regionx)
binsy = clrMboI.bins().fetch(regiony)
pixels = clrMboI.pixels().fetch(regionx)
pixels = pixels[(pixels['bin2_id']>=binsy.index[0]) & (pixels['bin2_id']<=binsy.index[-1])]


pixels['mid1'] = (binsx['start'][pixels['bin1_id']].values + binsx['end'][pixels['bin1_id']].values)//2
pixels['mid2'] = (binsy['start'][pixels['bin2_id']].values + binsy['end'][pixels['bin2_id']].values)//2
binxedges = numpy.append(binsx['start'].values, binsx['end'].values[-1])
binyedges = numpy.append(binsy['start'].values, binsy['end'].values[-1])
f, ax = matplotlib.pyplot.subplots()
seaborn.histplot(data=pixels, x='mid1', y='mid2', weights='count', bins=[list(binxedges),list(binyedges)], binrange=[binxedges[[0,-1]],binxedges[[0,-1]]], ax=ax)
f.savefig("test.pdf")
matplotlib.pyplot.close(f)