import cooler, numpy, pandas, bioframe

loop_file = 'results/express_and_loop_strength/WT_loop_ZZ.bedpe'
loops = bioframe.read_table(loop_file, schema=['chrom1','start1','end1','chrom2','start2','end2','score'])

dige = pandas.read_csv("MboI.bed", sep='\t')
len_dige = len(dige)
loops['mid1s'] = (loops['start1'] + loops['end1'])//2
loops['mid1e'] = loops['mid1s'] +1
loops['mid2s'] = (loops['start2'] + loops['end2'])//2
loops['mid2e'] = loops['mid2s'] +1
loops['bin1_id'] =  bioframe.closest(loops, dige, return_input=False, return_distance=False, return_index=True, cols1=['chrom1','mid1s','mid1e'])['index_'].astype('int64')
loops['bin2_id'] =  bioframe.closest(loops, dige, return_input=False, return_distance=False, return_index=True, cols1=['chrom2','mid2s','mid2e'])['index_'].astype('int64')

loops_back = loops.copy()

file = 'hic_matrices_MboI/WT-rep1.MboI.cool'
block_size = 100000
nbdige = 12

clr = cooler.Cooler(file)
loops = loops.copy()
loops['bin1_id_block'] = loops['bin1_id'] // block_size
loops['bin2_id_block'] = loops['bin2_id'] // block_size

enzyme_counts_all, a1covs_all, a2covs_all = numpy.zeros(len(loops)), numpy.zeros([len(loops),2*nbdige+1]), numpy.zeros([len(loops),2*nbdige+1])
for index, group in loops.groupby(['bin1_id_block','bin2_id_block']):
    bin1_id_block, bin2_id_block = index
    block_l1 = max(bin1_id_block*block_size-nbdige, 0)
    block_u1 = min((bin1_id_block+1)*block_size+nbdige, len_dige)
    block_l2 = max(bin2_id_block*block_size-nbdige, 0)
    block_u2 = min((bin2_id_block+1)*block_size+nbdige, len_dige)
    mat_block = clr.matrix(balance=False)[block_l1:block_u1,block_l2:block_u2]
    enzyme_counts, a1covs, a2covs = [], [], []
    for bin1_id, bin2_id in zip(group['bin1_id'], group['bin2_id']):
        l1 = max(bin1_id-block_l1-nbdige,0)
        u1 = min(bin1_id-block_l1+nbdige+1,block_u1-block_l1)
        l2 = max(bin2_id-block_l2-nbdige,0)
        u2 = min(bin2_id-block_l2+nbdige+1,block_u2-block_l2)
        submat = mat_block[l1:u1,l2:u2]
        enzyme_counts.append(numpy.sum(submat))
        a1covs.append(numpy.concatenate([[0]*(l1-bin1_id+block_l1+nbdige), numpy.sum(submat, axis=1), [0]*(bin1_id-block_l1+nbdige+1-u1)]))
        a2covs.append(numpy.concatenate([[0]*(l2-bin2_id+block_l2+nbdige), numpy.sum(submat, axis=0), [0]*(bin2_id-block_l2+nbdige+1-u2)]))
    enzyme_counts_all[group.index] = enzyme_counts
    a1covs_all[group.index] = a1covs
    a2covs_all[group.index] = a2covs
