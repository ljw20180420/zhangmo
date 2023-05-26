import os, ljwlib.hic_module, numpy, subprocess, shutil, pandas, bioframe

### renbin call tads matlab
def renbin_call_tads_matlab(file, binsize=10000, DIwindows=2000000, boundsizemin=3, boundthres=0.99):
    clr = ljwlib.hic_module.load_rename_add_normVec_cov_tot(file, binsize)
    os.makedirs('matlab_tads/tmp', exist_ok=True)
    for chr in clr.chromnames:
        numpy.savetxt('matlab_tads/tmp/tmp.matrix', numpy.nan_to_num(clr.matrix().fetch(chr), nan=0.), delimiter='\t')
        clr.bins().fetch(chr)[['chrom', 'start', 'end']].to_csv(path_or_buf='matlab_tads/tmp/tmp.bin', sep='\t', header=False, index=False)
        subprocess.check_output('paste matlab_tads/tmp/tmp.bin matlab_tads/tmp/tmp.matrix > matlab_tads/tmp/tmp.di.matrix', shell=True)
        subprocess.check_output(f'perl zmlib/DI_from_matrix.pl matlab_tads/tmp/tmp.di.matrix {clr.binsize} {DIwindows} zmlib/hg19.chrom.sizes > matlab_tads/{os.path.basename(file)}.{chr}.DI', shell=True)
    shutil.rmtree('matlab_tads/tmp')

    genomeDI = pandas.concat([bioframe.read_table(f'matlab_tads/{os.path.basename(file)}.{chr}.DI', schema=['chrom','start','end','DI']) for chr in clr.chromnames]).reset_index(drop=True)
    genomeDI['chrom'][genomeDI['chrom']=='M'] = 23
    genomeDI['chrom'][genomeDI['chrom']=='X'] = 24
    genomeDI['chrom'][genomeDI['chrom']=='Y'] = 25
    genomeDI.sort_values(by=['chrom','start']).to_csv(path_or_buf=f'matlab_tads/{os.path.basename(file)}.DI', sep='\t', na_rep='nan', header=False, index=False)
    subprocess.check_output(f'matlab -nodisplay -r "addpath(\'zmlib\');TADs(\'matlab_tads/{os.path.basename(file)}.DI\', \'matlab_tads/{os.path.basename(file)}.hmm\', \'zmlib/required_modules\');exit"', shell=True)
    subprocess.check_output(f'perl zmlib/file_ends_cleaner.pl matlab_tads/{os.path.basename(file)}.hmm matlab_tads/{os.path.basename(file)}.DI | perl zmlib/converter_7col.pl > matlab_tads/{os.path.basename(file)}.7colfile', shell=True)
    genome7col = bioframe.read_table(f'matlab_tads/{os.path.basename(file)}.7colfile', schema=['chrom','start','end','DI','aic','ind','states'])
    for i in range(1,26):
        if i < 23:
            chr = f'chr{i}'
        else:
            chr = 'chrM' if i==23 else 'chrX' if i==24 else 'chrY'
        genome7col[genome7col['chrom']==f'chr{i}'].to_csv(path_or_buf=f'matlab_tads/{os.path.basename(file)}.{chr}.7colfile', sep='\t', header=False, index=False)
        subprocess.check_output(f'perl zmlib/hmm_probablity_correcter.pl matlab_tads/{os.path.basename(file)}.{chr}.7colfile {boundsizemin} {boundthres} {clr.binsize} | perl zmlib/hmm-state_caller.pl faifile {chr} | perl zmlib/hmm-state_domains.pl > matlab_tads/{os.path.basename(file)}.{chr}.dm', shell=True)
    genomeTAD = pandas.concat([bioframe.read_table(f'matlab_tads/{os.path.basename(file)}.{chr}.dm', schema=['chrom','start','end']) for chr in clr.chromnames]).reset_index(drop=True)
    genomeTAD = genomeTAD[~genomeTAD['end'].isna()].reset_index(drop=True)
    genomeTAD[['start','end']] = genomeTAD[['start','end']].astype("int64")
    genomeTAD['chrom'][genomeTAD['chrom']=='chr23'] = 'chrM'
    genomeTAD['chrom'][genomeTAD['chrom']=='chr24'] = 'chrX'
    genomeTAD['chrom'][genomeTAD['chrom']=='chr25'] = 'chrY'
    genomeTAD = bioframe.sort_bedframe(genomeTAD, view_df = bioframe.make_viewframe(ljwlib.hic_module.chromsizes), df_view_col='chrom').drop_duplicates().to_csv(path_or_buf=f'matlab_tads/{os.path.basename(file)}.dm', sep='\t', header=False, index=False)
    for chr in clr.chromnames:
        os.remove(f'matlab_tads/{os.path.basename(file)}.{chr}.DI')
        os.remove(f'matlab_tads/{os.path.basename(file)}.{chr}.7colfile')
        os.remove(f'matlab_tads/{os.path.basename(file)}.{chr}.dm')

### compare tads
def compare_tads(tadfile1, tadfile2):
    tads1 = bioframe.read_table(tadfile1, schema=["chrom", "start", "end"])
    tads2 = bioframe.read_table(tadfile2, schema=["chrom", "start", "end"])
    tads1['sample'] = 1
    tads2['sample'] = 2
    tads_cluster = bioframe.cluster(pandas.concat([tads1, tads2]), min_dist=None)
    nums1, nums2, pairs = [], [], []
    for _, group in tads_cluster.groupby("cluster"):
        num1 = sum(group['sample']==1)
        num2 = sum(group['sample']==2)
        nums1.append(num1)
        nums2.append(num2)
        if num1==1 and num2==1:
            pairs.append(group.sort_values(by='sample'))

    return pandas.DataFrame({'num1' : nums1, 'num2' : nums2}), pandas.concat(pairs) 