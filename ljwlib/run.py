import functools, pickle, joblib, sys, pygenometracks.makeTracksFile, pygenometracks.plotTracks, seaborn, shutil, sklearn.ensemble
from ljwlib.hic_module import *
import ljwlib.softwares
import ljwlib.tads
import domaincaller.renbin
import HiCtool.scripts.HiCtool_TAD_analysis_lib
import TADtree.TADtreelib
import ljwlib.scrab_ucsc
import ljwlib.sundry

### to print which resolutions are stored in the mcool, use list_coolers
# cooler.fileops.list_coolers(file)
### initialization


filesWTmerge = [f"hic_matrices/{file}" for file in ['WT.allValidPairs.mcool', 'ZO1.allValidPairs.mcool', 'ZO2.allValidPairs.mcool']]
filesWTrep = [f"hic_matrices/{file}" for file in ['WT-rep1.allValidPairs.mcool', 'WT-rep2.allValidPairs.mcool', 'ZO1-rep1.allValidPairs.mcool', 'ZO1-rep2.allValidPairs.mcool', 'ZO2-rep1.allValidPairs.mcool', 'ZO2-rep2.allValidPairs.mcool']]
filesWT = filesWTmerge + filesWTrep
filesAIDmerge = [f"hic_matrices/{file}" for file in ['AID.allValidPairs.mcool', 'AID-ZO.allValidPairs.mcool', 'AID-CD.allValidPairs.mcool', 'AID-ZOCD.allValidPairs.mcool']]
filesAIDrep = [f"hic_matrices/{file}" for file in ['AID-rep1.allValidPairs.mcool', 'AID-rep2.allValidPairs.mcool', 'AID-ZO-rep1.allValidPairs.mcool', 'AID-ZO-rep2.allValidPairs.mcool', 'AID-CD-rep1.allValidPairs.mcool', 'AID-CD-rep2.allValidPairs.mcool', 'AID-ZOCD-rep1.allValidPairs.mcool', 'AID-ZOCD-rep2.allValidPairs.mcool']]
filesAID = filesAIDmerge + filesAIDrep
files = filesWT + filesAID


### print full genome
for file in files:
    clr = load_rename_add_normVec_cov_tot(file, 1_000_000, force=True)
    print_hic_map(clr, 'results/display', ylim=[0,1], cmap='bwr', norm=matplotlib.colors.LogNorm(vmin=1e-1,vmax=1e5))

    # print all chromosomes
    clr = load_rename_add_normVec_cov_tot(file, 100_000)
    for chr in clr.chromnames:
        print_hic_map(clr, 'results/display', region=(chr, 0, clr.chromsizes[chr]))

### print compartments and saddles for all chromosomes
for file in files:
    clr = load_rename_add_normVec_cov_tot(file, 100_000)
    for chr in clr.chromnames:
        compartmentalization(clr, genome, chr, 'results/compartment', cmap='bwr', norm=matplotlib.colors.Normalize(vmin=-0.5,vmax=0.5))
        saddle_plots(clr, genome, chr, 'results/saddle', vmin=-2, vmax=2, norm=matplotlib.colors.LogNorm(vmin=1e-1,vmax=1e1), cmap='bwr')

### slide correlation of compartments
# loop_file = 'zmloops/AID_loop_ZZ.bedpe'
file = "hic_matrices/AID.allValidPairs.hic"
binsize = 25000
clr = load_rename_add_normVec_cov_tot(f"{file[-3]}mcool", binsize)
file = file[-5]
chr = "chr1"
subprocess.check_output(f"java -jar juicer_tools_1.22.01.jar dump oe SCALE {file} {chr} {chr} BP {binsize} results/slide_correlation/{os.path.basename(file)}.{chr}", shell=True)
spamat = bioframe.read_table(f"results/slide_correlation/{os.path.basename(file)}.{chr}", schema=['start', 'end', 'count'])
spamat['bin1id'] = spamat['start']//binsize + 1
spamat['bin2id'] = spamat['end']//binsize + 1
spamat.reindex(columns=['bin1id', 'bin2id', 'count']).to_csv(f"results/slide_correlation/{os.path.basename(file)}.{chr}.sparse", index=False, header=False, sep='\t')
bins = clr.bins().fetch(chr)[['chrom', 'start', 'end']]
bins['binid'] = bins.index + 1
bins.to_csv(f'results/slide_correlation/{os.path.basename(file)}.{chr}.bed', index=False, header=False, sep='\t')
subprocess.check_output(f"python zmlib/sparseToDense.py -b results/slide_correlation/{os.path.basename(file)}.{chr}.bed results/slide_correlation/{os.path.basename(file)}.{chr}.sparse -o results/slide_correlation/{os.path.basename(file)}.{chr}.dense", shell=True)
denmat = numpy.loadtxt(f"results/slide_correlation/{os.path.basename(file)}.{chr}.dense")
for vmax in range(1,10,1):
    f, ax, _ = binheat(denmat, bins, norm=matplotlib.colors.Normalize(vmin=0,vmax=vmax), cmap="Oranges")
    f.savefig(f"results/slide_correlation/{os.path.basename(file)}.{chr}.{vmax}.oe.pdf")


### print P(s) for chromosomes
for file in files:
    clr = load_rename_add_normVec_cov_tot(file, 10_000)
    for chr in clr.chromnames:
        draw_Ps(clr, hg19_arms, chr, 'results/Ps')

### print average trans-frequencies
for file in files:
    clr = load_rename_add_normVec_cov_tot(file, 10_000)
    average_trans_frequency(clr, 'results/trans_avg', cmap='bwr', norm=matplotlib.colors.LogNorm(vmin=1e-4,vmax=1e3))

### tads related analyses of insulation score
clr = load_rename_add_normVec_cov_tot(file, 10_000)
window = 3*clr.binsize
insulation_table, threshold_li, threshold_otsu = tads_strength(clr, window,'results/tads_insul')
draw_region_tads(clr, ('chr15', 140100000, 140900000), window, insulation_table, 'results/tads_insul', norm=matplotlib.colors.Normalize(vmin=0,vmax=4), cmap='Oranges')
boundary_signal(clr, window, 'ChIP_Nexus_rawdata/WT_CTCF_ChIP-nexus.adp.bam.bw', insulation_table, threshold_li, threshold_otsu, 'results/tads_insul')
boundary_signal(clr, window, 'ChIP_Nexus_rawdata/WT_ZNF143_ChIP-nexus.adp.bam.bw', insulation_table, threshold_li, threshold_otsu, 'results/tads_insul')

### call loops
for file in files:
    clr = load_rename_add_normVec_cov_tot(file, 10_000)
    # draw_kernel(clr, 'results/loop')
    # _, _ = loop_calling(clr, ('chr15', 140100000, 140900000), 'results/loop')
    ## pile up
    cvd = cooltools.expected_cis(clr=clr, view_df=hg19_arms, ignore_diags=0, nproc=12)
    dots_df = cooltools.dots(clr, expected=cvd, view_df=hg19_arms, max_loci_separation=10_000_000, n_lambda_bins=70, nproc=12)
    _ = pile_up_clips(clr, dots_df, cvd, hg19_arms, 'dots', 'results/loop')
    ## pile up at CTCF peaks
    # peaks = bioframe.read_table('ChIP_Nexus_results/WT_CTCF_ChIP-nexus_peaks.narrowPeak', schema='bed').query(f'chrom in {clr.chromnames}')
    # pile_up_clips(clr, peaks, cvd, hg19_arms, 'peaks', 'result/loop')
    # paired_peaks = bioframe.pair_by_distance(peaks, min_sep=200_000, max_sep=1_000_000, suffixes=('1', '2'))
    # pile_up_clips(clr, paired_peaks, cvd, hg19_arms, 'paired_peaks', 'result/loop')

### pairwise compartment compare
## get summary of compartments in EigFrame
binsize = 100_000
clr = load_rename_add_normVec_cov_tot(files[0], binsize)
gc_cov = bioframe.frac_gc(clr.bins()[:][['chrom', 'start', 'end']], genome)
chrs, types, E1s, samples, chrindices = [], [], [], [], []
for file in files: 
    for chr in chromindices.keys():
        try:
            E1 = numpy.loadtxt(f'results/compartment/{os.path.basename(file)}.{chr}.{binsize}.eigen')
            GC = bioframe.select(gc_cov, (chr, 0, clr.chromsizes[chr])).GC
            try:
                if scipy.stats.pearsonr(numpy.nan_to_num(GC,nan=0.0), numpy.nan_to_num(E1,nan=0.0)).statistic<0:
                    E1 = -E1 
            except:
                raise     
        except Exception as err:
            print(repr(err))
            continue
        E1s.extend(E1)
        types.extend(['A' if va>0 else 'B' if va<0 else 'N' for va in E1])
        samples.extend([os.path.basename(file)] * len(E1))
        chrs.extend([chr] * len(E1))
        chrindices.extend([chromindices[chr]] * len(E1))
EigFrame = pandas.DataFrame({'sample' : samples, 'chr' : chrs, 'type' : types, 'E1' : E1s, 'chrindex' : chrindices})
EigFrame.to_csv(path_or_buf='results/compartment/EigFrame.tsv', sep='\t', header=True, index=False)

EigFrame = pandas.read_csv('results/compartment/EigFrame.tsv', sep='\t')
filesWTmergerep = filesWTmerge[:1] + filesWTrep[:2] + filesWTmerge[1:2] + filesWTrep[2:4] + filesWTmerge[2:] + filesWTrep[4:]
filesAIDmergerep = filesAIDmerge[:1] + filesAIDrep[:2] + filesAIDmerge[2:3] + filesAIDrep[4:6] + filesAIDmerge[1:2] + filesAIDrep[2:4] + filesAIDmerge[3:] + filesAIDrep[6:]
for groupname, group in zip(['WT', 'AID'], [filesWTmergerep, filesAIDmergerep]):
    f, axs = matplotlib.pyplot.subplots(figsize=(40,40), nrows=len(group), ncols=len(group), sharex=True)
    for i in range(len(group)):
        datai = EigFrame[EigFrame['sample']==os.path.basename(group[i])]
        seaborn.histplot(datai, x='chrindex', hue='type', hue_order=['A','B','N'], palette={'A':'r', 'B':'b', 'N':'g'}, weights='E1', discrete=True, ax=axs[i,i])
        axs[i,i].set(xticks=list(chromindices.values()), xticklabels=list(chromindices.keys()), xlim=[-0.5,len(chromindices)-0.5], ylabel=os.path.basename(group[i]))
        axs[i,i].tick_params(axis='x', rotation=90)
        for j in range(i+1,len(group)):
            dataj = EigFrame[EigFrame['sample']==os.path.basename(group[j])]
            pears = []
            for chr in chromindices.keys():
                try:
                    pear = scipy.stats.pearsonr(numpy.nan_to_num(datai[datai.chr==chr].E1,nan=0), numpy.nan_to_num(dataj[dataj.chr==chr].E1,nan=0)).statistic
                    pears.append(abs(pear))
                except:
                    pears.append(numpy.nan)
            for pe in [[i,j],[j,i]]:
                seaborn.barplot(x=list(chromindices.values()), y=pears, width=1, ax=axs[pe[0],pe[1]])
                axs[pe[0],pe[1]].set(xticks=list(chromindices.values()), xticklabels=list(chromindices.keys()))
                axs[pe[0],pe[1]].tick_params(axis='x', rotation=90)
    f.tight_layout()
    f.savefig(f'results/compartcompare/{groupname}.comcom.pdf')
    matplotlib.pyplot.close(f)


### count the non-HK and HK covered by ZNF143 and TAD boundary
ZNF143 = bioframe.read_table('ChIP_Nexus_results/WT_ZNF143_ChIP-nexus_peaks.narrowPeak', schema='bed')
tss = get_tss('hg19.ncbiRefSeq.gtf.gz')
tss = tss.drop_duplicates(subset=['start'], ignore_index=True)
HK_genes = bioframe.read_table('Cell_HK_genes.txt', schema=['gene_id', 'transcript_id'])
for i in range(len(HK_genes.gene_id)):
    HK_genes.gene_id[i] = HK_genes.gene_id[i].strip(' ')
tss['HK'] = ['HK' if HK else None for HK in tss.gene_id.isin(HK_genes.gene_id)]
tss = standard_assign(tss, ZNF143, '', 'ZNF143', select_col='distance', includes=['name'], select_mode='min')
tss = tss[tss.chrom.isin(list(chromindices.keys()))].reset_index(drop=True)
# deeptools-like plots
flank = 100
bins = 30
ZNF143_signal_NONHK = bbi.stackup('ChIP_Nexus_rawdata/WT_ZNF143_ChIP-nexus.adp.bam.bw', tss[tss.HK.isna()].chrom, tss[tss.HK.isna()].start-flank, tss[tss.HK.isna()].start+flank, bins=bins, summary='mean')
ZNF143_signal_HK = bbi.stackup('ChIP_Nexus_rawdata/WT_ZNF143_ChIP-nexus.adp.bam.bw', tss[~tss.HK.isna()].chrom, tss[~tss.HK.isna()].start-flank, tss[tss.HK.isna()].start+flank, bins=bins, summary='mean')
for tssname, ZNF143_signal in zip(['NONHK', 'HK'], [ZNF143_signal_NONHK, ZNF143_signal_HK]):
    deeptools_like(ZNF143_signal, f'results/enrichment/ZNF143.at.{tssname}.pdf', extent=(-flank, flank, ZNF143_signal.shape[0], 0))

# bar plot with chi2 test of independence
observed = pandas.DataFrame({'HK' : ['NONHK', 'HK'], 'NAZNF143' : [sum((tss.HK.isna()) & (tss.ZNF143_name.isna())), sum((~tss.HK.isna()) & (tss.ZNF143_name.isna()))], 'ZNF143' : [sum((tss.HK.isna()) & (~tss.ZNF143_name.isna())), sum((~tss.HK.isna()) & (~tss.ZNF143_name.isna()))]})
f, ax = matplotlib.pyplot.subplots()
observed.set_index('HK').plot(kind='bar', stacked=True, color=['steelblue', 'red'], ax=ax)
pvalue = scipy.stats.chi2_contingency(observed.set_index('HK').values).pvalue
ax.text(0.75, ax.get_ylim()[1], f'pvalue={pvalue}', ha='center', va='top')
f.tight_layout() 
f.savefig('results/enrichment/ZNF143.tss.HK.ctg.pdf')
matplotlib.pyplot.close(f)



### pile_up ZNF143 hichip loops
CTCF = bioframe.read_table('ChIP_Nexus_results/WT_CTCF_ChIP-nexus_peaks.narrowPeak', schema='bed')
ZNF143 = bioframe.read_table('ChIP_Nexus_results/WT_ZNF143_ChIP-nexus_peaks.narrowPeak', schema='bed')
RAD21 = bioframe.read_table('ChIP_Nexus_results/WT_RAD21_ChIP-nexus_peaks.narrowPeak', schema='bed')
tss = get_tss('hg19.ncbiRefSeq.gtf.gz')
tss = tss.drop_duplicates(subset=['start'], ignore_index=True)
HK_genes = bioframe.read_table('Cell_HK_genes.txt', schema=['gene_id', 'transcript_id'])
tss['HK'] = ['HK' if HK else None for HK in tss.gene_id.isin(HK_genes.gene_id)]
ZNF143 = standard_assign(ZNF143, tss, '', 'TSS', select_col='distance', includes=['name', 'gene_id', 'HK'], select_mode='min')



loop_files = ['hichip_results/WT_ZNF143_HiChIP/WT_ZNF143_HiChIP.filt.intra.loop_counts.bedpe', 'WT_RAD21_loop.bedpe', 'WT_CTCF_loop.bedpe', 'ZNF143KO_CTCF_loop.bedpe', 'ZNF143KO_RAD21_loop.bedpe']


# generate attach loops
for loop_file in loop_files:
    loop = bioframe.read_table(loop_file, sep='\t', schema='bedpe')

    loop = standard_assign(loop, ZNF143, 'a1', 'ZNF143', select_col='distance', includes=['name', 'TSS_gene_id', 'TSS_HK'], select_mode='min', cols1=('chrom1', 'start1', 'end1'))
    loop = standard_assign(loop, ZNF143, 'a2', 'ZNF143', select_col='distance', includes=['name', 'TSS_gene_id', 'TSS_HK'], select_mode='min', cols1=('chrom2', 'start2', 'end2'))
    loop = standard_assign(loop, CTCF, 'a1', 'CTCF', select_col='distance', includes=['name'], select_mode='min', cols1=('chrom1', 'start1', 'end1'))
    loop = standard_assign(loop, CTCF, 'a2', 'CTCF', select_col='distance', includes=['name'], select_mode='min', cols1=('chrom2', 'start2', 'end2'))
    loop = standard_assign(loop, tss, 'a1', 'TSS', select_col='distance', includes=['name', 'strand', 'HK'], select_mode='min', cols1=('chrom1', 'start1', 'end1'))
    loop = standard_assign(loop, tss, 'a2', 'TSS', select_col='distance', includes=['name', 'strand', 'HK'], select_mode='min', cols1=('chrom2', 'start2', 'end2'))

    loop.to_csv(path_or_buf=f'results/SBS_only_loop_type/{loop_file}.attach', sep='\t', header=True, index=False)
    

# analyze the relationship between loop and tss
for loop_file in loop_files:
    loop = pandas.read_csv(f'results/SBS_only_loop_type/{loop_file}.attach', sep='\t')
    # customize loop columns
    loop['a1_TSS_has'] = ['F' if na else 'T' for na in loop.a1_TSS_name.isna()]
    loop['a2_TSS_has'] = ['F' if na else 'T' for na in loop.a2_TSS_name.isna()]
    loop['TSS_has'] = 'TSS.'+loop['a1_TSS_has'] + '.TSS.' + loop['a2_TSS_has']
    loop['a1_ZNF143_TSS_has'] = ['F' if na else 'T' for na in loop.a1_ZNF143_TSS_gene_id.isna()]
    loop['a2_ZNF143_TSS_has'] = ['F' if na else 'T' for na in loop.a2_ZNF143_TSS_gene_id.isna()]
    loop['ZNF143_TSS_has'] = 'ZNFTSS.'+loop['a1_ZNF143_TSS_has'] + '.ZNFTSS.' + loop['a2_ZNF143_TSS_has']
    loop['a1_HK'] = ['F' if na else 'T' for na in loop.a1_TSS_HK.isna()]
    loop['a2_HK'] = ['F' if na else 'T' for na in loop.a2_TSS_HK.isna()]
    loop['HK'] = 'HK.' + loop['a1_HK'] + '.HK.' + loop['a2_HK']
    loop['a1_ZNF143_HK'] = ['F' if na else 'T' for na in loop.a1_ZNF143_TSS_HK.isna()]
    loop['a2_ZNF143_HK'] = ['F' if na else 'T' for na in loop.a2_ZNF143_TSS_HK.isna()]
    loop['ZNF143_HK'] = 'HK.' + loop['a1_ZNF143_HK'] + '.HK.' + loop['a2_ZNF143_HK']
    # initialize the filters
    boolESES = (~loop.a1_ZNF143_name.isna()) & (~loop.a2_ZNF143_name.isna())
    boolESNS = (~loop.a1_ZNF143_name.isna()) & (loop.a2_ZNF143_name.isna()) | (loop.a1_ZNF143_name.isna()) & (~loop.a2_ZNF143_name.isna())
    boolNSNS = (loop.a1_ZNF143_name.isna()) & (loop.a2_ZNF143_name.isna())
    boolNCNC = (loop.a1_CTCF_name.isna()) & (loop.a2_CTCF_name.isna())
    boolECNC = (~loop.a1_CTCF_name.isna()) & (loop.a2_CTCF_name.isna()) | (loop.a1_CTCF_name.isna()) & (~loop.a2_CTCF_name.isna())
    boolECEC = (~loop.a1_CTCF_name.isna()) & (~loop.a2_CTCF_name.isna())

    data = loop[boolESES & boolNCNC]

    ho1 = ['TSS.T.TSS.T', 'TSS.T.TSS.F', 'TSS.F.TSS.T', 'TSS.F.TSS.F']
    ho2 = ['TSS.T.TSS.T.HK.T.HK.T', 'TSS.T.TSS.T.HK.T.HK.F', 'TSS.T.TSS.T.HK.F.HK.T', 'TSS.T.TSS.T.HK.F.HK.F', 'TSS.T.TSS.F.HK.T.HK.F', 'TSS.T.TSS.F.HK.F.HK.F', 'TSS.F.TSS.T.HK.F.HK.T', 'TSS.F.TSS.T.HK.F.HK.F', 'TSS.F.TSS.F.HK.F.HK.F']
    ho = []
    zho1 = [ho.replace('TSS', 'ZNFTSS') for ho in ho1]
    zho2 = [ho.replace('TSS', 'ZNFTSS') for ho in ho2]

    f = seaborn.displot(data, x='TSS_has', hue='HK', hue_order=numpy.sort(data['HK'].unique()), col_order=numpy.sort(data['TSS_has'].unique()), multiple='stack')
    f.ax.tick_params(axis='x', rotation=90)
    f.savefig(f'results/SBS_only_loop_type/{os.path.basename(loop_file)}.TSS.count.pdf')

    f = seaborn.jointplot(data, x='a1_TSS_closest', y='a2_TSS_closest', hue=data['TSS_has']+'.'+data['HK'], hue_order=numpy.sort((data['TSS_has']+'.'+data['HK']).unique()), alpha=0.2)
    f.savefig(f'results/SBS_only_loop_type/{os.path.basename(loop_file)}.TSS.distance.pdf')

    f = seaborn.displot(data, x='score', hue=data['TSS_has']+'.'+data['HK'], hue_order=numpy.sort((data['TSS_has']+'.'+data['HK']).unique()), discrete=True, multiple='stack', binrange=[0,10])
    f.savefig(f'results/SBS_only_loop_type/{os.path.basename(loop_file)}.TSS.score.pdf')

    f = seaborn.displot(data, x='ZNF143_TSS_has', hue='ZNF143_HK', hue_order=numpy.sort(data['ZNF143_HK'].unique()), col_order=numpy.sort(data['ZNF143_TSS_has'].unique()), multiple='stack')
    f.ax.tick_params(axis='x', rotation=90)
    f.savefig(f'results/SBS_only_loop_type/{os.path.basename(loop_file)}.ZNFTSS.count.pdf')

    f = seaborn.jointplot(data, x='a1_TSS_closest', y='a2_TSS_closest', hue=data['ZNF143_TSS_has']+'.'+data['ZNF143_HK'], hue_order=numpy.sort((data['ZNF143_TSS_has']+'.'+data['ZNF143_HK']).unique()), alpha=0.2)
    f.savefig(f'results/SBS_only_loop_type/{os.path.basename(loop_file)}.ZNFTSS.distance.pdf')

    f = seaborn.displot(data, x='score', hue=data['ZNF143_TSS_has']+'.'+data['ZNF143_HK'], hue_order=numpy.sort((data['ZNF143_TSS_has']+'.'+data['ZNF143_HK']).unique()), discrete=True, multiple='stack', binrange=[0,10])
    f.savefig(f'results/SBS_only_loop_type/{os.path.basename(loop_file)}.ZNFTSS.score.pdf')

# pileup HiChiper loops in HiC
loops = []
samples = []
CVs = []
types = []
for file in files:
    clr = load_rename_add_normVec_cov_tot(file, 10_000)
    cvd = cooltools.expected_cis(clr=clr, view_df=hg19_arms, ignore_diags=0, nproc=12)
    for loop_file in loop_files:
        loop = pandas.read_csv(f'{loop_file}.attach', sep='\t')

        boolESES = (~loop.a1_ZNF143_name.isna()) & (~loop.a2_ZNF143_name.isna())
        boolESNS = (~loop.a1_ZNF143_name.isna()) & (loop.a2_ZNF143_name.isna()) | (loop.a1_ZNF143_name.isna()) & (~loop.a2_ZNF143_name.isna())
        boolNSNS = (loop.a1_ZNF143_name.isna()) & (loop.a2_ZNF143_name.isna())
        boolNCNC = (loop.a1_CTCF_name.isna()) & (loop.a2_CTCF_name.isna())
        boolECNC = (~loop.a1_CTCF_name.isna()) & (loop.a2_CTCF_name.isna()) | (loop.a1_CTCF_name.isna()) & (~loop.a2_CTCF_name.isna())
        boolECEC = (~loop.a1_CTCF_name.isna()) & (~loop.a2_CTCF_name.isna())

        mp = [['ESES','NCNC'], ['ESES','ECNC'], ['ESES','ECEC'], ['NSNS','ECEC'], ['ESNS','ECEC']]
        for SBS_state, CBS_state in mp:
            boolSBS = eval(f'bool{SBS_state}')
            boolCBS = eval(f'bool{CBS_state}')
            if any(boolSBS & boolCBS):
                stack = pile_up_clips(clr, loop[boolSBS & boolCBS].reset_index(drop=True), cvd, hg19_arms, f'{SBS_state}_{CBS_state}_{os.path.basename(loop_file)}', 'results/pileup_hichip_over_hic', flank=100_000)
                CV = [slice[slice.shape[0]//2, slice.shape[1]//2] for slice in stack]
                loops.extend([os.path.basename(loop_file)] * len(CV))
                samples.extend([os.path.basename(file)] * len(CV))
                types.extend([f'{SBS_state}_{CBS_state}'] * len(CV))
                CVs.extend(CV)

mpms = [f'{p[0]}_{p[1]}' for p in mp]
CVFrame = pandas.DataFrame({'loop' : loops, 'sample' : samples, 'type' : types, 'CV' : CVs})
CVFrame.to_csv(path_or_buf='results/pileup_hichip_over_hic/CVFrame', sep='\t', header=True, index=False)
# CVFrame = pandas.read_csv('results/pileup_hichip_over_hic/CVFrame', sep='\t')
for loop_file in loop_files:
    data = CVFrame[CVFrame.loop==os.path.basename(loop_file)]
    for mpm in mpms:
        data = CVFrame[(CVFrame.loop==os.path.basename(loop_file)) & (CVFrame.type==mpm)]
        if len(data) > 0:
            f, ax = matplotlib.pyplot.subplots(figsize=(20, 10))
            seaborn.violinplot(data=data, x='sample', y='CV', dodge=True, orient='v', ax=ax)
            ylim = ax.get_ylim()
            seaborn.stripplot(data=data, x='sample', y='CV', color='k', dodge=True, orient='v', ax=ax)
            ax.set_ylim(ylim)
            ax.tick_params(axis='x', rotation=90)
            f.tight_layout()
            f.savefig(f'results/pileup_hichip_over_hic/{os.path.basename(loop_file)}.{mpm}.centerValues.pdf')
            matplotlib.pyplot.close(f)
            Pvalue = numpy.tile(numpy.nan, (len(files), len(files)))
            for i in range(len(files)):
                CVi = data.loc[data.sample==os.path.basename(files[i]),'CV']
                for j in range(i+1,len(files)):
                    CVj = data.loc[data.sample==os.path.basename(files[j]),'CV']
                    Pvalue[i,j] = scipy.stats.ttest_ind(CVi, CVj, axis=0, nan_policy='omit', alternative='two-sided').pvalue
            with open(f'results/pileup_hichip_over_hic/{os.path.basename(loop_file)}.{mpm}.centerValues.pvalues', 'w') as fi:
                for file in files:
                    fi.write(f'\t&\t{os.path.basename(file)}')
                fi.write('\t\\\\\n')
                for i in range(len(files)):
                    fi.write(os.path.basename(file))
                    for j in range(len(files)):
                        fi.write(f'\t&\t{Pvalue[i,j]}')
                    fi.write('\t\\\\\n')
        

### call motifs
CTCF_wildtype_peak_files = ['ChIP_Nexus_results/WT_CTCF_ChIP-nexus_peaks.narrowPeak','ChIP_Nexus_results/AID_CTCF_ChIP-nexus_peaks.narrowPeak']
CTCF_pwm_files = ['motifs/CTCF_type1.txt','motifs/CTCF_type2.txt','motifs/CTCF_type3.txt']
for CTCF_wildtype_peak_file in CTCF_wildtype_peak_files:
    CTCF_motif_files = call_motifs(CTCF_wildtype_peak_file, CTCF_pwm_files)

ZNF143_wildtype_peak_files = ['ChIP_Nexus_results/WT_ZNF143_ChIP-nexus_peaks.narrowPeak','ChIP_Nexus_results/AID_ZNF143_ChIP-nexus_peaks.narrowPeak']
ZNF143_pwm_files = ['motifs/ZNF143_type1.txt','motifs/ZNF143_type2.txt','motifs/ZNF143_type3.txt']
for ZNF143_wildtype_peak_file in ZNF143_wildtype_peak_files:
    ZNF143_motif_files = call_motifs(ZNF143_wildtype_peak_file, ZNF143_pwm_files)


### ZNF143_peak SBS motif type
ZNF143_pwm_files = ['motifs/ZNF143_type1.txt','motifs/ZNF143_type2.txt','motifs/ZNF143_type3.txt']
ZNF143_wildtype_peak_files = ['ChIP_Nexus_results/WT_ZNF143_ChIP-nexus_peaks.narrowPeak','ChIP_Nexus_results/AID_ZNF143_ChIP-nexus_peaks.narrowPeak']
ZNF143_peak_filess = [['ChIP_Nexus_results/WT_ZNF143_ChIP-nexus_peaks.narrowPeak'],['ChIP_Nexus_results/AID_ZNF143_ChIP-nexus_peaks.narrowPeak','ChIP_Nexus_results/AID-CD_ZNF143_ChIP-nexus_peaks.narrowPeak']]
ZNF143_HiChIP_loop_filess = [['zmloops/WT_loop.bedpe'], ['zmloops/AID_loop_ZZ.bedpe', 'zmloops/AID-CD_loop_ZZ.bedpe']]
for ZNF143_wildtype_peak_file, ZNF143_peak_files, ZNF143_HiChIP_loop_files in zip(ZNF143_wildtype_peak_files, ZNF143_peak_filess, ZNF143_HiChIP_loop_filess):
    SBS = pandas.concat([bioframe.read_table(f"{ZNF143_pwm_file}.{os.path.basename(ZNF143_wildtype_peak_file)}.motifs", schema=['chrom', 'start', 'end', 'name', 'score', 'strand', 'pvalue', 'qvalue']) for ZNF143_pwm_file in ZNF143_pwm_files])
    for ZNF143_peak_file in ZNF143_peak_files:
        ZNF143 = bioframe.read_table(ZNF143_peak_file, schema='bed')
        ZNF143 = standard_assign(ZNF143, SBS, '', 'SBS', 'pvalue', includes=['name', 'strand'], select_mode='min')
        ZNF143['SBS_type'] = ['0' if not SBS_name else SBS_name[SBS_name.find('type')+4] for SBS_name in ZNF143['SBS_name']]
        f, ax = matplotlib.pyplot.subplots(figsize=(20,20))
        seaborn.histplot(data=ZNF143, x='SBS_type', ax=ax)
        for patch in ax.patches:
            ax.text(patch.get_x() + patch.get_width()/2., patch.get_height(), f'{int(patch.get_height())}', ha="center", va="bottom") 
        f.savefig(f'results/motif_type_count/ZNF143_peak.{os.path.basename(ZNF143_peak_file)}.motif.type.count.pdf')
        matplotlib.pyplot.close(f)
    for ZNF143_HiChIP_loop_file in ZNF143_HiChIP_loop_files:
        ZNF143_loops = bioframe.read_table(ZNF143_HiChIP_loop_file, schema=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'name', 'score'], sep=' ')
        ZNF143_loops = standard_assign(ZNF143_loops, SBS, 'a1', 'SBS', 'pvalue', includes=['name', 'strand'], select_mode='min', cols1=('chrom1', 'start1', 'end1'))
        ZNF143_loops = standard_assign(ZNF143_loops, SBS, 'a2', 'SBS', 'pvalue', includes=['name', 'strand'], select_mode='min', cols1=('chrom2', 'start2', 'end2'))
        ZNF143_loops['a1_SBS_type'] = ['0' if not SBS_name else SBS_name[SBS_name.find('type')+4] for SBS_name in ZNF143_loops['a1_SBS_name']]
        ZNF143_loops['a2_SBS_type'] = ['0' if not SBS_name else SBS_name[SBS_name.find('type')+4] for SBS_name in ZNF143_loops['a2_SBS_name']]
        ZNF143_loops['SBS_type'] = ["".join(numpy.sort([a1_SBS_type,a2_SBS_type])) for a1_SBS_type, a2_SBS_type in zip(ZNF143_loops['a1_SBS_type'], ZNF143_loops['a2_SBS_type'])]
        f, ax = matplotlib.pyplot.subplots(figsize=(20,20))
        seaborn.histplot(data=ZNF143_loops, x='SBS_type', ax=ax)
        for patch in ax.patches:
            ax.text(patch.get_x() + patch.get_width()/2., patch.get_height(), f'{int(patch.get_height())}', ha="center", va="bottom")
        f.savefig(f'results/motif_type_count/ZNF143_peak.{os.path.basename(ZNF143_HiChIP_loop_file)}.motif.type.count.pdf')
        matplotlib.pyplot.close(f)




### ctcf peak enrichment
CTCF = bioframe.read_table(CTCF_peak_file, schema='bed')
CTCF['mid'] = (CTCF.end+CTCF.start)//2
ZNF143 = bioframe.read_table(ZNF143_peak_file, schema='bed')
# assign CBS to CTCF
CBS = pandas.concat([bioframe.read_table(f"{CTCF_pwm_file}.{os.path.basename(CTCF_peak_file)}.motifs", schema=['chrom', 'start', 'end', 'name', 'score', 'strand', 'pvalue', 'qvalue']) for CTCF_pwm_file in CTCF_pwm_files])
CTCF = standard_assign(CTCF, CBS, '', 'CBS', 'pvalue', includes=['name', 'strand'], select_mode='min', cols1=('chrom', 'start', 'end'), cols2=('chrom', 'start', 'end'))
# assign SBS to ZNF143
SBS = pandas.concat([bioframe.read_table(f"{ZNF143_pwm_file}.{os.path.basename(ZNF143_peak_file)}.motifs", schema=['chrom', 'start', 'end', 'name', 'score', 'strand', 'pvalue', 'qvalue']) for ZNF143_pwm_file in ZNF143_pwm_files])
ZNF143 = standard_assign(ZNF143, SBS, '', 'SBS', 'pvalue', includes=['name', 'strand'], select_mode='min', cols1=('chrom', 'start', 'end'), cols2=('chrom', 'start', 'end'))
# assign ZNF143 to CTCF
CTCF = standard_assign(CTCF, ZNF143, '', 'ZNF143', 'distance', includes=['SBS_name', 'SBS_strand', 'SBS_closest'], select_mode='min', cols1=('chrom', 'start', 'end'), cols2=('chrom', 'start', 'end'))
CTCF['SBS_closest'] = CTCF['ZNF143_closest'] + CTCF['ZNF143_SBS_closest']
# add hue columns
CTCF['ZNF143_has'] = ['ZNF143NA' if ZNF143_per_==0 else 'ZNF143' for ZNF143_per_ in CTCF['ZNF143_per_']]
CTCF['CBS_has'] = ['CBSNA' if CBS_per_==0 else 'CBS' for CBS_per_ in CTCF['CBS_per_']]
CTCF['SBS_has'] = ['SBSNA' if not ZNF143_SBS_name else 'SBS' for ZNF143_SBS_name in CTCF['ZNF143_SBS_name']]
CTCF['total_has'] = CTCF['ZNF143_has'] + '_' + CTCF['CBS_has'] + '_' + CTCF['SBS_has']
CTCF['CBS_SBS_strand'] = 'CBS' + CTCF['CBS_strand'].astype('str') + 'SBS' + CTCF['ZNF143_SBS_strand'].astype('str')
# compute signal
CTCF['CTCF_signal'] = bbi.stackup('ChIP_Nexus_rawdata/WT_CTCF_ChIP-nexus.adp.bam.bw', CTCF.chrom, CTCF.mid-100, CTCF.mid+100, bins=1, summary='mean')
CTCF['ZNF143_signal'] = bbi.stackup('ChIP_Nexus_rawdata/WT_ZNF143_ChIP-nexus.adp.bam.bw', CTCF.chrom, CTCF.mid-100, CTCF.mid+100, bins=1, summary='mean')

# save CTCF
CTCF.to_csv(path_or_buf='scrab_ucsc/CTCF.bed', sep='\t', header=True, index=False)

# draw hue scatterplot of CTCF_signal and ZNF143_signal
f, ax = matplotlib.pyplot.subplots()
seaborn.scatterplot(data=CTCF, x='CTCF_signal', y='ZNF143_signal', hue='total_has', alpha=0.2, s=10, label='CTCF peak types', ax=ax)
ax.set(xscale='log', yscale='log')
f.savefig('results/playwith/CTCF_peak_type.pdf')
matplotlib.pyplot.close(f)
# regression version scatter
f = seaborn.lmplot(data=CTCF, x='CTCF_signal', y='ZNF143_signal', hue='total_has')
f.set(xscale='log', yscale='log')
f.savefig('results/playwith/CTCF_peak_type_regress.pdf')
# bicluster according to CBS strand, SBS strand, ZNF143
ZNF143_data = (CTCF['ZNF143_per_']==0).astype('int64')
CBS_strand_data = [1 if CS=='+' else -1 if CS=='-' else 0 for CS in CTCF['CBS_strand']]
SBS_strand_data = [1 if SS=='+' else -1 if SS=='-' else 0 for SS in CTCF['ZNF143_SBS_strand']]
data=pandas.DataFrame({'ZNF143' : ZNF143_data, 'CBS_strand' : CBS_strand_data, 'SBS_strand' : SBS_strand_data})
sys.setrecursionlimit(100000)
f = seaborn.clustermap(data=data, standard_scale=1, cmap='bwr')
f.savefig('results/playwith/annotation_bicluster.pdf')
# motif shift hue plot, to get scatter version, remove kind="kde" and add marker="+", alpha=0.2
f = seaborn.jointplot(x=CTCF['CBS_closest'].astype('float64'), y=CTCF['SBS_closest'].astype('float64'), hue=CTCF['CBS_SBS_strand'], kind='kde', xlim=[-60,60], ylim=[-100,100], dropna=True)
f.savefig('results/playwith/CBS_SBS_orientation.pdf')

### extract motif from ZNF143_CBSNA_SBSNA.CTCF.bed
ZNF143_CBSNA_SBSNA = pandas.read_csv('ZNF143_CBSNA_SBSNA.CTCF.bed', sep='\t')
with open('tmp.fa', 'w') as f:
    for record in ZNF143_CBSNA_SBSNA.itertuples():
        f.write(f'>{record.name}\n{genome[record.chrom].ff.fetch(record.chrom, record.start, record.end)}\n')
# de novo motif
subprocess.check_output('meme tmp.fa -oc meme_out -dna -w 20 -nmotifs 20 -revcomp', shell=True)
subprocess.check_output('tomtom -oc tomtom_out meme_out/meme.txt JASPAR2022_CORE_non-redundant_pfms_meme.txt', shell=True)
# clean
os.remove('tmp.fa')

### display znf143
files = [f'COV_CHECK/{basename}_{rep}_HiC-Pro.bw' for basename in ['WT', 'K2'] for rep in ['rep1', 'rep2']]
pygenometracks.makeTracksFile.main(['--trackFiles'] + files + ['hg19.ncbiRefSeq.gtf', '--out', 'tracks.ini'])
with open('tracks.ini', 'r') as f:
    replace_string = f.read().replace('\nlabels = false\n', '\nlabels = true\n')
with open('tracks.ini', 'w') as f:
    f.write(replace_string)
pygenometracks.plotTracks.main(['--tracks', 'tracks.ini', '--region', 'chr11:9,480,000-9,552,000', '--outFileName', 'check_genome/hic_cover.pdf'])

### loop overlap
extend = 2000
fpkm_files = ['single_cut_RNA/cufflink/AID/genes.fpkm_tracking','single_cut_RNA/cufflink/AID-CD/genes.fpkm_tracking']
loops0_file = 'zmloops/AID_loop_ZZ.bedpe'
loops1_file = 'zmloops/AID-CD_loop_ZZ.bedpe'
loopss = [bioframe.read_table(file, schema=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'name', 'score'], sep=' ') for file in [loops0_file, loops1_file]]

# assign ctcf znf143 to anchor
CTCF = bioframe.read_table('ChIP_Nexus_results/AID_CTCF_ChIP-nexus_peaks.narrowPeak', schema='bed')
ZNF143 = bioframe.read_table('ChIP_Nexus_results/AID_ZNF143_ChIP-nexus_peaks.narrowPeak', schema='bed')
for i in range(len(loopss)):
    loopss[i] = standard_assign(loopss[i], ZNF143, 'a1', 'ZNF143', select_col='distance', includes=[], select_mode='min', cols1=('chrom1', 'start1', 'end1'))
    loopss[i] = standard_assign(loopss[i], ZNF143, 'a2', 'ZNF143', select_col='distance', includes=[], select_mode='min', cols1=('chrom2', 'start2', 'end2'))
    loopss[i] = standard_assign(loopss[i], CTCF, 'a1', 'CTCF', select_col='distance', includes=[], select_mode='min', cols1=('chrom1', 'start1', 'end1'))
    loopss[i] = standard_assign(loopss[i], CTCF, 'a2', 'CTCF', select_col='distance', includes=[], select_mode='min', cols1=('chrom2', 'start2', 'end2'))

# get overlap information
loopss[0], loopss[1] = get_overlap_loops(loopss[0], loopss[1])

# # get fpkm for tss
# tss = get_tss('hg19.ncbiRefSeq.gtf.gz')
# tss = tss.drop_duplicates(subset=['start'], ignore_index=True)
# coding = bioframe.read_table('hg19_gene_id_name.bed', schema=['chrom', 'start', 'end', 'strand', 'name', 'gene_id', 'type'])
# coding = coding[coding['type']=='protein_coding'].drop_duplicates(subset=['gene_id'], ignore_index=True)
# tss = pandas.merge(tss, coding['gene_id'], on='gene_id', how='inner')
# for fpkm_file in fpkm_files:
#     fpkm_name = os.path.basename(os.path.dirname(fpkm_file))
#     fpkm = pandas.read_csv(fpkm_file, sep='\t')[["gene_id",'FPKM']].groupby('gene_id')
#     sumFPKM = fpkm['FPKM'].sum().values
#     fpkm = fpkm.nth(0).reset_index()
#     fpkm['FPKM'] = sumFPKM
#     tss = pandas.merge(tss, fpkm[['gene_id', 'FPKM']], on='gene_id', how='left')
#     tss = tss.rename(columns={'FPKM' : f'FPKM_{fpkm_name}'}) 

# # assign fpkm to anchor
# for i in range(len(loopss)):
#     includes = [f'FPKM_{os.path.basename(os.path.dirname(fpkm_file))}' for fpkm_file in fpkm_files]
#     loopsecp = loopss[i][['start1','end1','start2','end2']]
#     loopss[i]['start1'], loopss[i]['end1'] = (loopss[i]['start1'] + loopss[i]['end1'])//2-extend, (loopss[i]['start1'] + loopss[i]['end1'])//2+extend
#     loopss[i]['start2'], loopss[i]['end2'] = (loopss[i]['start2'] + loopss[i]['end2'])//2-extend, (loopss[i]['start2'] + loopss[i]['end2'])//2+extend
#     loopss[i] = standard_assign(loopss[i], tss, 'a1', 'TSS', select_col='distance', includes=includes, select_mode='min', cols1=('chrom1', 'start1', 'end1'))
#     loopss[i] = standard_assign(loopss[i], tss, 'a2', 'TSS', select_col='distance', includes=includes, select_mode='min', cols1=('chrom2', 'start2', 'end2'))
#     loopss[i][['start1','end1','start2','end2']] = loopsecp

# assign ctcf to loop span
shrink = 2000
for i in range(len(loopss)):
    span_starts, span_ends = [], []
    for span_start, span_end in zip((loopss[i]['start1']+loopss[i]['end1'])//2+shrink, (loopss[i]['start2']+loopss[i]['end2'])//2-shrink):
        if span_start>=span_end:
            tmp = (span_start + span_end) // 2
            span_start, span_end = tmp, tmp+1
        span_starts.append(span_start)
        span_ends.append(span_end)
    loopss[i]['span_start'] = span_starts
    loopss[i]['span_end'] = span_ends
    loopss[i] = standard_assign(loopss[i], CTCF, 'span', 'ctcf', select_col='distance', includes=[], select_mode='min', cols1=('chrom1', 'span_start', 'span_end'))
    loopss[i]['span_ctcf'] = ['ctcf' if cps > 0 else 'noctcf' for cps in loopss[i]['ctcf_per_span']]

loops_names = ['AID', 'AID-CD']
arrs = {}
for i in range(len(loopss)):
    for gp in [['yes'], ['no','half']]:
        data = loopss[i][loopss[i]['over'].isin(gp)]
        f, ax = matplotlib.pyplot.subplots(figsize=(10,10))
        seaborn.histplot(data=data, x='span_ctcf', ax=ax)
        arrs[f"{loops_names[i]}.{gp}"] = []
        for patch in ax.patches:
            ax.text(patch.get_x() + patch.get_width()/2., patch.get_height(), f'{int(patch.get_height())}', ha="center", va="bottom")
            arrs[f"{loops_names[i]}.{gp}"].append(int(patch.get_height()))
        f.savefig(f'results/loop_counts/middle_ctcf.{loops_names[i]}.{gp}.pdf')
        matplotlib.pyplot.close(f)

title = ""
for key1, key2 in itertools.combinations(arrs.keys(), 2):
    pvalue = ljwlib.sundry.chi_square_two_multinomial_same(numpy.array(arrs[key1]), numpy.array(arrs[key2]))
    title += f"{key1}|{key2}={pvalue:.10f}\n"

for i in range(len(loopss)):
    for gp in [['yes'], ['no','half']]:
        data = loopss[i][loopss[i]['over'].isin(gp)]
        f, ax = matplotlib.pyplot.subplots(figsize=(10,10))
        seaborn.histplot(data=data, x='span_ctcf', stat='percent', ax=ax)
        arrs[f"{loops_names[i]}.{gp}"] = []
        for patch in ax.patches:
            ax.text(patch.get_x() + patch.get_width()/2., patch.get_height(), f'{int(patch.get_height())}', ha="center", va="bottom")
            arrs[f"{loops_names[i]}.{gp}"].append(int(patch.get_height()))
        ax.set(title=title)
        f.savefig(f'results/loop_counts/middle_ctcf.{loops_names[i]}.{gp}.percent.pdf')
        matplotlib.pyplot.close(f)

# loopss[0].to_csv(path_or_buf=f'results/loop_counts/{os.path.basename(loops0_file)}.expand', sep='\t', header=True, index=False)
# loopss[1].to_csv(path_or_buf=f'results/loop_counts/{os.path.basename(loops1_file)}.expand', sep='\t', header=True, index=False)

# for loop_file in [loops0_file, loops1_file]:
#     for sbs in [['yes'],['no'],['half'],['no','half'],['yes','no','half']]:
#         loops = pandas.read_csv(f'results/loop_counts/{os.path.basename(loop_file)}.expand', sep='\t')
#         loops = loops[loops['over'].isin(sbs)]
#         loops['FPKM_AID'] = loops['a1_TSS_FPKM_AID'].fillna(0.) + loops['a2_TSS_FPKM_AID'].fillna(0.)
#         loops['FPKM_AID-CD'] = loops['a1_TSS_FPKM_AID-CD'].fillna(0.) + loops['a2_TSS_FPKM_AID-CD'].fillna(0.)
#         loops = loops[loops['FPKM_AID']+loops['FPKM_AID-CD']>0]
#         pvalue = scipy.stats.ttest_rel(loops['FPKM_AID'], loops['FPKM_AID-CD'], alternative='greater').pvalue
#         f, ax = matplotlib.pyplot.subplots(figsize=(20,20))
#         seaborn.barplot(data=loops[['FPKM_AID','FPKM_AID-CD']], ax=ax)
#         for patch in ax.patches:
#             ax.text(patch.get_x() + patch.get_width()/2., patch.get_height(), f'{patch.get_height():.5f}', ha="center", va="bottom") 
#         ax.set(title=f'pvalue={pvalue}')
#         f.savefig(f'results/loop_counts/{os.path.basename(loop_file)}.{sbs}.fpkm.comp.pdf')
#         matplotlib.pyplot.close(f)

# for fil in [None, 'ZZ', 'ZO', 'OO']:
#     if not fil:
#         mask0, mask1 = [True] * len(loopss[0]), [True] * len(loopss[1])
#     else:
#         if fil=='ZZ':
#             mask0, mask1 = (loopss[0]['CTCF_per_a1']==0) & (loopss[0]['CTCF_per_a2']==0), (loopss[1]['CTCF_per_a1']==0) & (loopss[1]['CTCF_per_a2']==0)
#         elif fil=='OO':
#             mask0, mask1 = (loopss[0]['CTCF_per_a1']>0) & (loopss[0]['CTCF_per_a2']>0), (loopss[1]['CTCF_per_a1']>0) & (loopss[1]['CTCF_per_a2']>0)
#         elif fil=='ZO':
#             mask0, mask1 = (loopss[0]['CTCF_per_a1']>0) & (loopss[0]['CTCF_per_a2']==0) | (loopss[0]['CTCF_per_a1']==0) & (loopss[0]['CTCF_per_a2']>0), (loopss[1]['CTCF_per_a1']>0) & (loopss[1]['CTCF_per_a2']==0) | (loopss[1]['CTCF_per_a1']==0) & (loopss[1]['CTCF_per_a2']>0)
#         else:
#             raise Exception("fil must be ZZ, ZO, OO")
#     loops0 = loopss[0][mask0]
#     loops1 = loopss[1][mask1]
#     names = [f'{loops_names[0]}_only', f'{loops_names[0]}_half', 'over', f'{loops_names[1]}_only', f'{loops_names[1]}_half']
#     counts = [sum(loops0['over']=='no'), sum(loops0['over']=='half'), len(loops0.loc[loops0['over']=='yes', 'comp'].unique()), sum(loops1['over']=='no'), sum(loops1['over']=='half')]
#     data = pandas.DataFrame({'name' : names, 'count' : counts})
#     f, ax = matplotlib.pyplot.subplots()
#     seaborn.barplot(data=data, x='name', y='count', ax=ax)
#     for patch in ax.patches:
#         ax.text(patch.get_x() + patch.get_width()/2., patch.get_height(), f'{int(patch.get_height())}', ha="center", va="bottom") 
#     f.savefig(f'results/loop_counts/{os.path.basename(loops0_file)}_{os.path.basename(loops1_file)}.{fil}.overlap.count.pdf')
#     matplotlib.pyplot.close(f)

### run HiCPro
hichip_filepairs = [[f"hichip_rawdata/{file}" for file in pair] for pair in [['WT_ZNF143_HiChIP_R1.fq.gz','WT_ZNF143_HiChIP_R2.fq.gz'], ['AID_ZNF143_HiChIP_R1.fq.gz','AID_ZNF143_HiChIP_R2.fq.gz'], ['AID-CD_ZNF143_HiChIP_R1.fq.gz','AID-CD_ZNF143_HiChIP_R2.fq.gz'], ['WT_ZNF143_HiChIP_rep1_R1.fq.gz','WT_ZNF143_HiChIP_rep1_R2.fq.gz'], ['WT_ZNF143_HiChIP_rep2_R1.fq.gz','WT_ZNF143_HiChIP_rep2_R2.fq.gz'], ['AID_ZNF143_HiChIP_rep1_R1.fq.gz','AID_ZNF143_HiChIP_rep1_R2.fq.gz'], ['AID_ZNF143_HiChIP_rep2_R1.fq.gz','AID_ZNF143_HiChIP_rep2_R2.fq.gz'], ['AID-CD_ZNF143_HiChIP_rep1_R1.fq.gz','AID-CD_ZNF143_HiChIP_rep1_R2.fq.gz'], ['AID-CD_ZNF143_HiChIP_rep2_R1.fq.gz','AID-CD_ZNF143_HiChIP_rep2_R2.fq.gz']]]
hicpro_inpath = 'hichip_results/rawdata'
hicpro_outpath = 'hichip_results/hicpro'
hicpro_config = 'hichip_results/config-hicpro.txt'
ljwlib.softwares.run_HiC_Pro(hichip_filepairs, hicpro_inpath, hicpro_outpath, hicpro_config)

### run macs2
nexus_reps = [[f"ChIP_Nexus_rawdata/{file}" for file in rep] for rep in [['WT_CTCF_ChIP-nexus_rep1.fq.gz','WT_CTCF_ChIP-nexus_rep2.fq.gz'], ['AID_CTCF_ChIP-nexus_rep1.fq.gz','AID_CTCF_ChIP-nexus_rep2.fq.gz'], ['AID-CD_CTCF_ChIP-nexus_rep1.fq.gz','AID-CD_CTCF_ChIP-nexus_rep2.fq.gz'], ['WT_ZNF143_ChIP-nexus_rep1.fq.gz','WT_ZNF143_ChIP-nexus_rep2.fq.gz'], ['AID_ZNF143_ChIP-nexus_rep1.fq.gz','AID_ZNF143_ChIP-nexus_rep2.fq.gz'], ['AID-CD_ZNF143_ChIP-nexus_rep1.fq.gz','AID-CD_ZNF143_ChIP-nexus_rep2.fq.gz'], ['WT_RAD21_ChIP-nexus_rep1.fq.gz','WT_RAD21_ChIP-nexus_rep2.fq.gz']]]
macs2_outpath = 'ChIP_Nexus_results'
adapter = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
for rep in nexus_reps:
    ljwlib.softwares.run_macs2(rep, macs2_outpath, adapter, minlen=20, utrim=10, dtrim=0, threads=6, genome_index='/home/ljw/hg19_with_bowtie2_index/hg19', qvalue=0.001)

### run hichipper
hichipper_outpath = 'hichip_results'
resfrags_file = 'hichip_results/MboI.bed'
ljwlib.softwares.run_hichipper(hichipper_outpath, hicpro_outpath, macs2_outpath, resfrags_file)

### compare hichipper loops
hichipper_loop_files = [os.path.join("/media/ljw/0f85f321-d9f4-4fab-8153-c0f4a87a67271/ZNF143_loop/",file) for file in ["WT_ZNF143_HiChIP_rep1/WT_ZNF143_HiChIP_rep1.filt.intra.loop_counts.bedpe", "WT_ZNF143_HiChIP_rep2/WT_ZNF143_HiChIP_rep2.filt.intra.loop_counts.bedpe", "AID_ZNF143_HiChIP_rep1/AID_ZNF143_HiChIP_rep1.filt.intra.loop_counts.bedpe", "AID_ZNF143_HiChIP_rep2/AID_ZNF143_HiChIP_rep2.filt.intra.loop_counts.bedpe", "AID-CD_ZNF143_HiChIP_rep1/AID-CD_ZNF143_HiChIP_rep1.filt.intra.loop_counts.bedpe", "AID-CD_ZNF143_HiChIP_rep2/AID-CD_ZNF143_HiChIP_rep2.filt.intra.loop_counts.bedpe"]]
count_threses = [1,2,3]
nbdiges = [0,1,2]
dige = pandas.read_csv('MboI.bed', sep='\t')
for nbdige in nbdiges:
    mdige = merge_neighbour_bins(dige, nbdige)
    for count_thres in count_threses:
        f, axs = matplotlib.pyplot.subplots(figsize=(40,40), nrows=len(hichipper_loop_files), ncols=len(hichipper_loop_files), sharex=False)
        loopss = []
        for i in range(len(hichipper_loop_files)):
            loopss.append(bioframe.read_table(hichipper_loop_files[i], schema="bedpe", sep=" "))
            seaborn.histplot(data=loopss[-1], x='score', ax=axs[i,i])
            axs[i,i].tick_params(axis='x', rotation=90)
            axs[i,i].set(title=os.path.basename(hichipper_loop_files[i]).split('.')[0])
            loopss[-1] = loopss[-1][loopss[-1]['score']>=count_thres].reset_index(drop=True)
            loopss[-1] = align_to(loopss[-1], mdige, cols1=('chrom1', 'start1', 'end1'), cols2=('chrom', 'start', 'end'))
            loopss[-1] = align_to(loopss[-1], mdige, cols1=('chrom2', 'start2', 'end2'), cols2=('chrom', 'start', 'end'))
        for i in range(len(hichipper_loop_files)):
            for j in range(i+1,len(hichipper_loop_files)):
                loops1, loops2 = get_overlap_loops(loopss[i], loopss[j])
                df_over = pandas.DataFrame({'yes1' : sum(loops1['over']=='yes'), 'no1' : sum(loops1['over']=='no'), 'half1' : sum(loops1['over']=='half'), 'yes2' : sum(loops2['over']=='yes'), 'no2' : sum(loops2['over']=='no'), 'half2' : sum(loops2['over']=='half'), 'comp': len(loops1.loc[loops1['over']=='yes', 'comp'].unique())}, index=[0])
                title = f"{os.path.basename(hichipper_loop_files[i]).split('.')[0]}\n{os.path.basename(hichipper_loop_files[j]).split('.')[0]}"
                for pe in [[i,j],[j,i]]:
                    seaborn.barplot(data=df_over, ax=axs[pe[0],pe[1]])
                    for p in axs[pe[0],pe[1]].patches:
                        axs[pe[0],pe[1]].text(p.get_x() + p.get_width()/2., p.get_height(), f'{int(p.get_height())}', ha="center", va="bottom") 
                    axs[pe[0],pe[1]].tick_params(axis='x', rotation=90)
                    axs[pe[0],pe[1]].set(title=title)
        f.tight_layout()
        f.savefig(f"results/hichipper/loops_compare.{nbdige}.{count_thres}.pdf")
        matplotlib.pyplot.close(f)
            
### compare chip-nexus peaks
chip_nexus_peak_files = [f"ChIP_Nexus_results/{file}" for file in ["WT_ZNF143_ChIP-nexus_peaks.narrowPeak", "WT_ZNF143_ChIP-nexus_rep1_peaks.narrowPeak", "WT_ZNF143_ChIP-nexus_rep2_peaks.narrowPeak", "AID_ZNF143_ChIP-nexus_peaks.narrowPeak", "AID_ZNF143_ChIP-nexus_rep1_peaks.narrowPeak", "AID_ZNF143_ChIP-nexus_rep2_peaks.narrowPeak", "AID-CD_ZNF143_ChIP-nexus_peaks.narrowPeak", "AID-CD_ZNF143_ChIP-nexus_rep1_peaks.narrowPeak", "AID-CD_ZNF143_ChIP-nexus_rep2_peaks.narrowPeak"]]
f, axs = matplotlib.pyplot.subplots(figsize=(40,40), nrows=len(chip_nexus_peak_files), ncols=len(chip_nexus_peak_files), sharex=False)
for i in range(len(chip_nexus_peak_files)):
    peaks1 = bioframe.read_table(chip_nexus_peak_files[i], schema="bed")
    seaborn.kdeplot(x=peaks1['end']-peaks1['start'], ax=axs[i,i])
    axs[i,i].set(title=os.path.basename(chip_nexus_peak_files[i]).split('.')[0])
    for j in range(i+1,len(chip_nexus_peak_files)):
        peaks2 = bioframe.read_table(chip_nexus_peak_files[j], schema="bed")
        peaks1_peaks2 = standard_assign(peaks1, peaks2, '', 'assign', select_col='distance', includes=['name'], select_mode='min')
        peaks2_peaks1 = standard_assign(peaks2, peaks1, '', 'assign', select_col='distance', includes=['name'], select_mode='min')
        peaks1_only_count = sum(peaks1_peaks2['assign_name'].isna())
        peaks2_only_count = sum(peaks2_peaks1['assign_name'].isna())
        comp_count = len(bioframe.merge(pandas.concat([peaks1[~peaks1_peaks2['assign_name'].isna()],peaks2[~peaks2_peaks1['assign_name'].isna()]]), min_dist=None))
        df_over = pandas.DataFrame({f'only1' : peaks1_only_count, f'over1' : len(peaks1)-peaks1_only_count, f'only2' : peaks2_only_count, f'over2' : len(peaks2)-peaks2_only_count, 'comp' : comp_count}, index=[0])
        title = f"{os.path.basename(chip_nexus_peak_files[i]).split('.')[0]}\n{os.path.basename(chip_nexus_peak_files[j]).split('.')[0]}"
        for pe in [[i,j],[j,i]]:
            seaborn.barplot(data=df_over, ax=axs[pe[0],pe[1]])
            for p in axs[pe[0],pe[1]].patches:
                axs[pe[0],pe[1]].text(p.get_x() + p.get_width()/2., p.get_height(), f'{int(p.get_height())}', ha="center", va="bottom") 
            axs[pe[0],pe[1]].tick_params(axis='x', rotation=90)
            axs[pe[0],pe[1]].set(title=title)
f.tight_layout()
f.savefig(f"results/hichipper/peaks_compare.pdf")
matplotlib.pyplot.close(f)

### count validpairs and create MboI matrices
vpfiles = ['hic_matrices/WT-rep1.allValidPairs.gz', 'hic_matrices/WT-rep2.allValidPairs.gz', 'hic_matrices/ZO1-rep1.allValidPairs.gz', 'hic_matrices/ZO1-rep2.allValidPairs.gz', 'hic_matrices/ZO2-rep1.allValidPairs.gz', 'hic_matrices/ZO2-rep2.allValidPairs.gz', 'hic_matrices/AID-rep1.allValidPairs.gz', 'hic_matrices/AID-rep2.allValidPairs.gz', 'hic_matrices/AID-CD-rep1.allValidPairs.gz', 'hic_matrices/AID-CD-rep2.allValidPairs.gz', 'hic_matrices/AID-ZO-rep1.allValidPairs.gz', 'hic_matrices/AID-ZO-rep2.allValidPairs.gz', 'hic_matrices/AID-ZOCD-rep1.allValidPairs.gz', 'hic_matrices/AID-ZOCD-rep2.allValidPairs.gz']
dige = pandas.read_csv("MboI.bed", sep='\t')
lns = {}
for vpfile in vpfiles:
    cool_uri = f'{os.path.basename(vpfile).split(".")[0]}.MboI.cool'
    VPs = bioframe.read_table(vpfile, schema=['chrom1', 'start1', 'chrom2', 'start2'], usecols=[1,2,4,5])
    lns[cool_uri] = len(VPs)
    VPs['end1'] = VPs['start1']+1
    VPs['bin1_id'] = bioframe.closest(VPs, dige, return_input=False, return_distance=False, return_index=True, cols1=['chrom1','start1','end1']).set_index('index')['index_'].astype('int64')
    VPs = VPs.drop(columns = ['chrom1','start1','end1'])
    VPs['end2'] = VPs['start2']+1
    VPs['bin2_id'] = bioframe.closest(VPs, dige, return_input=False, return_distance=False, return_index=True, cols1=['chrom2','start2','end2']).set_index('index')['index_'].astype('int64')
    VPs = VPs.drop(columns = ['chrom2','start2','end2'])
    VPs['count'] = 1
    mask_rev = VPs['bin1_id'] > VPs['bin2_id']
    VPs.loc[mask_rev, 'bin1_id'], VPs.loc[mask_rev, 'bin2_id'] = VPs.loc[mask_rev, 'bin2_id'], VPs.loc[mask_rev, 'bin1_id']
    VPs = VPs.groupby(['bin1_id','bin2_id'])
    sumcount = VPs['count'].sum().values
    VPs = VPs.nth(0).reset_index()
    VPs['count'] = sumcount
    VPs = VPs.sort_values(['bin1_id','bin2_id']).reset_index(drop=True)
    cooler.create_cooler(f'hic_matrices_MboI/{cool_uri}', dige, VPs)
with open('hic_matrices/allValidPairs.ln', 'wb') as f:
    pickle.dump(lns, f)


### hist to MboI bins
def get_group_greater_pvalues(data=None, filesWT=None, filesKO=None):
    if not isinstance(data,pandas.DataFrame):
        return ['ttest_ind','ranksums','nrm']
    WTbases = [os.path.basename(file) for file in filesWT]
    KObases = [os.path.basename(file) for file in filesKO]
    WT_counts = pandas.concat([data[base] for base in filesWT])
    KO_counts = pandas.concat([data[base] for base in filesKO])
    pvalues = {}
    pvalues['ttest_ind'] = scipy.stats.ttest_ind(WT_counts, KO_counts, axis=0, nan_policy='omit', alternative='greater').pvalue
    pvalues['ranksums'] = scipy.stats.ranksums(WT_counts, KO_counts, alternative='greater').pvalue
    pvalues['nrm'] = ljwlib.sundry.nested_repeated_measure(data[WTbases+KObases], [WTbases,KObases])
    return pvalues

# loop_groups = [['hichip_results/WT_ZNF143_HiChIP_rep1/WT_ZNF143_HiChIP_rep1.filt.intra.loop_counts.bedpe','/hichip_results/WT_ZNF143_HiChIP_rep2/WT_ZNF143_HiChIP_rep2.filt.intra.loop_counts.bedpe'], 'hichip_results/WT_ZNF143_HiChIP/WT_ZNF143_HiChIP.filt.intra.loop_counts.bedpe']
loop_groups = [['hichip_results/WT_ZNF143_HiChIP_rep1/WT_ZNF143_HiChIP_rep1.filt.intra.loop_counts.bedpe','hichip_results/WT_ZNF143_HiChIP_rep2/WT_ZNF143_HiChIP_rep2.filt.intra.loop_counts.bedpe']]

with open('hic_matrices/allValidPairs.ln', 'rb') as f:
    lns = pickle.load(f) 
f, ax = matplotlib.pyplot.subplots(figsize=(20, 10))
seaborn.barplot(data=pandas.DataFrame(lns, index=[0]), orient='v', ax=ax)
for p in ax.patches:
    ax.text(p.get_x() + p.get_width()/2., p.get_height(), f'{int(p.get_height())}', ha="center", va="bottom") 
ax.tick_params(axis='x', rotation=90)
f.tight_layout()
f.savefig(f'results/playwith/allValidPairs.count.pdf')
matplotlib.pyplot.close(f)

dige = pandas.read_csv('MboI.bed', sep='\t')
CTCF = bioframe.read_table('ChIP_Nexus_results/WT_CTCF_ChIP-nexus_peaks.narrowPeak', schema='bed')
ZNF143 = bioframe.read_table('ChIP_Nexus_results/WT_ZNF143_ChIP-nexus_peaks.narrowPeak', schema='bed')
inter = True
count_threses = range(1,32,3)
nbdiges = range(0,22,3)
for loop_group in loop_groups:
    # get loop
    if isinstance(loop_group, list):
        loopss = []
        for file in loop_group:
            loopss.append(bioframe.read_table(file, schema='bedpe', sep=' '))
        loops1, loops2 = get_overlap_loops(loopss[0], loopss[1])
        if inter:
            loops1 = loops1[loops1['over']=='yes'].reset_index(drop=True)
            loops2 = loops2[loops2['over']=='yes'].reset_index(drop=True)
        loops_all = pandas.concat([loops1, loops2]).reset_index(drop=True).groupby('comp')
        chroms1, chroms2, scores, starts1, ends1, starts2, ends2 = loops_all['chrom1'].nth(0).values.flatten(), loops_all['chrom2'].nth(0).values.flatten(), loops_all['score'].sum().values.flatten(), loops_all['start1'].min().values.flatten(), loops_all['end1'].max().values.flatten(), loops_all['start2'].min().values.flatten(), loops_all['end2'].max().values.flatten()
        loops_all = pandas.DataFrame({'chrom1' : chroms1, 'start1' : starts1, 'end1' : ends1, 'chrom2' : chroms2, 'start2' : starts2, 'end2' : ends2, 'name' : '.', 'score' : scores})
        out_name = f"{os.path.commonprefix([os.path.basename(loop_group[0]),os.path.basename(loop_group[1])])}.{inter}"
    else:
        loops_all = bioframe.read_table(loop_group, schema='bedpe', sep=' ')
        out_name = os.path.basename(loop_group).split('.')[0]
    # attach CTCF ZNF143
    loops_all = standard_assign(loops_all, ZNF143, 'a1', 'ZNF143', select_col='distance', includes=['name'], select_mode='min', cols1=('chrom1', 'start1', 'end1'))
    loops_all = standard_assign(loops_all, ZNF143, 'a2', 'ZNF143', select_col='distance', includes=['name'], select_mode='min', cols1=('chrom2', 'start2', 'end2'))
    loops_all = standard_assign(loops_all, CTCF, 'a1', 'CTCF', select_col='distance', includes=['name'], select_mode='min', cols1=('chrom1', 'start1', 'end1'))
    loops_all = standard_assign(loops_all, CTCF, 'a2', 'CTCF', select_col='distance', includes=['name'], select_mode='min', cols1=('chrom2', 'start2', 'end2'))
    # filter SBS only loops
    boolESES = (~loops_all['a1_ZNF143_name'].isna()) & (~loops_all['a2_ZNF143_name'].isna())
    boolNCNC = (loops_all['a1_CTCF_name'].isna()) & (loops_all['a2_CTCF_name'].isna())
    loops_all = loops_all[boolESES & boolNCNC].reset_index(drop=True)
    # align to emzyme fragment
    loops_all['mid1s'] = (loops_all['start1'] + loops_all['end1'])//2
    loops_all['mid1e'] = loops_all['mid1s'] +1
    loops_all['mid2s'] = (loops_all['start2'] + loops_all['end2'])//2
    loops_all['mid2e'] = loops_all['mid2s'] +1
    loops_all['bin1_id'] =  bioframe.closest(loops_all, dige, return_input=False, return_distance=False, return_index=True, cols1=['chrom1','mid1s','mid1e'])['index_'].astype('int64')
    loops_all['bin2_id'] =  bioframe.closest(loops_all, dige, return_input=False, return_distance=False, return_index=True, cols1=['chrom2','mid2s','mid2e'])['index_'].astype('int64')
    with open(f"results/playwith/{out_name}.pvalues", 'w') as pf:
        pf.write("nbdige\tcount_thres")
        for key in get_group_greater_pvalues():
            for pre in ['WT','AID_merge','AID','AID_CD']:
                pf.write(f"\t{pre}_{key}")
        pf.write("\n")
        for nbdige in nbdiges:
            # count
            MboI_cool_files = ['hic_matrices_MboI/WT-rep1.MboI.cool', 'hic_matrices_MboI/WT-rep2.MboI.cool', 'hic_matrices_MboI/ZO1-rep1.MboI.cool', 'hic_matrices_MboI/ZO1-rep2.MboI.cool', 'hic_matrices_MboI/ZO2-rep1.MboI.cool', 'hic_matrices_MboI/ZO2-rep2.MboI.cool', 'hic_matrices_MboI/AID-rep1.MboI.cool', 'hic_matrices_MboI/AID-rep2.MboI.cool', 'hic_matrices_MboI/AID-CD-rep1.MboI.cool', 'hic_matrices_MboI/AID-CD-rep2.MboI.cool', 'hic_matrices_MboI/AID-ZO-rep1.MboI.cool', 'hic_matrices_MboI/AID-ZO-rep2.MboI.cool', 'hic_matrices_MboI/AID-ZOCD-rep1.MboI.cool', 'hic_matrices_MboI/AID-ZOCD-rep2.MboI.cool']
            jobpacks = joblib.Parallel(n_jobs=len(MboI_cool_files))(joblib.delayed(count_loop_valid)(file, loops_all, nbdige, len(dige)) for file in MboI_cool_files)
            for i, file in enumerate(MboI_cool_files):
                loops_all[os.path.basename(file)] = jobpacks[i][0]
                loops_all[f"{os.path.basename(file)}.a1cov"] = jobpacks[i][1]
                loops_all[f"{os.path.basename(file)}.a2cov"] = jobpacks[i][2]
            for file in MboI_cool_files:
                loops_all[f"{os.path.basename(file)}.norm"] = loops_all[os.path.basename(file)]/lns[os.path.basename(file)]*1e6
            loops_all.to_csv(f'results/playwith/{out_name}.{nbdige}.loops.count', sep='\t', index=False)
            # loops_all = pandas.read_csv(f'results/playwith/{out_name}.{nbdige}.loops.count', sep='\t')

            for count_thres in count_threses:
                # filter by count_thres
                loops = loops_all[loops_all['score']>=count_thres].reset_index(drop=True)
                # normalize
                data = loops[[os.path.basename(file) for file in MboI_cool_files]]
                for file in MboI_cool_files:
                    data[os.path.basename(file)] = data[os.path.basename(file)]/lns[os.path.basename(file)]*1e6
                # draw
                # f, ax = matplotlib.pyplot.subplots(figsize=(20, 10))
                # seaborn.violinplot(data=data, orient='v', ax=ax)
                # ax.set(title=f"{len(loops)}")
                # ax.tick_params(axis='x', rotation=90)
                # f.tight_layout()
                # f.savefig(f'results/playwith/{out_name}.{nbdige}.{count_thres}.MboI.count.volin.pdf')
                # matplotlib.pyplot.close(f)

                # f, ax = matplotlib.pyplot.subplots(figsize=(20, 10))
                # seaborn.boxplot(data=data, orient='v', ax=ax)
                # ax.set(title=f"{len(loops)}")
                # ax.tick_params(axis='x', rotation=90)
                # f.tight_layout()
                # f.savefig(f'results/playwith/{out_name}.{nbdige}.{count_thres}.MboI.count.box.pdf')
                # matplotlib.pyplot.close(f)

                MboI_cool_files_groups = [['hic_matrices_MboI/WT-rep1.MboI.cool','hic_matrices_MboI/WT-rep2.MboI.cool'], ['hic_matrices_MboI/ZO1-rep1.MboI.cool','hic_matrices_MboI/ZO1-rep2.MboI.cool', 'hic_matrices_MboI/ZO2-rep1.MboI.cool','hic_matrices_MboI/ZO2-rep2.MboI.cool'], ['hic_matrices_MboI/AID-rep1.MboI.cool','hic_matrices_MboI/AID-rep2.MboI.cool'], ['hic_matrices_MboI/AID-CD-rep1.MboI.cool','hic_matrices_MboI/AID-CD-rep2.MboI.cool'], ['hic_matrices_MboI/AID-ZO-rep1.MboI.cool','hic_matrices_MboI/AID-ZO-rep2.MboI.cool'], ['hic_matrices_MboI/AID-ZOCD-rep1.MboI.cool','hic_matrices_MboI/AID-ZOCD-rep2.MboI.cool']]
                width = 5
                sideunit = 0.03
                unit = (1-sideunit)/len(MboI_cool_files_groups)
                topunit = 0.3
                botunit = 1-topunit
                heiratio = len(loops)/len(loops_all)
                
                heatss = []
                heatflatten = numpy.array([])
                for i, group in enumerate(MboI_cool_files_groups):
                    heatss.append([])
                    for j, acov in enumerate(['a1cov','a2cov']):
                        heat = 0
                        for file in group:
                            heat += numpy.vstack(loops[f"{os.path.basename(file)}.{acov}"])
                        heat /= sum([lns[os.path.basename(file)] for file in group])
                        heat *= 1e6
                        maxheat = max(maxheat, numpy.max(heat))
                        heatss[-1].append(heat)
                        heatflatten = numpy.concatenate([heatflatten, heat.flatten()])
                norm=matplotlib.colors.Normalize(vmin=0,vmax=0.8*maxheat)

                maxcov = -numpy.inf
                maxheat = numpy.percentile(heatflatten, 80)
                f, axs = matplotlib.pyplot.subplots(figsize=(width*len(MboI_cool_files_groups), 20), ncols=len(MboI_cool_files_groups)+1, nrows=4)
                for i, group in enumerate(MboI_cool_files_groups):
                    gname = os.path.commonprefix([os.path.basename(file) for file in group])
                    for j, acov in enumerate(['a1cov','a2cov']):
                        heatcov = numpy.sum(heatss[i][j],axis=0)
                        maxcov = max(maxcov, numpy.max(heatcov))
                        axs[2*j,i].plot(numpy.arange(-nbdige,nbdige+1), heatcov, marker='.')
                        im = axs[2*j+1,i].matshow(heatss[i][j], norm=norm, extent=(-nbdige-0.5,nbdige+0.5,len(loops),0), cmap='bwr', aspect='auto')
                        axs[2*j,i].set_position([unit*0.1+unit*i,(botunit+topunit*0.1+j)*0.5,unit*0.8,topunit*0.8*0.5])
                        axs[2*j+1,i].set(xticks=[-nbdige,0,nbdige], yticks=[0, len(loops)], title=f"{gname}.{acov}")
                        axs[2*j+1,i].set_position([unit*0.1+unit*i,(botunit*(0.9-0.8*heiratio)+j)*0.5,unit*0.8,botunit*0.8*heiratio*0.5])
                        axs[2*j+1,i].xaxis.set_ticks_position('bottom')
                for i, group in enumerate(MboI_cool_files_groups):
                    for j, acov in enumerate(['a1cov','a2cov']):
                        axs[2*j,i].set(ylim=[0,int(maxcov*1.2)],yticks=[0,maxcov*1.2],xlim=[-nbdige-0.5,nbdige+0.5],xticks=numpy.arange(-nbdige,nbdige+1))
                axs[1,-1].set_position([unit*len(MboI_cool_files_groups),botunit*0.1*0.5,sideunit*0.3,botunit*0.8*0.5])
                f.colorbar(im, cax=axs[1,-1])
                for i in [0,2,3]:
                    axs[i,-1].set(visible=False)
                f.savefig(f'results/playwith/{out_name}.{nbdige}.{count_thres}.2Dheatmap.pdf')
                matplotlib.pyplot.close(f)

                f, ax = matplotlib.pyplot.subplots(figsize=(20, 10))
                seaborn.barplot(data=data, orient='v', ax=ax)
                for p in ax.patches:
                    ax.text(p.get_x() + p.get_width()/2., p.get_height(), f'{p.get_height():.10f}', ha="center", va="bottom") 
                ax.tick_params(axis='x', rotation=90)
                ax.set(title=f"{len(loops)}")
                f.tight_layout()
                f.savefig(f'results/playwith/{out_name}.{nbdige}.{count_thres}.MboI.count.barerr.pdf')
                matplotlib.pyplot.close(f)

                dataargr = pandas.DataFrame({'sample' : ['WT','WT','ZO','ZO','ZO','ZO','AID','AID','AID','AID','AID-ZO','AID-ZO','AID-ZO','AID-ZO'], 'count' : [sum(data[os.path.basename(file)]) for file in MboI_cool_files]})
                WTpvalues = get_group_greater_pvalues(data=data, filesWT=['WT-rep1.MboI.cool','WT-rep2.MboI.cool'], filesKO=['ZO1-rep1.MboI.cool','ZO1-rep2.MboI.cool','ZO2-rep1.MboI.cool','ZO2-rep2.MboI.cool'])
                AIDmergepvalues = get_group_greater_pvalues(data=data, filesWT=['AID-rep1.MboI.cool','AID-rep2.MboI.cool','AID-CD-rep1.MboI.cool','AID-CD-rep2.MboI.cool'], filesKO=['AID-ZO-rep1.MboI.cool','AID-ZO-rep2.MboI.cool','AID-ZOCD-rep1.MboI.cool','AID-ZOCD-rep2.MboI.cool'])
                f, ax = matplotlib.pyplot.subplots(figsize=(20, 10))
                seaborn.barplot(data=dataargr, x='sample', y ='count', orient='v', ax=ax)
                for p in ax.patches: 
                    ax.text(p.get_x() + p.get_width()/2., p.get_height(), f'{p.get_height():.2f}', ha="center", va="bottom")
                ylim = ax.get_ylim()
                psplit = (ylim[1] - ylim[0])/30
                for i, pvalues in enumerate([WTpvalues,AIDmergepvalues]):
                    for j, key in enumerate(pvalues.keys()):
                        ax.text((ax.patches[2*i].get_x()+ax.patches[2*i].get_width()/2+ax.patches[2*i+1].get_x()+ax.patches[2*i+1].get_width()/2)/2, max(ax.patches[2*i].get_height(),ax.patches[2*i+1].get_height())+j*psplit, f'{key}={pvalues[key]}', ha="center", va="bottom")
                ax.set(title=f"{len(loops)}")
                f.tight_layout()
                f.savefig(f'results/playwith/{out_name}.{nbdige}.{count_thres}.MboI.count.barargr.merge.pdf')
                matplotlib.pyplot.close(f)

                dataargr = pandas.DataFrame({'sample' : ['WT','WT','ZO','ZO','ZO','ZO','AID','AID','AID-CD','AID-CD','AID-ZO','AID-ZO','AID-ZOCD','AID-ZOCD'], 'count' : [sum(data[os.path.basename(file)]) for file in MboI_cool_files]})
                AIDpvalues = get_group_greater_pvalues(data=data, filesWT=['AID-rep1.MboI.cool','AID-rep2.MboI.cool'], filesKO=['AID-ZO-rep1.MboI.cool','AID-ZO-rep2.MboI.cool'])
                AIDCDpvalues = get_group_greater_pvalues(data=data, filesWT=['AID-CD-rep1.MboI.cool','AID-CD-rep2.MboI.cool'], filesKO=['AID-ZOCD-rep1.MboI.cool','AID-ZOCD-rep2.MboI.cool'])
                f, ax = matplotlib.pyplot.subplots(figsize=(20, 10))
                seaborn.barplot(data=dataargr, x='sample', y='count', orient='v', hue_order=['WT','ZO','AID','AID-ZO','AID-CD','AID-ZOCD'], ax=ax)
                for p in ax.patches: 
                    ax.text(p.get_x() + p.get_width()/2., p.get_height(), f'{p.get_height():.2f}', ha="center", va="bottom")
                ylim = ax.get_ylim()
                psplit = (ylim[1] - ylim[0])/30
                for i, pvalues in enumerate([WTpvalues,AIDpvalues,AIDCDpvalues]):
                    for j, key in enumerate(pvalues.keys()):
                        ax.text((ax.patches[2*i].get_x()+ax.patches[2*i].get_width()/2+ax.patches[2*i+1].get_x()+ax.patches[2*i+1].get_width()/2)/2, max(ax.patches[2*i].get_height(),ax.patches[2*i+1].get_height())+j*psplit, f'{key}={pvalues[key]}', ha="center", va="bottom")
                ax.set(title=f"{len(loops)}")
                f.tight_layout()
                f.savefig(f'results/playwith/{out_name}.{nbdige}.{count_thres}.MboI.count.barargr.pdf')
                matplotlib.pyplot.close(f)

                pf.write(f"{nbdige}\t{count_thres}")
                for key in get_group_greater_pvalues():
                    for pvalues in [WTpvalues, AIDmergepvalues, AIDpvalues, AIDCDpvalues]:
                        pf.write(f"\t{pvalues[key]}")
                pf.write("\n")

### chi-square two multinomial distributions are the same
arr1, arr2 = numpy.array([9226, 12571, 5915]), numpy.array([1609, 6442, 7199])
chi2 = ljwlib.sundry.chi_square_two_multinomial_same(arr1, arr2)

### pcdh cluster RNA-seq
f, axs = matplotlib.pyplot.subplots(figsize=(10,10), nrows=3, ncols=1, sharex=True)
rpkmss = {}
for i, name in enumerate(['WT', 'ZO']):
    rpkmss[name] = []
    for suf in ['-rep1', '-rep2']:
        rpkm = pandas.read_csv(f'single_cut_RNA/cufflink/{name}{suf}/genes.fpkm_tracking', sep='\t')
        rpkm = rpkm[[(tid.startswith("PCDHA") or tid.startswith("PCDHB") or tid.startswith("PCDHG")) and not tid.endswith("P") for tid in rpkm['tracking_id']]]
        rpkmss[name].append(rpkm)
    rpkmss[name] = pandas.concat(rpkmss[name])
    seaborn.barplot(data=rpkmss[name], x='tracking_id', y='FPKM', color='k', ax=axs[i])
    for p in axs[i].patches:
        axs[i].text(p.get_x() + p.get_width()/2., p.get_height(), f'{p.get_height(): .2f}', rotation=90, ha="center", va="bottom") 
    axs[i].set(title=name)
    # axs[i].tick_params(axis='x', rotation=90)
rpkmdiff = pandas.DataFrame({'tracking_id' : rpkmss['WT']['tracking_id'].values, 'FPKM' : rpkmss['ZO']['FPKM'].values-rpkmss['WT']['FPKM'].values})
seaborn.barplot(data=rpkmdiff, x='tracking_id', y='FPKM', color='k', ax=axs[2])
for p, lab in zip(axs[2].patches, axs[2].get_xticklabels()):
    hc, hp = p.get_height(), 0
    vac, vap = "bottom", "top"
    if p.get_height()<0:
        hc, hp = hp, hc
    axs[2].text(p.get_x() + p.get_width()/2., hc, f'{p.get_height(): .5f}', rotation=90, ha="center", va=vac)
    pvalue = scipy.stats.ttest_ind(rpkmss['WT'].loc[rpkmss['WT']['tracking_id']==lab._text,'FPKM'].values, rpkmss['ZO'].loc[rpkmss['ZO']['tracking_id']==lab._text,'FPKM'].values, axis=0, nan_policy='omit', alternative='two-sided').pvalue
    axs[2].text(p.get_x() + p.get_width()/2., hp, f'{pvalue: .5f}', rotation=90, ha="center", va=vap)
axs[2].set(title='ZO-WT')
axs[2].tick_params(axis='x', rotation=90)
f.tight_layout()
f.savefig(f'single_cut_RNA/fpkm.pdf')
matplotlib.pyplot.close(f)



### SBS loop based analysis
extend = 2000
loop_files = ['results/express_and_loop_strength/AID_loop_ZZ.bedpe', 'results/express_and_loop_strength/AID-CD_loop_ZZ.bedpe']
fpkm_filess = [['single_cut_RNA/cufflink/AID/genes.fpkm_tracking','single_cut_RNA/cufflink/AID-ZO/genes.fpkm_tracking','single_cut_RNA/cufflink/AID-CD/genes.fpkm_tracking','single_cut_RNA/cufflink/AID-ZOCD/genes.fpkm_tracking']] * 2
hic_filesss = [[['hic_matrices_MboI/AID-rep1.MboI.cool', 'hic_matrices_MboI/AID-rep2.MboI.cool'],['hic_matrices_MboI/AID-ZO-rep1.MboI.cool', 'hic_matrices_MboI/AID-ZO-rep2.MboI.cool'],['hic_matrices_MboI/AID-CD-rep1.MboI.cool', 'hic_matrices_MboI/AID-CD-rep2.MboI.cool'],['hic_matrices_MboI/AID-ZOCD-rep1.MboI.cool', 'hic_matrices_MboI/AID-ZOCD-rep2.MboI.cool']]] * 2
H3K4me3_bwfiless = [['???','???']] * 2
H3K27ac_files = ['results/express_and_loop_strength/H3K27ac_peaks.broadPeak'] * 2


for loop_file, fpkm_files, hic_filess, H3K4me3_bwfiles, H3K27ac_file in zip(loop_files, fpkm_filess, hic_filesss, H3K4me3_bwfiless, H3K27ac_files):
    loop = bioframe.read_table(loop_file, schema=['chrom1','start1','end1','chrom2','start2','end2','score'])
    loopsecp = loop[['start1','end1','start2','end2']]
    loop['start1'], loop['end1'] = (loop['start1'] + loop['end1'])//2-extend, (loop['start1'] + loop['end1'])//2+extend
    loop['start2'], loop['end2'] = (loop['start2'] + loop['end2'])//2-extend, (loop['start2'] + loop['end2'])//2+extend

    ## count hic
    # dige = pandas.read_csv("MboI.bed", sep='\t')
    # loop['mid1s'] = (loop['start1'] + loop['end1'])//2
    # loop['mid1e'] = loop['mid1s'] +1
    # loop['mid2s'] = (loop['start2'] + loop['end2'])//2
    # loop['mid2e'] = loop['mid2s'] +1
    # loop['bin1_id'] =  bioframe.closest(loop, dige, return_input=False, return_distance=False, return_index=True, cols1=['chrom1','mid1s','mid1e'])['index_'].astype('int64')
    # loop['bin2_id'] =  bioframe.closest(loop, dige, return_input=False, return_distance=False, return_index=True, cols1=['chrom2','mid2s','mid2e'])['index_'].astype('int64')

    # nbdige = 12
    # with open('hic_matrices/allValidPairs.ln', 'rb') as f:
    #     lns = pickle.load(f) 
    # MboI_cool_files = [hic_file for hic_files in hic_filess for hic_file in hic_files]
    # jobpacks = joblib.Parallel(n_jobs=len(MboI_cool_files))(joblib.delayed(count_loop_valid)(hic_file, loop, nbdige, len(dige)) for hic_file in MboI_cool_files)
    # total_reads = {}
    # for hic_files in hic_filess:
    #     hic_name = os.path.commonprefix([os.path.basename(hic_file) for hic_file in hic_files])
    #     loop[hic_name] = 0
    #     total_reads[hic_name] = 0
    #     for i, hic_file in enumerate(MboI_cool_files):
    #         if hic_file in hic_files:
    #             loop[hic_name] += jobpacks[i][0]
    #             total_reads[hic_name] += lns[os.path.basename(file)]
    #     loop[f'{hic_name}_norm'] = loop[hic_name]/total_reads[hic_name]*1e6

    # start assign annoations here
    # flank = 2000
    # for H3K4me3_bwfile in H3K4me3_bwfiles:
    #     loop[f'a1_{H3K4me3_bwfile}'] = bbi.stackup(H3K4me3_bwfile, loop['chrom1'], loop['start1']-flank, loop['end1']+flank, bins=1, summary='mean').flatten()
    #     loop[f'a2_{H3K4me3_bwfile}'] = bbi.stackup(H3K4me3_bwfile, loop['chrom2'], loop['start2']-flank, loop['end2']+flank, bins=1, summary='mean').flatten()
    #     loop[H3K4me3_bwfile] = loop[f'a1_{H3K4me3_bwfile}'] + loop[f'a2_{H3K4me3_bwfile}']

    tss = get_tss('hg19.ncbiRefSeq.gtf.gz')
    tss = tss.drop_duplicates(subset=['start'], ignore_index=True)
    coding = bioframe.read_table('hg19_gene_id_name.bed', schema=['chrom', 'start', 'end', 'strand', 'name', 'gene_id', 'type'])
    coding = coding[coding['type']=='protein_coding'].drop_duplicates(subset=['gene_id'], ignore_index=True)
    tss = pandas.merge(tss, coding['gene_id'], on='gene_id', how='inner')
    for fpkm_file in fpkm_files:
        fpkm_name = os.path.basename(os.path.dirname(fpkm_file))
        fpkm = pandas.read_csv(fpkm_file, sep='\t')[["gene_id",'FPKM']].groupby('gene_id')
        sumFPKM = fpkm['FPKM'].sum().values
        fpkm = fpkm.nth(0).reset_index()
        fpkm['FPKM'] = sumFPKM
        tss = pandas.merge(tss, fpkm[['gene_id', 'FPKM']], on='gene_id', how='left')
        tss = tss.rename(columns={'FPKM' : f'FPKM_{fpkm_name}'}) 

    for anchor, cols1 in zip(['a1','a2'],[('chrom1', 'start1', 'end1'),('chrom2', 'start2', 'end2')]):
        includes = [f'FPKM_{os.path.basename(os.path.dirname(fpkm_file))}' for fpkm_file in fpkm_files]
        loop = standard_assign(loop, tss, anchor, 'TSS', select_col='distance', includes=includes, select_mode='min', cols1=cols1)

    for fpkm_file in fpkm_files:
        fpkm_name = os.path.basename(os.path.dirname(fpkm_file))
        loop[f'FPKM_{fpkm_name}'] = (loop[f'a1_TSS_FPKM_{fpkm_name}'].fillna(0.0) + loop[f'a2_TSS_FPKM_{fpkm_name}'].fillna(0.0))/2

    H3K27ac = bioframe.read_table(H3K27ac_file, schema='bed')
    loop = standard_assign(loop, H3K27ac, 'a1', 'H3K27ac', select_col='distance', includes=[], select_mode='min', cols1=('chrom1', 'start1', 'end1'))
    loop = standard_assign(loop, H3K27ac, 'a2', 'H3K27ac', select_col='distance', includes=[], select_mode='min', cols1=('chrom2', 'start2', 'end2'))
    # end assign annoations here

    loop[['start1','end1','start2','end2']] = loopsecp
    # loop.to_csv(f'{loop_file}.attach', sep='\t', index=False)

    for fpkm_file in fpkm_files:
        fpkm_name = os.path.basename(os.path.dirname(fpkm_file))
        loop['a1_type'] = ['A' if not numpy.isnan(fpa) else 'I' if ta>0 else 'E' if ha>0 else 'O' for fpa, ta, ha in zip(loop[f'a1_TSS_FPKM_{fpkm_name}'], loop['TSS_per_a1'], loop['H3K27ac_per_a1'])]
        loop['a2_type'] = ['A' if not numpy.isnan(fpa) else 'I' if ta>0 else 'E' if ha>0 else 'O' for fpa, ta, ha in zip(loop[f'a2_TSS_FPKM_{fpkm_name}'], loop['TSS_per_a2'], loop['H3K27ac_per_a2'])]
        loop['a1_a2_type'] = ["_".join(sorted([a1t,a2t])) for a1t, a2t in zip(loop['a1_type'], loop['a2_type'])]
        loop['a1_a2_type_collapse'] = ['P_P' if tp in ['A_A','A_I','I_I'] else 'E_P' if tp in ['A_E','E_I'] else 'O_P' if tp in ['I_O','A_O'] else tp for tp in loop['a1_a2_type']]

        AI = numpy.concatenate([loop.loc[loop['a1_type']=='A','a1_type'].values, loop.loc[loop['a2_type']=='A','a2_type'].values, loop.loc[loop['a1_type']=='I','a1_type'].values, loop.loc[loop['a2_type']=='I','a2_type'].values])

        f, ax = matplotlib.pyplot.subplots(figsize=(10,10))
        seaborn.histplot(x=AI, stat='percent', ax=ax)
        for p in ax.patches:
            ax.text(p.get_x() + p.get_width()/2., p.get_height(), f'{p.get_height(): .2f}', ha="center", va="bottom")
        f.savefig(f'results/express_and_loop_strength/{os.path.basename(loop_file)}.AI.pdf')
        f.tight_layout()
        matplotlib.pyplot.close(f)

        for tycol in ['a1_a2_type', 'a1_a2_type_collapse']:
            f, ax = matplotlib.pyplot.subplots(figsize=(10,10))
            seaborn.histplot(data=loop, x=tycol, stat='percent', ax=ax)
            for p in ax.patches:
                ax.text(p.get_x() + p.get_width()/2., p.get_height(), f'{p.get_height(): .2f}', rotation=90, ha="center", va="bottom")
            ax.tick_params(axis='x', rotation=90)
            f.savefig(f'results/express_and_loop_strength/{os.path.basename(loop_file)}.{tycol}.pdf')
            f.tight_layout()
            matplotlib.pyplot.close(f)

        with open('results/express_and_loop_strength/znf143_loop.fit.summary', 'w') as fw:
            for types, name in zip([['A_A','A_I','I_I','A_E','E_I'],['A_A','A_I','I_I'],['A_E','E_I']],['pppe','pp','pe']):
                for met in ['lin', 'log']:
                    mask = loop['a1_a2_type'].isin(types)
                    score = loop.loc[mask,'score'].values
                    mean_FPKM = loop.loc[mask,f'FPKM_{fpkm_name}'].values
                    if met=='log':
                        score = numpy.log10(score)
                        mean_FPKM = numpy.log10(mean_FPKM)
                    regmask = (~numpy.isinf(score)) & (~numpy.isinf(mean_FPKM)) & (~numpy.isnan(score)) & (~numpy.isnan(mean_FPKM))
                    score , mean_FPKM = score[regmask], mean_FPKM[regmask]
                    slope, intercept, rvalue, pvalue, stderr = scatter_linearreg_cmap(score, mean_FPKM, f'results/express_and_loop_strength/znf143_loop.{name}.{met}.{fpkm_name}.pdf', alternative='greater')
                    fw.write(f'znf143_loop.{name}.{met}.pdf\t{slope}\t{intercept}\t{rvalue}\t{pvalue}\t{stderr}\n')

    # xs = ['mean_FPKM_WT','HiC_WT_norm']
    # ys = ['mean_FPKM_ZO','HiC_ZO_norm']
    # zs = ['score','score']
    # for x, y, z in zip(xs, ys, zs):
    #     f, ax = matplotlib.pyplot.subplots(figsize=(20,20))
    #     ax.scatter(loop[x], loop[y], c=loop[z], cmap='bwr', alpha=0.2, norm=matplotlib.colors.Normalize(vmin=0, vmax=numpy.percentile(loop[z],90)))
    #     ax.set(aspect=1)
    #     xlim = ax.get_xlim()
    #     ylim = ax.get_ylim()
    #     xylim = [min(xlim[0],ylim[0]),max(xlim[1],ylim[1])]
    #     ax.plot(xylim, xylim)
    #     ax.tick_params(axis='x', rotation=90)
    #     ax.set(xlim=xylim, ylim=xylim)
    #     f.savefig(f'results/express_and_loop_strength/{os.path.basename(loop_file)}.{x}.{y}.{z}.pdf')
    #     f.tight_layout()
    #     matplotlib.pyplot.close(f)



### tss based analysis
fpkm_filess = [['single_cut_RNA/cufflink/AID/genes.fpkm_tracking','single_cut_RNA/cufflink/AID-ZO/genes.fpkm_tracking','single_cut_RNA/cufflink/AID-CD/genes.fpkm_tracking','single_cut_RNA/cufflink/AID-ZOCD/genes.fpkm_tracking']] * 2
H3K4me3_bwfiless = [['???',"???"]] * 2
loop_files = ['zmloops/AID_loop_ZZ.bedpe', 'zmloops/AID-CD_loop_ZZ.bedpe']

for fpkm_files, H3K4me3_bwfiles, loop_file in zip(fpkm_filess, H3K4me3_bwfiless, loop_files):
    tss = get_tss('hg19.ncbiRefSeq.gtf.gz')
    tss = tss.drop_duplicates(subset=['start'], ignore_index=True)
    coding = bioframe.read_table('hg19_gene_id_name.bed', schema=['chrom', 'start', 'end', 'strand', 'name', 'gene_id', 'type'])
    coding = coding[coding['type']=='protein_coding'].drop_duplicates(subset=['gene_id'], ignore_index=True)
    tss = pandas.merge(tss, coding['gene_id'], on='gene_id', how='inner')
    tss = tss[[chr.find('_')==-1 for chr in tss['chrom']]]

    for fpkm_file in fpkm_files:
        fpkm_name = os.path.basename(os.path.dirname(fpkm_file))
        fpkm = pandas.read_csv(fpkm_file, sep='\t')[["gene_id",'FPKM']].groupby('gene_id')
        sumFPKM = fpkm['FPKM'].sum().values
        fpkm = fpkm.nth(0).reset_index()
        fpkm['FPKM'] = sumFPKM
        tss = pandas.merge(tss, fpkm[['gene_id', 'FPKM']], on='gene_id', how='left')
        tss = tss.rename(columns={'FPKM' : f'FPKM_{fpkm_name}'})
        tss[f'FPKM_{fpkm_name}'] = tss[f'FPKM_{fpkm_name}'].fillna(0.0)
        tss[f'log_FPKM_{fpkm_name}'] = numpy.log10(tss[f'FPKM_{fpkm_name}']+1)

    # flank = 2000
    # for H3K4me3_bwfile in H3K4me3_bwfiles:
    #     tss[H3K4me3_bwfile] = bbi.stackup(H3K4me3_bwfiles, tss['chrom'], tss['start']-flank, tss['end']+flank, bins=1, summary='mean').flatten()

    # if len(fpkm_files)==2 and len(H3K4me3_bwfiles)==2:
    #     H3K4me3_WmZ = tss[H3K4me3_bwfiles[1]]-tss[H3K4me3_bwfiles[0]]
    #     fpkm_name0, fpkm_name1 = os.path.basename(os.path.dirname(fpkm_files[0])), os.path.basename(os.path.dirname(fpkm_files[1]))
    #     FPKM_WmD = tss[f'FPKM_{fpkm_name1}']-tss[f'FPKM_{fpkm_name0}']
    #     iso = sklearn.ensemble.IsolationForest()
    #     yhat = iso.fit_predict(FPKM_WmD.values.reshape(-1,1))
    #     H3K4me3_WmZ = H3K4me3_WmZ[yhat!=-1]
    #     FPKM_WmD = FPKM_WmD[yhat!=-1]
    #     _ = scatter_linearreg_cmap(H3K4me3_WmZ, FPKM_WmD, f'results/express_and_loop_strength/FPKM_2_H3K4me3.pdf', alternative='greater')

    extend = 2000
    loop = bioframe.read_table(loop_file, schema=['chrom1','start1','end1','chrom2','start2','end2','score'])
    loop['start1'], loop['end1'] = (loop['start1'] + loop['end1'])//2-extend, (loop['start1'] + loop['end1'])//2+extend
    loop['start2'], loop['end2'] = (loop['start2'] + loop['end2'])//2-extend, (loop['start2'] + loop['end2'])//2+extend

    tss = standard_assign(tss, loop, '', 'a1', select_col='distance', includes=[], select_mode='min', cols2=('chrom1', 'start1', 'end1'))
    tss = standard_assign(tss, loop, '', 'a2', select_col='distance', includes=[], select_mode='min', cols2=('chrom2', 'start2', 'end2'))
    tss['anchor'] = ['ANC' if a1+a2>0 else 'NANC' for a1, a2 in zip(tss['a1_per_'], tss['a2_per_'])]

    for ycol in ['FPKM', 'log_FPKM']:
        ctc = [f"{ycol}_{os.path.basename(os.path.dirname(fpkm_file))}" for fpkm_file in fpkm_files]
        for anchor in ['ANC', 'NANC']:
            f, ax = matplotlib.pyplot.subplots(figsize=(20,20))
            
            seaborn.barplot(data=tss[tss['anchor']==anchor][ctc], ax=ax)
            ylim = list(ax.get_ylim())
            title = ""
            for fpkm_file_pair in itertools.combinations(fpkm_files, 2):
                fpkm_name0, fpkm_name1 = os.path.basename(os.path.dirname(fpkm_file_pair[0])), os.path.basename(os.path.dirname(fpkm_file_pair[1]))
                pvalue = scipy.stats.ttest_rel(tss.loc[tss['anchor']==anchor,f'FPKM_{fpkm_name0}'], tss.loc[tss['anchor']==anchor,f'FPKM_{fpkm_name1}'], axis=0, nan_policy='omit').pvalue
                title += f"{fpkm_name0}_{fpkm_name1}={pvalue:.5f}\n"
            ax.set(title=title)
            f.savefig(f'results/express_and_loop_strength/RNAseq.bar.{os.path.basename(loop_file)}.{ycol}.{anchor}.pdf')
            matplotlib.pyplot.close(f)

        f, ax = matplotlib.pyplot.subplots(figsize=(20,20))
        seaborn.barplot(data=tss[ctc], ax=ax)
        ylim = list(ax.get_ylim())
        title = ""
        for fpkm_file_pair in itertools.combinations(fpkm_files, 2):
            fpkm_name0, fpkm_name1 = os.path.basename(os.path.dirname(fpkm_file_pair[0])), os.path.basename(os.path.dirname(fpkm_file_pair[1]))
            pvalue = scipy.stats.ttest_rel(tss[f'FPKM_{fpkm_name0}'], tss[f'FPKM_{fpkm_name1}'], axis=0, nan_policy='omit').pvalue
            title += f"{fpkm_name0}_{fpkm_name1}={pvalue:.5f}\n"
        ax.set(title=title)
        f.savefig(f'results/express_and_loop_strength/RNAseq.bar..{os.path.basename(loop_file)}.{ycol}.pdf')
        matplotlib.pyplot.close(f)


### H3K4me3 based analysis
samples = ['WT', 'ZO', 'WT_rep1', 'WT_rep2', 'K2_rep1', 'K2_rep2', 'K45_rep1', 'K45_rep2']
counts = {}
H3K4me3s = {}
for sample in samples:
    H3K4me3s[sample] = bioframe.read_table(f'H3K4me3/peak/{sample}_peaks.narrowPeak', schema=['chrom','start','end','name','score','strand','signal','pvalue','qvalue'])
    counts[sample]=len(H3K4me3s[sample])

f, ax = matplotlib.pyplot.subplots(figsize=(10,10))
seaborn.barplot(data=pandas.DataFrame(counts, index=[0]), ax=ax)
for p in ax.patches:
    ax.text(p.get_x() + p.get_width()/2., p.get_height(), f'{int(p.get_height())}', ha="center", va="bottom")
WTK2p = scipy.stats.ttest_ind([counts['WT_rep1'],counts['WT_rep2']], [counts['K2_rep1'],counts['K2_rep2']]).pvalue
WTK45p = scipy.stats.ttest_ind([counts['WT_rep1'],counts['WT_rep2']], [counts['K45_rep1'],counts['K45_rep2']]).pvalue
ax.set(title=f"WTK2p={WTK2p:.5f}, WTK45p={WTK45p:.5f}")
f.savefig(f'results/express_and_loop_strength/H3K4me3.count.pdf')
matplotlib.pyplot.close(f)

inter_counts = {}
for spe in [['WT','ZO'],['ZO','WT']]:
    inter_counts[f'{spe[0]}_inter'] = sum(standard_assign(H3K4me3s[spe[0]], H3K4me3s[spe[1]], '', spe[1], select_col='distance', includes=[], select_mode='min')[f'{spe[1]}_per_']>0)
    inter_counts[f'{spe[0]}_only'] = counts[spe[0]] - inter_counts[f'{spe[0]}_inter']
f, ax = matplotlib.pyplot.subplots(figsize=(10,10))
seaborn.barplot(data=pandas.DataFrame(inter_counts, index=[0]), ax=ax)
for p in ax.patches:
    ax.text(p.get_x() + p.get_width()/2., p.get_height(), f'{int(p.get_height())}', ha="center", va="bottom")
f.savefig(f'results/express_and_loop_strength/H3K4me3.inter.count.pdf')
matplotlib.pyplot.close(f)

extend = 2000
loop_file = 'results/express_and_loop_strength/ZNF143_loop_ZZ.bedpe'
loop = bioframe.read_table(loop_file, schema=['chrom1','start1','end1','chrom2','start2','end2','score'])
loop['start1'], loop['end1'] = (loop['start1'] + loop['end1'])//2-extend, (loop['start1'] + loop['end1'])//2+extend
loop['start2'], loop['end2'] = (loop['start2'] + loop['end2'])//2-extend, (loop['start2'] + loop['end2'])//2+extend

for sample in samples:
    H3K4me3s[sample] = standard_assign(H3K4me3s[sample], loop, '', 'a1', select_col='distance', includes=[], select_mode='min', cols2=('chrom1', 'start1', 'end1'))
    H3K4me3s[sample] = standard_assign(H3K4me3s[sample], loop, '', 'a2', select_col='distance', includes=[], select_mode='min', cols2=('chrom2', 'start2', 'end2'))
    H3K4me3s[sample]['anchor'] = ['ANC' if a1+a2>0 else 'NANC' for a1, a2 in zip(H3K4me3s[sample]['a1_per_'], H3K4me3s[sample]['a2_per_'])]

for anchor in ['ALL', 'NANC', 'ANC']:
    if anchor!='ALL':
        subH3K4me3s = {}
        for sample in samples:
            subH3K4me3s[sample] = H3K4me3s[sample][H3K4me3s[sample]['anchor']==anchor]
    else:
        subH3K4me3s = H3K4me3s
    sampleslong = []
    slong = {'signal' : [], 'pvalue' : [], 'qvalue' : []}
    for sample in samples:
        sampleslong += [sample]*len(subH3K4me3s[sample])
        for ycol in slong.keys():
            slong[ycol] += subH3K4me3s[sample][ycol].values.tolist()

    ycoldf = pandas.DataFrame({'sample' : sampleslong, 'signal' : slong['signal'], 'pvalue' : slong['pvalue'], 'qvalue' : slong['qvalue']})

    for ycol in slong.keys():
        f, ax = matplotlib.pyplot.subplots(figsize=(20,10))
        seaborn.violinplot(data=ycoldf, x='sample', y=ycol, ax=ax)
        WTZOp=scipy.stats.ttest_ind(subH3K4me3s['WT'][ycol].values.tolist(), subH3K4me3s['ZO'][ycol].values.tolist()).pvalue
        WTK2p=scipy.stats.ttest_ind(subH3K4me3s['WT_rep1'][ycol].values.tolist()+subH3K4me3s['WT_rep2'][ycol].values.tolist(), subH3K4me3s['K2_rep1'][ycol].values.tolist()+subH3K4me3s['K2_rep2'][ycol].values.tolist()).pvalue
        WTK45p=scipy.stats.ttest_ind(subH3K4me3s['WT_rep1'][ycol].values.tolist()+subH3K4me3s['WT_rep2'][ycol].values.tolist(), subH3K4me3s['K45_rep1'][ycol].values.tolist()+subH3K4me3s['K45_rep2'][ycol].values.tolist()).pvalue
        WTKp=scipy.stats.ttest_ind(subH3K4me3s['WT_rep1'][ycol].values.tolist()+subH3K4me3s['WT_rep2'][ycol].values.tolist(), subH3K4me3s['K2_rep1'][ycol].values.tolist()+subH3K4me3s['K2_rep2'][ycol].values.tolist()+subH3K4me3s['K45_rep1'][ycol].values.tolist()+subH3K4me3s['K45_rep2'][ycol].values.tolist()).pvalue
        ax.set(title=f"WTZOp={WTZOp}, WTK2p={WTK2p}, WTK45p={WTK45p}, WTKp={WTKp}")
        f.savefig(f'results/express_and_loop_strength/H3K4me3.{anchor}.{ycol}.pdf')
        matplotlib.pyplot.close(f)

### Delta_ZNF143 Delta_CTCF
loop = bioframe.read_table('AID_HiC_loop_domain/AID_loop.bedpe', schema=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2'])
samples = ['AID','AID-ZO','AID-CD','AID-ZOCD']
clrs, stacks = {}, {}
for sample in samples:
    clr = load_rename_add_normVec_cov_tot(f'hic_matrices/{sample}.allValidPairs.mcool', 10_000)
    cvd = cooltools.expected_cis(clr=clr, view_df=hg19_arms, ignore_diags=0, nproc=4)
    loop[f"HiC_{sample}"] = cooltools.pileup(clr, loop, view_df=hg19_arms, expected_df=cvd, flank=0, nproc=4).flatten()

CTCF = bioframe.read_table('ChIP_Nexus_results/WT_CTCF_ChIP-nexus_peaks.narrowPeak', schema='bed')
ZNF143 = bioframe.read_table('ChIP_Nexus_results/WT_ZNF143_ChIP-nexus_peaks.narrowPeak', schema='bed')
loop = standard_assign(loop, ZNF143, 'a1', 'ZNF143', select_col='distance', includes=[], select_mode='min', cols1=('chrom1', 'start1', 'end1'))
loop = standard_assign(loop, ZNF143, 'a2', 'ZNF143', select_col='distance', includes=[], select_mode='min', cols1=('chrom2', 'start2', 'end2'))
loop = standard_assign(loop, CTCF, 'a1', 'CTCF', select_col='distance', includes=[], select_mode='min', cols1=('chrom1', 'start1', 'end1'))
loop = standard_assign(loop, CTCF, 'a2', 'CTCF', select_col='distance', includes=[], select_mode='min', cols1=('chrom2', 'start2', 'end2'))

loop['SBS'] = ['ESES' if za1>0 and za2>0 else 'NSNS' if za1==0 and za2==0 else 'ESNS' for za1, za2 in zip(loop['ZNF143_per_a1'], loop['ZNF143_per_a2'])]
loop['CBS'] = ['ECEC' if ca1>0 and ca2>0 else 'NCNC' if ca1==0 and ca2==0 else 'ECNC' for ca1, ca2 in zip(loop['CTCF_per_a1'], loop['CTCF_per_a2'])]
loop['SBS_CBS'] = loop['SBS'] + '_' + loop['CBS']
loop['ALL'] = 'ALL'

f, ax = matplotlib.pyplot.subplots(figsize=(20,20))
seaborn.histplot(data=loop, x='SBS_CBS', ax=ax)
for p in ax.patches:
    ax.text(p.get_x() + p.get_width()/2., p.get_height(), f'{int(p.get_height())}', ha="center", va="bottom")
f.savefig('results/AID_ZOCD/AID_ZOCD.loop.count.pdf')
matplotlib.pyplot.close(f)

for colcla in ['ALL', 'SBS_CBS', 'CBS']:
    for type in loop[colcla].unique():
        data = loop[loop[colcla]==type]
        for pairxy in itertools.combinations(['HiC_AID','HiC_AID-CD','HiC_AID-ZOCD'], 2):
            f, ax = matplotlib.pyplot.subplots(figsize=(20,20))
            seaborn.scatterplot(data=data, x=pairxy[0], y=pairxy[1], hue='SBS_CBS', alpha=0.2, s=20, label='loop types', ax=ax)
            ax.set(aspect=1)
            xlim, ylim = ax.get_xlim(), ax.get_ylim()
            xylim = [min(xlim[0],ylim[0]),max(xlim[1],ylim[1])]
            ax.plot(xylim, xylim)
            ax.tick_params(axis='x', rotation=90)
            ax.set(xlim=xylim, ylim=xylim)
            f.savefig(f'results/AID_ZOCD/AID_ZOCD.loop.express.{pairxy[0]}.{pairxy[1]}.{type}.pdf')
            f.tight_layout()
            matplotlib.pyplot.close(f)

        f, ax = matplotlib.pyplot.subplots(figsize=(20,20))
        seaborn.boxplot(data=data[[f"HiC_{sample}" for sample in samples]], ax=ax)
        title = ""
        for pairxy in itertools.combinations([f"HiC_{sample}" for sample in samples], 2):
            pvalue = scipy.stats.ttest_rel(data[pairxy[0]], data[pairxy[1]], axis=0, nan_policy='omit', alternative='greater').pvalue
            title += f"{pairxy[0]}_{pairxy[1]}={pvalue:.5f}, "
        ax.set(title=title)
        f.savefig(f'results/AID_ZOCD/AID_ZOCD.loop.express.box.{type}.pdf')
        f.tight_layout()
        matplotlib.pyplot.close(f)

        f, ax = matplotlib.pyplot.subplots(figsize=(20,20))
        seaborn.barplot(data=data[[f"HiC_{sample}" for sample in samples]], ax=ax)
        for p in ax.patches:
            ax.text(p.get_x() + p.get_width()/2., p.get_height(), f'{p.get_height():.5f}', ha="center", va="bottom")
        title = ""
        for pairxy in itertools.combinations(['HiC_AID','HiC_AID-CD','HiC_AID-ZOCD'], 2):
            pvalue = scipy.stats.ttest_rel(data[pairxy[0]], data[pairxy[1]], axis=0, nan_policy='omit', alternative='greater').pvalue
            title += f"{pairxy[0]}_{pairxy[1]}={pvalue:.5f}, "
        ax.set(title=title)
        f.savefig(f'results/AID_ZOCD/AID_ZOCD.loop.express.bar.{type}.pdf')
        f.tight_layout()
        matplotlib.pyplot.close(f)

for colcla in ['ALL', 'SBS_CBS', 'CBS']:
    for type in loop[colcla].unique():
        data = loop[loop[colcla]==type]
        for pairxy in itertools.combinations([f"HiC_{sample}" for sample in samples], 2):
            f, ax = matplotlib.pyplot.subplots(figsize=(20,20))
            seaborn.scatterplot(data=data, x=pairxy[0], y=pairxy[1], hue='SBS_CBS', alpha=0.2, s=20, label='loop types', ax=ax)
            ax.set(aspect=1)
            xlim, ylim = ax.get_xlim(), ax.get_ylim()
            xylim = [min(xlim[0],ylim[0]),max(xlim[1],ylim[1])]
            ax.plot(xylim, xylim)
            ax.tick_params(axis='x', rotation=90)
            ax.set(xlim=xylim, ylim=xylim)
            f.savefig(f'results/AID_ZOCD/AID_ZOCD.loop.express.{pairxy[0]}.{pairxy[1]}.{type}.pdf')
            f.tight_layout()
            matplotlib.pyplot.close(f)

### AID hichip loop

with open('hic_matrices/allValidPairs.ln', 'rb') as f:
    lns = pickle.load(f) 
dige = pandas.read_csv('MboI.bed', sep='\t')
CTCF = bioframe.read_table('ChIP_Nexus_results/AID_CTCF_ChIP-nexus_peaks.narrowPeak', schema='bed')
ZNF143 = bioframe.read_table('ChIP_Nexus_results/AID_ZNF143_ChIP-nexus_peaks.narrowPeak', schema='bed')
inter = True
block_size = 40000
nbdiges = [12]

loop_groups = [['results/express_and_loop_strength/AID_loop_ZZ.bedpe','results/express_and_loop_strength/AID_loop_ZO.bedpe','results/express_and_loop_strength/AID_loop_OO.bedpe'], ['results/express_and_loop_strength/AID-CD_loop_ZZ.bedpe','results/express_and_loop_strength/AID-CD_loop_ZO.bedpe','results/express_and_loop_strength/AID-CD_loop_OO.bedpe']]

for out_name, loop_group in zip(['AID','AID-CD'], loop_groups):
    # get loop
    loops_all = pandas.concat([bioframe.read_table(file, schema=['chrom1','start1','end1','chrom2','start2','end2','score']) for file in loop_group]).reset_index(drop=True)
    # attach CTCF ZNF143
    loops_all = standard_assign(loops_all, ZNF143, 'a1', 'ZNF143', select_col='distance', includes=[], select_mode='min', cols1=('chrom1', 'start1', 'end1'))
    loops_all = standard_assign(loops_all, ZNF143, 'a2', 'ZNF143', select_col='distance', includes=[], select_mode='min', cols1=('chrom2', 'start2', 'end2'))
    loops_all = standard_assign(loops_all, CTCF, 'a1', 'CTCF', select_col='distance', includes=[], select_mode='min', cols1=('chrom1', 'start1', 'end1'))
    loops_all = standard_assign(loops_all, CTCF, 'a2', 'CTCF', select_col='distance', includes=[], select_mode='min', cols1=('chrom2', 'start2', 'end2'))
    # # filter SBS only loops
    # boolESES = (loops_all['ZNF143_per_a1']>0) & (loops_all['ZNF143_per_a2']>0)
    # boolNCNC = (loops_all['CTCF_per_a1']==0) & (loops_all['CTCF_per_a2']==0)
    # loops_all = loops_all[boolESES & boolNCNC].reset_index(drop=True)
    # align to emzyme fragment
    loops_all['mid1s'] = (loops_all['start1'] + loops_all['end1'])//2
    loops_all['mid1e'] = loops_all['mid1s'] +1
    loops_all['mid2s'] = (loops_all['start2'] + loops_all['end2'])//2
    loops_all['mid2e'] = loops_all['mid2s'] +1
    loops_all['bin1_id'] =  bioframe.closest(loops_all, dige, return_input=False, return_distance=False, return_index=True, cols1=['chrom1','mid1s','mid1e'])['index_'].astype('int64')
    loops_all['bin2_id'] =  bioframe.closest(loops_all, dige, return_input=False, return_distance=False, return_index=True, cols1=['chrom2','mid2s','mid2e'])['index_'].astype('int64')
    for nbdige in nbdiges:
        # count
        MboI_cool_files = ['hic_matrices_MboI/WT-rep1.MboI.cool', 'hic_matrices_MboI/WT-rep2.MboI.cool', 'hic_matrices_MboI/ZO1-rep1.MboI.cool', 'hic_matrices_MboI/ZO1-rep2.MboI.cool', 'hic_matrices_MboI/ZO2-rep1.MboI.cool', 'hic_matrices_MboI/ZO2-rep2.MboI.cool', 'hic_matrices_MboI/AID-rep1.MboI.cool', 'hic_matrices_MboI/AID-rep2.MboI.cool', 'hic_matrices_MboI/AID-CD-rep1.MboI.cool', 'hic_matrices_MboI/AID-CD-rep2.MboI.cool', 'hic_matrices_MboI/AID-ZO-rep1.MboI.cool', 'hic_matrices_MboI/AID-ZO-rep2.MboI.cool', 'hic_matrices_MboI/AID-ZOCD-rep1.MboI.cool', 'hic_matrices_MboI/AID-ZOCD-rep2.MboI.cool']
        jobpacks = joblib.Parallel(n_jobs=len(MboI_cool_files)//2)(joblib.delayed(count_loop_valid2)(file, loops_all, nbdige, len(dige), block_size) for file in MboI_cool_files)
        for i, file in enumerate(MboI_cool_files):
            loops_all[os.path.basename(file)] = jobpacks[i][0]
            # loops_all[f"{os.path.basename(file)}.a1cov"] = jobpacks[i][1]
            # loops_all[f"{os.path.basename(file)}.a2cov"] = jobpacks[i][2]
        for file in MboI_cool_files:
            loops_all[f"{os.path.basename(file)}.norm"] = loops_all[os.path.basename(file)]/lns[os.path.basename(file)]*1e6
        loops_all.to_csv(f'results/playwith/{out_name}.{nbdige}.loops.count', sep='\t', index=False)


### TADs analysis
tadfile1 = '/media/ljw/f3b85364-e45d-4166-8db8-1cca425f188e1/HiC/output/TAD/WT_10000_iced.matrix_Domain_10Kb.bed'
tadfile2 = '/media/ljw/f3b85364-e45d-4166-8db8-1cca425f188e1/HiC/output/TAD/ZO1_10000_iced.matrix_Domain_10Kb.bed'
df_counts, df_pairs = ljwlib.tads.compare_tads(tadfile1, tadfile2)
# histcount of tads clusters
f = seaborn.displot(x=df_counts['num1'].astype(str)+'_'+df_counts['num2'].astype(str))
for p in f.ax.patches:
    f.ax.text(p.get_x() + p.get_width()/2., p.get_height(), f'{int(p.get_height())}', rotation=90, ha="center", va="bottom") 
f.ax.tick_params(axis='x', rotation=90)
f.savefig('results/tads/pair_num_counts.pdf')
# 1to1 tad-pair shift
f = seaborn.JointGrid(x=df_pairs['start'][1::2].values-df_pairs['start'][::2].values, y=df_pairs['end'][1::2].values-df_pairs['end'][::2].values, marginal_ticks=True)
f.plot_joint(seaborn.scatterplot, alpha=0.2)
f.plot_marginals(seaborn.histplot, kde=True)
f.savefig('results/tads/pair_shifts.pdf')


### renbin call tads matlab
for file in files:
    ljwlib.tads.renbin_call_tads_matlab(file, binsize=10000, DIwindows=2000000, boundsizemin=3, boundthres=0.99)

### renbin call tads domaincaller # pomegranate has been updated, but the min_cov in Normal is not implemented
for file in files:
    args = ['--uri', f'{file}::resolutions/10000', '--output', f'domaincaller_tads/{os.path.basename(file)}.tads', '--DI-output', f'domaincaller_tads/{os.path.basename(file)}.DIs', '--window-size', f'{2000000}', '--weight-col', 'weight', '--cpu-core', f'{12}', '--probs', f'{0.99}', '--minsize', f'{3}', '--logFile', f'domaincaller_tads/{os.path.basename(file)}.log']
    args, _ = domaincaller.renbin.getargs(args)
    domaincaller.renbin.call_tads(args)

### hictool call tads
exclude = ['chrM', 'chrY']
output_path = "hictool_tads"
os.makedirs(os.path.join(output_path,"tmp"), exist_ok=True)
for file in files:
    clr = load_rename_add_normVec_cov_tot(file, 10_000)
    bins = clr.bins()[:][['chrom','start','end']]
    bins = bins[[chrom not in exclude for chrom in bins['chrom']]]
    DIs, hmm_states, domains = [], [], []
    for chr in clr.chromnames:
        if chr not in exclude:
            HiCtool.scripts.HiCtool_TAD_analysis_lib.hictools_call_tads(action='full_tad_analysis', contact_matrix_or_DI=numpy.nan_to_num(clr.matrix().fetch(chr), nan=0.0), chromSizes_path='HiCtool/scripts/chromSizes', isGlobal=False, tab_sep=True, chrs=[f"{chromindices[chr]+1}"], species='hg19', data_type='normalized', full_chromosome=True, coord=None, input_file_hmm=None, plot_legend=True, plot_grid=True, bin_size=clr.binsize, zero_threshold=0.1, output_path=os.path.join(output_path,"tmp"))
            DIs.extend(numpy.loadtxt(os.path.join(output_path,"tmp",f"HiCtool_chr{chromindices[chr]+1}_DI.txt")).tolist())
            hmm_states.extend(numpy.loadtxt(os.path.join(output_path,"tmp",f"HiCtool_chr{chromindices[chr]+1}_hmm_states.txt")).tolist())
            domains.append(bioframe.read_table(os.path.join(output_path,"tmp",f"HiCtool_chr{chromindices[chr]+1}_topological_domains.txt"), schema=['start', 'end']))
            domains[-1].insert(0, 'chrom', chr)
    pandas.concat([bins, pandas.DataFrame({'DI' : DIs})], axis=1).to_csv(os.path.join(output_path,f"{os.path.basename(file)}.DI"), index=False, header=False, sep='\t', na_rep='nan')
    pandas.concat([bins, pandas.DataFrame({'hmm_state' : hmm_states})], axis=1).to_csv(os.path.join(output_path,f"{os.path.basename(file)}.hmm"), index=False, header=False, sep='\t', na_rep='nan')
    pandas.concat(domains).to_csv(os.path.join(output_path,f"{os.path.basename(file)}.dm"), index=False, header=False, sep='\t', na_rep='nan')
shutil.rmtree(os.path.join(output_path,"tmp"))

### tadtree call tads # extremely slow
for file in files:
    clr = load_rename_add_normVec_cov_tot(file, 10_000)
    TADtree.TADtreelib.use_TADtree(clr, 'TADtree_tads', M=10, gamma=500)

### scrabe tracks from ucsc, and compare it with ZNF143_CBSNA_SBSNA
ZNF143_CBSNA_SBSNA = pandas.read_csv('scrab_ucsc/ZNF143_CBSNA_SBSNA.CTCF.bed', sep='\t')
ljwlib.scrab_ucsc.scrab_ucsc('scrab_ucsc/ZNF143_CBSNA_SBSNA.log', 'scrab_ucsc/ZNF143_CBSNA_SBSNA.err', functools.partial(ljwlib.scrab_ucsc.colocalize_count, ZNF143_CBSNA_SBSNA))