import bioframe, matplotlib.pyplot, seaborn, pandas, os, numpy, Bio.Seq, pysam, ljwlib.hic_module, subprocess, re, string, pysam, more_itertools, shutil, itertools, sys
sys.path.append('/home/ljw/wuqiang/lierlib')
import lierlib

##########################################################
# KMP HIC start
##########################################################

### extract ZNF143 region reference
region = ('chr11', 9_480_000, 9_552_000)
genome = bioframe.load_fasta('/home/ljw/hg19_with_random/hg19.fa')
with open('znf143.fa', 'w') as f:
    f.write(f'>znf143\n{genome[region[0]].ff.fetch(region[0], region[1], region[2])}\n')

### KMP map HIC to ZNF143 region
# subprocess.check_output(f'bash ljwlib/kmp.sh {region[0]} {region[1]} {region[2]}', shell=True)

### initialize cuts
exp_up = 9499952
exp_down = 9522663
K2_inv_up = [9499952, 9499952]
K2_inv_down = [9522663, 9522663]
K45_inv_up = [9499953, 9499952]
K45_inv_down = [9522663, 9522663]
K45_del1_up = 9499952
K45_del1_down = 9499959
K45_del2_up = 9522641
K45_del2_down = 9522729
AID_ZO_inv1_up = [9499952, 9499952]
AID_ZO_inv1_down = [9522664, 9522663]
AID_ZO_inv2_up = [9499952, 9499952]
AID_ZO_inv2_down = [9522670, 9522670]
AID_ZO_del_up = 9499952
AID_ZO_del_down = 9522663

### bin edges at the region of ZNF143
dige = pandas.read_csv('MboI.bed', sep='\t')
dige = bioframe.select(dige, region, cols=None).reset_index(drop=True)
bins = numpy.append(dige.start.values, dige.end.iloc[-1]) # binedges

### expected junctions
sg_K2 = [[K2_inv_up[0],K2_inv_down[0]], [K2_inv_up[1],K2_inv_down[1]]]
sg_K45 = [[K45_inv_up[0],K45_inv_down[0]], [K45_inv_up[1],K45_inv_down[1]], [K45_del1_up,K45_del1_down], [K45_del2_up,K45_del2_down]]
sg_AID = [[AID_ZO_inv1_up[0],AID_ZO_inv1_down[0]], [AID_ZO_inv1_up[1],AID_ZO_inv1_down[1]], [AID_ZO_inv2_up[0],AID_ZO_inv2_down[0]], [AID_ZO_inv2_up[1],AID_ZO_inv2_down[1]], [AID_ZO_del_up,AID_ZO_del_down]]

### heatmap based on KMP HIC
ext = 100
basenames = ['K2', 'K45', 'AID_ZO']
sg_all = [sg_K2, sg_K45, sg_AID]
for sg, basename in zip(sg_all, basenames):
    files = [f'check_genome_run/{basename}_rep{i}_R{j}.fq.gz.LPM.for' for i, j in [[1,1],[1,2],[2,1],[2,2]]]
    chimer = pandas.concat([bioframe.read_table(file, schema='bedpe')[['end1','start2','score']] for file in files]).reset_index(drop=True)
    Eboolend1 = numpy.full((len(chimer),), False)
    Eboolstart2 = numpy.full((len(chimer),), False)
    for edge in bins:
        Eboolend1 = numpy.logical_or(Eboolend1, numpy.isclose(chimer.end1.values, edge, rtol=0, atol=ext))
        Eboolstart2 = numpy.logical_or(Eboolstart2, numpy.isclose(chimer.start2.values, edge, rtol=0, atol=ext))
    Ebool = numpy.logical_or(Eboolend1, Eboolstart2) 
    binrange = [bins[0], bins[-1]]
    for mode in ['total','nonemzyme']:
        data = chimer if mode=='total' else chimer[~Ebool]
        f = seaborn.JointGrid(data=data, x="end1", y="start2", xlim=binrange, ylim=binrange, marginal_ticks=True)
        f.plot_joint(seaborn.histplot, weights=data.score, bins=list(bins), binrange=binrange)
        f.plot_marginals(seaborn.histplot, weights=data.score, bins=list(bins), binrange=binrange, kde=True)
        for sgp in sg:
            f.ax_joint.plot(sgp, sgp[::-1], c='k', ls='', mew=0.1, marker="+")
        f.ax_joint.set(xticks=binrange, yticks=binrange, xticklabels=binrange, yticklabels=binrange, xlabel='chr11', ylabel='chr11')
        f.ax_marg_x.xaxis.set(visible=False)
        f.ax_marg_y.yaxis.set(visible=False)
        f.savefig(f'check_genome/{basename}.{mode}.check_genome.pdf')

    for sgp in sg:
        rx = numpy.where(bins>sgp[0])[0][0]
        ry = numpy.where(bins>sgp[1])[0][0]
        binranges = [[bins[rx-1],bins[rx]],[bins[ry-1],bins[ry]]]
        
        binBool = (chimer.end1>=bins[rx-1]) & (chimer.end1<bins[rx]) & (chimer.start2>=bins[ry-1]) & (chimer.start2<bins[ry])
        f = seaborn.JointGrid(xlim=binranges[0], ylim=binranges[1], marginal_ticks=True)
        seaborn.histplot(data=chimer[binBool], x='end1', y='start2', weights='score', binrange=binranges, bins=300, ax=f.ax_joint)
        seaborn.histplot(data=chimer[binBool], x='end1', weights='score', binrange=binranges[0], bins=300, ax=f.ax_marg_x)
        seaborn.histplot(data=chimer[binBool], y='start2', weights='score', binrange=binranges[1], bins=300, ax=f.ax_marg_y)
        f.ax_joint.plot(sgp[0], sgp[1], c='k', ls='', mew=0.1, marker="+")
        f.ax_joint.set(xticks=binranges[0], yticks=binranges[1], xticklabels=binranges[0], yticklabels=binranges[1], xlabel='chr11', ylabel='chr11')
        f.ax_marg_x.xaxis.set(visible=False)
        f.ax_marg_y.yaxis.set(visible=False)
        f.savefig(f'check_genome/{basename}.{sgp}.updown.singlebin.check_genome.pdf')

        binBool = (chimer.end1>=bins[ry-1]) & (chimer.end1<bins[ry]) & (chimer.start2>=bins[rx-1]) & (chimer.start2<bins[rx])
        f = seaborn.JointGrid(xlim=binranges[1], ylim=binranges[0], marginal_ticks=True)
        seaborn.histplot(data=chimer[binBool], x='end1', y='start2', weights='score', binrange=binranges[::-1], bins=300, ax=f.ax_joint)
        seaborn.histplot(data=chimer[binBool], x='end1', weights='score', binrange=binranges[1], bins=300, ax=f.ax_marg_x)
        seaborn.histplot(data=chimer[binBool], y='start2', weights='score', binrange=binranges[0], bins=300, ax=f.ax_marg_y)
        f.ax_joint.plot(sgp[1], sgp[0], c='k', ls='', mew=0.1, marker="+")
        f.ax_joint.set(xticks=binranges[1], yticks=binranges[0], xticklabels=binranges[1], yticklabels=binranges[0], xlabel='chr11', ylabel='chr11')
        f.ax_marg_x.xaxis.set(visible=False)
        f.ax_marg_y.yaxis.set(visible=False)
        f.savefig(f'check_genome/{basename}.{sgp}.downup.singlebin.check_genome.pdf')
    
##########################################################
# KMP HIC end
##########################################################


##########################################################
# bowtie2 HIC start
##########################################################

### get junc fasta files
def get_junc_reference(chrom, cutup, cutdown, genome, extend=50, type='delete'):
    if type == 'inverse':
        refup = genome[chrom].ff.fetch(chrom, cutup[0]-extend, cutup[0]).upper() + Bio.Seq.Seq(genome[chrom].ff.fetch(chrom, cutdown[0]-extend, cutdown[0])).reverse_complement().__str__().upper()
        refdown = Bio.Seq.Seq(genome[chrom].ff.fetch(chrom, cutup[1], cutup[1]+extend)).reverse_complement().__str__().upper() + genome[chrom].ff.fetch(chrom, cutdown[1], cutdown[1]+extend).upper()
        return refup, refdown
    elif type == 'delete':
        ref = genome[chrom].ff.fetch(chrom, cutup-extend, cutup).upper() + genome[chrom].ff.fetch(chrom, cutdown, cutdown+extend).upper()
        return ref
    else:
        raise('type must be delete or inverse')


shutil.rmtree('check_genome_run/juncs_indices')
os.makedirs(f'check_genome_run/juncs_indices', exist_ok=True)
extend = 150
refname2cuts = {}
with open('check_genome_run/juncs_indices/AID_ZO_juncs.fa', 'w') as f:
    refup = get_junc_reference('chr11', exp_up, exp_up, genome, extend=extend, type='delete')
    f.write(f'>AID_ZO_exp_up\n{refup}\n')
    refdown = get_junc_reference('chr11', exp_down, exp_down, genome, extend=extend, type='delete')
    f.write(f'>AID_ZO_exp_down\n{refdown}\n')
    refname2cuts['AID_ZO_exp_up'] = [exp_up, exp_up]
    refname2cuts['AID_ZO_exp_down'] = [exp_down, exp_down]
    refup, refdown = get_junc_reference('chr11', AID_ZO_inv1_up, AID_ZO_inv1_down, genome, extend=extend, type='inverse')
    f.write(f'>AID_ZO_inv1_up\n{refup}\n')
    f.write(f'>AID_ZO_inv1_down\n{refdown}\n')
    refname2cuts['AID_ZO_inv1_up'] = [AID_ZO_inv1_up[0], AID_ZO_inv1_down[0]]
    refname2cuts['AID_ZO_inv1_down'] = [AID_ZO_inv1_up[1], AID_ZO_inv1_down[1]]
    refup, refdown = get_junc_reference('chr11', AID_ZO_inv2_up, AID_ZO_inv2_down, genome, extend=extend, type='inverse')
    f.write(f'>AID_ZO_inv2_up\n{refup}\n')
    f.write(f'>AID_ZO_inv2_down\n{refdown}\n')
    refname2cuts['AID_ZO_inv2_up'] = [AID_ZO_inv2_up[0], AID_ZO_inv2_down[0]]
    refname2cuts['AID_ZO_inv2_down'] = [AID_ZO_inv2_up[1], AID_ZO_inv2_down[1]]
    ref = get_junc_reference('chr11', AID_ZO_del_up, AID_ZO_del_down, genome, extend=extend, type='delete')
    f.write(f'>AID_ZO_del\n{ref}\n')
    refname2cuts['AID_ZO_del'] = [AID_ZO_del_up, AID_ZO_del_down]

with open('check_genome_run/juncs_indices/K2_juncs.fa', 'w') as f:
    refup = get_junc_reference('chr11', exp_up, exp_up, genome, extend=extend, type='delete')
    f.write(f'>K2_exp_up\n{refup}\n')
    refdown = get_junc_reference('chr11', exp_down, exp_down, genome, extend=extend, type='delete')
    f.write(f'>K2_exp_down\n{refdown}\n')
    refname2cuts['K2_exp_up'] = [exp_up, exp_up]
    refname2cuts['K2_exp_down'] = [exp_down, exp_down]
    refup, refdown = get_junc_reference('chr11', K2_inv_up, K2_inv_down, genome, extend=extend, type='inverse')
    f.write(f'>K2_inv_up\n{refup}\n')
    f.write(f'>K2_inv_down\n{refdown}\n')
    refname2cuts['K2_inv_up'] = [K2_inv_up[0], K2_inv_down[0]]
    refname2cuts['K2_inv_down'] = [K2_inv_up[1], K2_inv_down[1]]

with open('check_genome_run/juncs_indices/K45_juncs.fa', 'w') as f:
    refup = get_junc_reference('chr11', exp_up, exp_up, genome, extend=extend, type='delete')
    f.write(f'>K45_exp_up\n{refup}\n')
    refdown = get_junc_reference('chr11', exp_down, exp_down, genome, extend=extend, type='delete')
    f.write(f'>K45_exp_down\n{refdown}\n')
    refname2cuts['K45_exp_up'] = [exp_up, exp_up]
    refname2cuts['K45_exp_down'] = [exp_down, exp_down]
    refup, refdown = get_junc_reference('chr11', K45_inv_up, K45_inv_down, genome, extend=extend, type='inverse')
    f.write(f'>K45_inv_up\n{refup}\n')
    f.write(f'>K45_inv_down\n{refdown}\n')
    refname2cuts['K45_inv_up'] = [K45_inv_up[0], K45_inv_down[0]]
    refname2cuts['K45_inv_down'] = [K45_inv_up[1], K45_inv_down[1]]
    ref = get_junc_reference('chr11', K45_del1_up, K45_del1_down, genome, extend=extend, type='delete')
    f.write(f'>K45_del1\n{ref}\n')
    refname2cuts['K45_del1'] = [K45_del1_up, K45_del1_down]
    ref = get_junc_reference('chr11', K45_del2_up, K45_del2_down, genome, extend=extend, type='delete')
    f.write(f'>K45_del2\n{ref}\n')
    refname2cuts['K45_del2'] = [K45_del2_up, K45_del2_down]

### bowtie2 map
for basename in basenames:
	subprocess.check_output(f'bowtie2-build check_genome_run/juncs_indices/{basename}_juncs.fa check_genome_run/juncs_indices/{basename}_juncs', shell=True)
	subprocess.check_output(f'bowtie2 -p 12 -x check_genome_run/juncs_indices/{basename}_juncs -U hic_rawdata/{basename}_rep1_R1.fq.gz hic_rawdata/{basename}_rep1_R2.fq.gz hic_rawdata/{basename}_rep2_R1.fq.gz hic_rawdata/{basename}_rep2_R2.fq.gz | samtools view -@ 12 -b -o check_genome_run/{basename}_juncs.bam', shell=True)

### filter mapped reads, sort and index bam
for basename in basenames:
	subprocess.check_output(f'samtools view -@ 12 -h check_genome_run/{basename}_juncs.bam | ' + 'awk \'{if ($0~/^@/) print; else if ($2%8<4) print;}\'' + f' | samtools view -b | samtools sort -o check_genome_run/{basename}_juncs.map.bam', shell=True)
	subprocess.check_output(f'samtools index -b check_genome_run/{basename}_juncs.map.bam', shell=True)

### tview reads
for basename in basenames:
    ref = bioframe.load_fasta(f'check_genome_run/juncs_indices/{basename}_juncs.fa')
    for refname in ref.keys():
        refseq = ref[refname].ff.fetch(refname)
        try:
            subprocess.check_output(f'samtools tview -p {refname} -d T -w {len(refseq)+100} check_genome_run/{basename}_juncs.map.bam check_genome_run/juncs_indices/{basename}_juncs.fa > check_genome_run/{refname}.display', shell=True)
        except Exception as err:
            print(repr(err))

### show map counts
for basename in basenames:
    refs = []
    with pysam.AlignmentFile(f'check_genome_run/{basename}_juncs.map.bam') as sf:
        for read in sf.fetch():
            refs.append(read.reference_name)
    f = seaborn.displot(x=ref, kind='hist')
    f.ax.tick_params(axis='x', rotation=90)
    f.savefig(f'check_genome/{basename}.align.count.pdf')

### show reads heatmap
def split_map(read, cutup, cutdown, extend=150):
    if read.reference_name.find('inv')>=0:
        if read.reference_name.find('up')>=0:
            mode = 'upstream_inverse'
        else:
            mode = 'downstream_inverse'
    else:
        mode = 'delete'
    starts, ends, refs, strands = [], [], [], []
    if mode=='delete':
        if read.reference_start >= extend:
            starts.append(cutdown-extend+read.reference_start)
            ends.append(cutdown-extend+read.reference_end)
        elif read.reference_end <= extend:
            starts.append(cutup-extend+read.reference_start)
            ends.append(cutup-extend+read.reference_end)
        else:
            starts.extend([cutup-extend+read.reference_start, cutdown])
            ends.extend([cutup, cutdown-extend+read.reference_end])
        if read.is_reverse:
            strands.extend(['-']*len(starts))
        else:
            strands.extend(['+']*len(starts))
    elif mode=='upstream_inverse':
        if read.reference_start >= extend:
            starts.append(cutdown+extend-read.reference_end)
            ends.append(cutdown+extend+read.reference_start)
            if read.is_reverse:
                strands.append('+')
            else:
                strands.append('-')
        elif read.reference_end <= extend:
            starts.append(cutup-extend+read.reference_start)
            ends.append(cutup-extend+read.reference_end)
            if read.is_reverse:
                strands.append('-')
            else:
                strands.append('+')
        else:
            starts.extend([cutup-extend+read.reference_start, cutdown+extend-read.reference_end])
            ends.extend([cutup, cutdown])
            if read.is_reverse:
                strands.extend(['-','+'])
            else:
                strands.extend(['+','-'])
    elif mode=='downstream_inverse':
        if read.reference_start >= extend:
            starts.append(cutdown-extend+read.reference_start)
            ends.append(cutdown-extend+read.reference_end)
            if read.is_reverse:
                strands.append('-')
            else:
                strands.append('+')
        elif read.reference_end <= extend:
            starts.append(cutup+extend-read.reference_end)
            ends.append(cutup+extend-read.reference_start)
            if read.is_reverse:
                strands.append('+')
            else:
                strands.append('-')
        else:
            starts.extend([cutup, cutdown])
            ends.extend([cutup+extend-read.reference_start, cutdown+read.reference_end-extend])
            if read.is_reverse:
                strands.extend(['+','-'])
            else:
                strands.extend(['-','+'])
    else:
        raise('wrong type')
    refs.extend([read.reference_name] * len(starts))
    return refs, starts, ends, strands

for basename in basenames:
    with pysam.AlignmentFile(f'check_genome_run/{basename}_juncs.map.bam') as f:
        refss, startss, endss, strandss = [], [], [], []
        for read in f.fetch():
            refs, starts, ends, strands = split_map(read, refname2cuts[read.reference_name][0], refname2cuts[read.reference_name][1], extend=150)
            refss.extend(refs)
            startss.extend(starts)
            endss.extend(ends)
            strandss.extend(strands)
    df = pandas.DataFrame({'ref' : refss, 'start' : startss, 'end' : endss, 'count' : [1]*len(refss), 'strand' : strandss})
    for refname in df['ref'].unique():
        dfsubfor = df[(df['ref']==refname) & (df['strand']=='+')]
        dfsubrev = df[(df['ref']==refname) & (df['strand']=='-')]
        f, ax = matplotlib.pyplot.subplots()
        if len(dfsubfor)>0:
            [fx,fy] = lierlib.fill_coverage(dfsubfor[['start','end']], dfsubfor['count'], '+')
            ax.fill(fx, fy, fc='k', ec=None)
        if len(dfsubrev)>0:
            [fx,fy] = lierlib.fill_coverage(dfsubrev[['start','end']], dfsubrev['count'], '-')
            ax.fill(fx, fy, fc='k', ec=None)
        ylim = list(ax.get_ylim())
        for cut in refname2cuts[refname]:
            ax.plot([cut,cut], ylim, c='k', ls='--')
        ax.set(ylim=ylim, yticks=list(set(ylim+[0])), yticklabels=list(set(ylim+[0])))
        f.tight_layout()
        f.savefig(f'check_genome/{basename}.{refname}.align.display.pdf')



##########################################################
# bowtie2 HIC end
##########################################################



##########################################################
# RNA-seq hisat2 start
##########################################################

### construct novel splice
genes = bioframe.read_table('hg19.ncbiRefSeq.gtf.gz', schema='gtf')
genes143 = genes[[attr.find('ZNF143')>=0 for attr in genes.attributes]].reset_index(drop=True)[['feature', 'start', 'end', 'strand', 'attributes']]
genes143 = genes143[[fea in ['transcript', 'exon'] for fea in genes143.feature]].reset_index(drop=True)
genes143['sa'] = ['transcript' if row.feature=='transcript' else row.attributes.split(';')[2].split()[1].strip('"') for row in genes143.itertuples()]
exon5 = genes143.loc[(genes143.feature=='exon') & (genes143.sa=='5'), ['start', 'end', 'strand', 'sa']].iloc[-1]
exon12 = genes143.loc[(genes143.feature=='exon') & (genes143.sa=='12'), ['start', 'end', 'strand', 'sa']].iloc[-1]
for basename in basenames:
    with open(f'hg19_hisat2/{basename}.novel', 'w') as f:
        f.write(f'chr11\t{exon5.end-1}\t{exon12.start-1}\t+\n')

### build HFM, feasible
subprocess.check_output('hisat2-build -p 12 hg19_hisat2/hg19.fa hg19_hisat2/hg19', shell=True)

### map reads
for basename in basenames:
    for rep in ['rep1', 'rep2']:
        subprocess.check_output(f'hisat2 -p 12 --pen-noncansplice 0 --known-splicesite-infile hg19_hisat2/splice_sites --novel-splicesite-infile hg19_hisat2/${basename}.novel --novel-splicesite-outfile hg19_hisat2/${basename}_${rep}.novel.out -x hg19_hisat2/hg19 -1 RNA-seq/${basename}_RNA-seq_${rep}_R1.fq -2 RNA-seq/${basename}_RNA-seq_${rep}_R2.fq | samtools view -@ 12 -b -o check_genome_run/${basename}_RNA-seq_${rep}.bam', shell=True)

### sort and index
for basename in basenames:
    for rep in ['rep1', 'rep2']:
        subprocess.check_output(f'samtools sort -@ 12 check_genome_run/{basename}_RNA-seq_{rep}.bam > check_genome_run/{basename}_RNA-seq_{rep}.sort.bam', shell=True)
        subprocess.check_output(f'samtools index -b check_genome_run/{basename}_RNA-seq_{rep}.sort.bam', shell=True)

### tview
for basename in basenames:
    for rep in ['rep1', 'rep2']:
        subprocess.check_output(f'samtools tview -p chr11:9496000 -d T -w 34500 check_genome_run/{basename}_RNA-seq_{rep}.sort.bam hg19_hisat2/hg19.fa > check_genome_run/{basename}_RNA-seq_{rep}.sort.bam.display', shell=True)

### igv
# /home/ljw/IGV_Linux_2.16.0_WithJava/IGV_Linux_2.16.0/igv.sh
# 9496097-9530393
# 9496000-9530500

### get exons
genes = bioframe.read_table('hg19.ncbiRefSeq.gtf.gz', schema='gtf')
exons = genes[genes.feature=='exon']
exons.chrom[exons.chrom=='chrMT'] = 'chrM'
genes143 = genes[[attr.find('ZNF143')>=0 for attr in genes.attributes]].reset_index(drop=True)[['feature', 'start', 'end', 'strand', 'attributes']]
genes143 = genes143[[fea in ['transcript', 'exon'] for fea in genes143.feature]].reset_index(drop=True)
exon2attr = {}
for i in range(1, 17):
    if i < 3:
        exon2attr[f'exon{i}'] = genes143.attributes[[i,16+i,33+i]].values
    elif i == 3:
        exon2attr[f'exon{i}'] = genes143.attributes[[16+i,33+i]].values
    else:
        exon2attr[f'exon{i}'] = genes143.attributes[[i-1,16+i,33+i]].values
genesRPL = genes[[attr.find('RPL23AP65')>=0 for attr in genes.attributes]].reset_index(drop=True)[['feature', 'start', 'end', 'strand', 'attributes']]
genesRPL = genesRPL[[fea in ['transcript', 'exon'] for fea in genesRPL.feature]].reset_index(drop=True)
RPL = genesRPL.attributes[[1]].values

### get blocks in the range and count exon5 an exon12 neighbours
def get_skip_blocks(read):
    cigars = read.cigarstring.split('N')
    skips = []
    for i in range(len(cigars)-1):
        skips.append(int(re.findall(r'\d+$', cigars[i])[0]))
        cigars[i] = cigars[i].rstrip(string.digits)
    tread = pysam.libcalignedsegment.AlignedSegment()
    blocks = [[read.reference_start, 0]]
    for i in range(len(cigars)-1):
        tread.cigarstring = cigars[i]
        blocks[i][1] = blocks[i][0] + tread.reference_length
        blocks.append([blocks[i][1]+skips[i], 0])
    tread.cigarstring = cigars[-1]
    blocks[-1][1] = blocks[-1][0] + tread.reference_length
    return blocks

def set_up_total_blocks(file):
    with pysam.AlignmentFile(f'check_genome_run/{file}') as sf:
        chromss, startss, endss, namess, scoress, strandss, strandidss = [], [], [], [], [], [], []
        for read in sf.fetch(region[0], region[1], region[2]):
            if read.is_unmapped or read.is_secondary:
                continue
            blocks = get_skip_blocks(read)
            starts = [block[0] for block in blocks]
            ends = [block[1] for block in blocks]
            strandids = [id for id in range(len(blocks))]
            if read.is_forward:
                strand = '+'
            else:
                strand = '-'
                strandids.reverse()
            chromss.extend([read.reference_name] * len(blocks))
            startss.extend(starts)
            endss.extend(ends)
            namess.extend([read.query_name + strand] * len(blocks))
            scoress.extend([read.mapping_quality] * len(blocks))
            strandss.extend([strand] * len(blocks))
            strandidss.extend(strandids)
            
    blocks = pandas.DataFrame({'chrom' : chromss, 'start' : startss, 'end' : endss, 'name' : namess, 'score' : scoress, 'strand' : strandss, 'strandid' : strandidss})
    blocks = blocks.sort_values(by = ['name', 'strandid']).reset_index(drop=True)

    return blocks

# count splice
files = ['K2_RNA-seq_rep1.sort.bam', 'K2_RNA-seq_rep2.sort.bam', 'K45_RNA-seq_rep1.sort.bam', 'K45_RNA-seq_rep2.sort.bam', 'AID_ZO_RNA-seq_rep1.sort.bam', 'AID_ZO_RNA-seq_rep2.sort.bam']
juncss, exonNamess = [], []
for file in files:
    blocks = set_up_total_blocks(file)
    blocks = ljwlib.hic_module.standard_assign(blocks, exons, '', 'exon', select_col='distance', includes=['attributes'], select_mode='min')
    juncs, exonNames = [], []
    for _, group in blocks.groupby(['name']):
        for i in range(len(group)):
            for k in [5,12]:
                if group['exon_attributes'].iloc[i] in exon2attr[f'exon{k}']:
                    for shif, dire in zip([-1,+1], ['pre','next']):
                        exonNames.append(f'exon{k}_{dire}{group.strand.iloc[i]}')
                        if i+shif < len(group) and i+shif >= 0:
                            junc = f'{group.exon_attributes.iloc[i+shif]}{group.strand.iloc[i+shif]}'
                            if not junc:
                                junc = 'No'
                            if group['exon_attributes'].iloc[i+shif] in RPL:
                                junc = f'RPL{group.strand.iloc[i+shif]}'
                            else:
                                for j in range(1,17):
                                    if group['exon_attributes'].iloc[i+shif] in exon2attr[f'exon{j}']:
                                        junc = f'exon{j}{group.strand.iloc[i+shif]}'
                                        break
                        else:
                            junc = 'No'
                        juncs.append(junc)
    juncss.append(juncs)
    exonNamess.append(exonNames)

hue_order = set([junc for juncs in juncss for junc in juncs])
for i in range(len(files)):
    f = seaborn.displot(x=exonNamess[i], hue=juncss[i], hue_order=hue_order, kind='hist', multiple="stack")
    f.ax.tick_params(axis='x', rotation=90)
    f.savefig(f'check_genome/{files[i]}.exon5_exon12_juncs.count.pdf')

### extract raw RNA-seq reads mapped by hisat2 to ZNF143 region
for basename in basenames:
    for rep in ['rep1', 'rep2']:
        with pysam.AlignmentFile(f'check_genome_run/{basename}_RNA-seq_{rep}.sort.bam') as sf:
            query_names1 = set([read.query_name for read in sf.fetch(region[0], region[1], region[2]) if read.is_read1])
            query_names2 = set([read.query_name for read in sf.fetch(region[0], region[1], region[2]) if read.is_read2])
        for i in [1,2]:
            query_names = eval(f'query_names{i}')
            with open(f'RNA-seq/{basename}_RNA-seq_{rep}_R{i}.fq', 'r') as f:
                os.makedirs(f'check_genome_run/{basename}_RNA-seq_{rep}_R{i}', exist_ok=True)
                fw = open(f'check_genome_run/{basename}_RNA-seq_{rep}_R{i}/{basename}_RNA-seq_{rep}_R{i}.hisat2.znf143.fq', 'w')
                for query_name, seq, plus, qual in more_itertools.batched(f, 4):
                    query_name = query_name.split()[0]
                    if query_name[1:] in query_names:
                        fw.write(f'{query_name}\n{seq}{plus}{qual}')        
                fw.close()

### do chimeric map
path = '/media/ljw/f3b85364-e45d-4166-8db8-1cca425f188e1/zhangmo'
for basename in basenames:
    for rep in ['rep1', 'rep2']:
        for i in [1,2]:
            run = f'{basename}_RNA-seq_{rep}_R{i}'
            folder = f'{path}/check_genome_run/{run}'
            subprocess.check_output(f'align.sh {folder} 4 "" "" "" "" 10 "" 12', shell=True)

### transform chimeric to sam
for basename in basenames:
    for rep in ['rep1', 'rep2']:
        for i in [1,2]:
            run = f'{basename}_RNA-seq_{rep}_R{i}'
            folder = f'{path}/check_genome_run/{run}'
            subprocess.check_output(f'to_sam.sh {folder} {run}', shell=True)

### draw chimeric heatmap
thres = 50
for sg, basename in zip(sg_all, basenames):
    for rep in ['rep1', 'rep2']:
        for i in [1,2]:
            file = f'{basename}_RNA-seq_{rep}_R{i}/{basename}_RNA-seq_{rep}_R{i}.bam'
            blocks = set_up_total_blocks(file)
            chroms1, starts1, ends1, chroms2, starts2, ends2, names, scores  = [], [], [], [], [], [], [], []
            for _, group in blocks.groupby(['name']):
                for j in range(len(group)-1):
                    if group['end'].iloc[j]-group['start'].iloc[j]<thres or group['end'].iloc[j+1]-group['start'].iloc[j+1]<thres:
                        pass
                    chroms1.append(group['chrom'].iloc[j])
                    chroms2.append(group['chrom'].iloc[j+1])
                    names.append(group['name'].iloc[j] + '_' + group['name'].iloc[j+1])
                    scores.append(1)
                    if group['strand'].iloc[j]=='+':
                        starts1.append(group['start'].iloc[j])
                        ends1.append(group['end'].iloc[j])
                        starts2.append(group['start'].iloc[j+1])
                        ends2.append(group['end'].iloc[j+1])
                    else:
                        starts1.append(group['end'].iloc[j])
                        ends1.append(group['start'].iloc[j])
                        starts2.append(group['end'].iloc[j+1])
                        ends2.append(group['start'].iloc[j+1])
                    
            chipairs = pandas.DataFrame({'chrom1' : chroms1, 'start1' : starts1, 'end1' : ends1, 'chrom2' : chroms2, 'start2' : starts2, 'end2' : ends2, 'name' : names, 'score' : scores})
            data = chipairs[(chipairs['chrom1']=='chr11') & (chipairs['chrom2']=='chr11')]
            nbins = 300
            binrange = [region[1],region[2]]
            f = seaborn.JointGrid(data=data, x="end1", y="start2", xlim=binrange, ylim=binrange, marginal_ticks=True)
            f.plot_joint(seaborn.histplot, weights=data.score, bins=nbins, binrange=binrange)
            f.plot_marginals(seaborn.histplot, weights=data.score, bins=nbins, binrange=binrange, kde=True)
            for sgp in sg:
                f.ax_joint.plot(sgp, sgp[::-1], c='k', ls='', mew=0.1, marker="+")
            for j in itertools.chain(range(1,15), range(17,32), range(34,49)):
                f.ax_joint.plot(genes143['end'][j], genes143['start'][j+1], c='r', ls='', mew=0.1, marker="+")
                f.ax_joint.plot(genes143['start'][j+1], genes143['end'][j], c='r', ls='', mew=0.1, marker="+")
            f.ax_joint.plot(9496180, 9530166, c='r', ls='', mew=0.1, marker="+")
            f.ax_joint.plot(9530166, 9496180, c='r', ls='', mew=0.1, marker="+")
            f.ax_joint.set(xticks=binrange, yticks=binrange, xticklabels=binrange, yticklabels=binrange, xlabel='chr11', ylabel='chr11')
            f.ax_joint.tick_params(axis='x', rotation=90)
            f.ax_marg_x.xaxis.set(visible=False)
            f.ax_marg_y.yaxis.set(visible=False)
            f.savefig(f'check_genome/{basename}_{rep}_R{i}.RNA_chimeric.pdf')

##########################################################
# RNA-seq hisat2 end
##########################################################