import matplotlib.pyplot, numpy, pandas, cooler, cooltools, bioframe, mpl_toolkits.axes_grid1, os, matplotlib.ticker, matplotlib.colors, itertools, skimage.filters, bbi, matplotlib.patches, scipy.stats, pymemesuite, pymemesuite.common, pymemesuite.fimo, subprocess, scipy, networkx

chromsizes = {'chr1' : 249250621, 'chr2' : 243199373, 'chr3' : 198022430, 'chr4' : 191154276, 'chr5' : 180915260, 'chr6' : 171115067, 'chr7' : 159138663, 'chr8' : 146364022, 'chr9' : 141213431, 'chr10' : 135534747, 'chr11' : 135006516, 'chr12' : 133851895, 'chr13' : 115169878, 'chr14' : 107349540, 'chr15' : 102531392, 'chr16' : 90354753, 'chr17' : 81195210, 'chr18' : 78077248, 'chr19' : 59128983, 'chr20' : 63025520, 'chr21' : 48129895, 'chr22' : 51304566, 'chrM' : 16571, 'chrX' : 155270560, 'chrY' : 59373566}
chromindices = {'chr1' : 0, 'chr2' : 1, 'chr3' : 2, 'chr4' : 3, 'chr5' : 4, 'chr6' : 5, 'chr7' : 6, 'chr8' : 7, 'chr9' : 8, 'chr10' : 9, 'chr11' : 10, 'chr12' : 11, 'chr13' : 12, 'chr14' : 13, 'chr15' : 14, 'chr16' : 15, 'chr17' : 16, 'chr18' : 17, 'chr19' : 18, 'chr20' : 19, 'chr21' : 20, 'chr22' : 21, 'chrM' : 22, 'chrX' : 23, 'chrY' : 24}
hg19_arms = bioframe.make_chromarms(bioframe.fetch_chromsizes('hg19'), bioframe.fetch_centromeres('hg19'))
hg19_arms = bioframe.sort_bedframe(hg19_arms, view_df = bioframe.make_viewframe(chromsizes), df_view_col='chrom')
genome = bioframe.load_fasta('/home/ljw/hg19_with_random/hg19.fa')

### convert validpairs to .hic and .mcool
def vp2hicmcool(pair_file, hicpro2juicebox='/home/ljw/bin/HiC-Pro-master/bin/utils/hicpro2juicebox.sh', chromsizesfile='/home/ljw/hg19_with_bowtie2_index/hg19.chrom.sizes'):
    subprocess.check_output(f'bash {hicpro2juicebox} -i {pair_file} -o {os.path.dirname(pair_file)} -g <(grep -v "chr.*_" {chromsizesfile} | sort -k1,1V) -j juicer_tools_1.22.01.jar', shell=True)
    subprocess.check_output(f'hic2cool convert {pair_file}.hic {pair_file}.mcool', shell=True)

### load and pre-process of cooler
def load_rename_add_normVec_cov_tot(file, resolution, minDis=0, force=False):
    chr_dict = {'X' : 'chrX', 'Y' : 'chrY', 'MT' : 'chrM', 'M' : 'chrM'}
    for i in range(1, 23):
        chr_dict[f'{i}'] = f'chr{i}'

    clr = cooler.Cooler(f'{file}::/resolutions/{resolution}')
    cooler.rename_chroms(clr, chr_dict)
    if force and ('weight' in clr.bins().columns):
        with clr.open('r+') as f:
            del f['bins']['weight']
    if 'weight' not in clr.bins().columns:
        ctmp = cooler.Cooler(f'{file}::{cooler.fileops.list_coolers(file)[-1]}')
        btmp = ctmp.bins()[:]
        ptmp = ctmp.pixels()[:]
        totalreads = sum(ptmp.loc[(ptmp.bin2_id-ptmp.bin1_id>=minDis) & (btmp.chrom[ptmp.bin1_id].reset_index(drop=True)==btmp.chrom[ptmp.bin2_id].reset_index(drop=True)), 'count'])
        with clr.open('r+') as f:
            f["bins"].create_dataset("weight", data=1./f["bins/SCALE"][:]/(totalreads**0.5)*1e4, compression="gzip", compression_opts=6)
    if 'cov_cis_raw' not in clr.bins().columns:
        cooltools.coverage(clr, ignore_diags=2, store=True)
    return clr

### draw bins heatmap
def getBinRange(bins):
    starts = []
    ends = []
    chri = bins.chrom.iloc[0]
    end = 0
    for bin in bins.itertuples():
        if bin.chrom!=chri:
            chri = bin.chrom
            end += lastend
        starts.append(bin.start+end)
        ends.append(bin.end+end)
        lastend = bin.end
    end += lastend
    return (bins.start.iloc[0], end), numpy.array(starts), numpy.array(ends)

def binheat(matrix, bins, annot=False, *args, **kwargs):
    f, ax = matplotlib.pyplot.subplots(figsize=(12,12))
    binrange, _, _ = getBinRange(bins)
    im = ax.matshow(matrix, extent=(binrange[0], binrange[1], binrange[1], binrange[0]), *args, **kwargs)
    if annot:
        pos = numpy.linspace(binrange[0], binrange[1], matrix.shape[0]+1)
        pos = (pos[:-1]+pos[1:])/2
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                ax.text(pos[i], pos[j], f'{matrix[i][j]:.2f}', ha='center', va='center')
    ax.set(xticks=binrange, xticklabels=binrange, yticks=binrange, yticklabels=binrange)
    divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    try:
        f.colorbar(im, cax=cax)
        return f, ax, divider
    except Exception as err:
        matplotlib.pyplot.close(f)
        raise err
    


### print region hic
def print_hic_map(clr, output, region=None, ylim=None, *args, **kwargs):
    os.makedirs(output, exist_ok=True)

    if not region:
        matrix = clr.matrix()[:]
        bins = clr.bins()[:]
    else:
        matrix = clr.matrix().fetch(region)
        bins = clr.bins().fetch(region)
    try:    
        f, ax, divider = binheat(matrix, bins, *args, **kwargs)
    except Exception as err:
        print(repr(err))
        return
    
    _, starts, ends = getBinRange(bins)
    cax = divider.append_axes('bottom', size='15%', pad=0.1, sharex=ax)
    cax.plot((starts+ends)//2, bins['cov_cis_raw']/bins['cov_tot_raw'])
    cax.set(ylabel='cis/tot', ylim=ylim, yticks=ylim, yticklabels=ylim)

    f.tight_layout()
    f.savefig(os.path.join(output, f'{os.path.basename(clr.filename)}.{clr.binsize}.{region}.pdf'))
    matplotlib.pyplot.close(f)

### get nan corrcoef
def nancorrcoef(matrix):
    matrix = matrix.copy()
    std = numpy.nanstd(matrix,axis=1)
    nanmat = (~numpy.isnan(matrix)).astype('float64')
    matrix = numpy.nan_to_num(matrix, nan=0.0)
    matrix = matrix - numpy.mean(matrix, axis=1)
    pearson = numpy.matmul(matrix.T, matrix)
    pearson /= numpy.matmul(nanmat, nanmat)
    pearson /= std
    pearson /= numpy.reshape(std,[-1,1])
    return pearson


### calculate compartment for region
def compartmentalization(clr, genome, chr, output, *args, **kwargs):
    os.makedirs(output, exist_ok=True)

    outpearson = os.path.join(output, f'{os.path.basename(clr.filename)}.{chr}.{clr.binsize}.pearson')
    if not os.path.isfile(outpearson):
        subprocess.check_output(f'java -jar juicer_tools_1.22.01.jar pearsons -p SCALE {clr.filename[:-5]}hic {chr} BP {clr.binsize} {outpearson}', shell=True)
    pearson = numpy.loadtxt(outpearson)
    bins = clr.bins().fetch(chr)
    try:
        f, ax, divider = binheat(pearson, bins, *args, **kwargs)
    except Exception as err:
        print(repr(err))
        return
    _, starts, ends = getBinRange(bins)

    gc_cov = bioframe.frac_gc(clr.bins().fetch(chr)[['chrom', 'start', 'end']], genome)
    cis_eigs = cooltools.eigs_cis(clr, gc_cov, view_df=bioframe.make_viewframe((chr, 0, clr.chromsizes[chr])), n_eigs=1, ignore_diags=0)
    cax = divider.append_axes('bottom', size='15%', pad=0.1, sharex=ax)
    cax.bar(starts, numpy.maximum(cis_eigs[1].E1,0.0), width=ends-starts, align='edge', color='r', label='A')
    cax.bar(starts, numpy.minimum(cis_eigs[1].E1,0.0), width=ends-starts, align='edge', color='b', label='B')
    cax.set(ylabel='E1')

    outeigen = os.path.join(output, f'{os.path.basename(clr.filename)}.{chr}.{clr.binsize}.eigen')
    if not os.path.isfile(outeigen):
        subprocess.check_output(f'java -jar juicer_tools_1.22.01.jar eigenvector -p SCALE {clr.filename[:-5]}hic {chr} BP {clr.binsize} {outeigen}', shell=True)
    jE1 = numpy.loadtxt(outeigen)
    if scipy.stats.pearsonr(numpy.nan_to_num(cis_eigs[1].E1, nan=0.0), numpy.nan_to_num(jE1, nan=0.0)).statistic<0:
        jE1 = -jE1
    cax = divider.append_axes('bottom', size='15%', pad=0.1, sharex=ax)
    cax.bar(starts, numpy.maximum(jE1,0.0), width=ends-starts, align='edge', color='r', label='A')
    cax.bar(starts, numpy.minimum(jE1,0.0), width=ends-starts, align='edge', color='b', label='B')
    cax.set(ylabel='jE1')

    cax = divider.append_axes('bottom', size='15%', pad=0.1, sharex=ax)
    cax.bar(starts, gc_cov.GC, width=ends-starts, align='edge', color='r', label='GC')
    cax.set(ylabel='GC')

    # for i in numpy.where(numpy.diff((cis_eigs[1].E1>0).astype(int)))[0]:
    #     ax.axhline(ends[i], region[1], region[2], c='k', lw=0.25)
    #     ax.axvline(ends[i], region[1], region[2], c='k', lw=0.25)

    f.tight_layout()
    f.savefig(os.path.join(output, f'{os.path.basename(clr.filename)}.{chr}.{clr.binsize}.compartment.pdf'))
    matplotlib.pyplot.close(f)


### Saddleplots
def saddle_plots(clr, genome, chrom, output, n_bins=38, QL=0.025, QU=0.975, vmin=None, vmax=None, *args, **kwargs):
    os.makedirs(output, exist_ok=True)

    gc_cov = bioframe.frac_gc(clr.bins().fetch(chrom)[['chrom', 'start', 'end']], genome)
    try:
        cis_eigs = cooltools.eigs_cis(clr, gc_cov, view_df=bioframe.make_viewframe((chrom,0,clr.chromsizes[chrom])), n_eigs=1, ignore_diags=0)
    except Exception as err:
        print(repr(err))
        return
    track = bioframe.select(cis_eigs[1][['chrom','start','end','E1']], (chrom,0,clr.chromsizes[chrom]))
    digitized_track, _ = cooltools.digitize(track, n_bins, vrange=None, qrange=(QL, QU))
    groupmean = track[track.columns[3]].groupby(digitized_track[digitized_track.columns[3]]).mean()
    regionFrame = pandas.DataFrame([{'chrom' : chrom, 'start' : 0, 'end' : clr.chromsizes[chrom], 'name' : chrom}])
    cvd = cooltools.expected_cis(clr=clr, view_df=regionFrame, ignore_diags=2, nproc=12)
    try:
        interaction_sum, interaction_count =  cooltools.saddle(clr, cvd, track, 'cis', n_bins=n_bins, qrange=(QL, QU), view_df=regionFrame)
    except Exception as err:
        print(repr(err))
        return
    binedges = numpy.linspace(QL, QU, n_bins + 1)

    groupmean = groupmean[2:-1]
    interaction_sum = interaction_sum[1:-1, 1:-1]
    interaction_count = interaction_count[1:-1, 1:-1]
    binedges = binedges[:-1]
    
    ### saddle
    f, ax = matplotlib.pyplot.subplots(figsize=(12, 10))
    im = ax.matshow(interaction_sum/interaction_count, extent=(QL, QU, QU, QL), *args, **kwargs)
    ax.xaxis.set(visible=False)
    ax.yaxis.set(visible=False)
    divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    f.colorbar(im, cax=cax, label='interaction mean')
    cax = divider.append_axes("top", size="20%", pad=0.25, sharex=ax)
    cax.bar(binedges, height=groupmean, width=(QU-QL)/n_bins, align="edge")
    cax.set(ylim=[vmin,vmax], yticks=[vmin,vmax], yticklabels=[vmin,vmax])
    cax.xaxis.set(ticks_position='top')
    cax = divider.append_axes("left", size="20%", pad=0.25, sharey=ax)
    cax.barh(binedges, height=(QU-QL)/n_bins, width=groupmean, align="edge")
    cax.set(xlim=[vmin,vmax], xticks=[vmin,vmax], xticklabels=[vmin,vmax])
    f.tight_layout()
    f.savefig(os.path.join(output, f'{os.path.basename(clr.filename)}.{clr.binsize}.{chrom}.saddle.pdf'))
    matplotlib.pyplot.close(f)

    ### saddle_strength
    stren = cooltools.api.saddle.saddle_strength(interaction_sum, interaction_count)
    f, ax = matplotlib.pyplot.subplots(figsize=(12, 10))
    ax.plot(numpy.arange(n_bins+1), numpy.append(stren, stren[-1]), drawstyle='steps-post')
    ax.axhline(0, c='black', ls='--', lw=1)
    ax.set(xlim=[0, n_bins], xlabel='extent', ylabel='(AA + BB) / (AB + BA)', title='saddle strength profile')
    f.tight_layout()
    f.savefig(os.path.join(output, f'{os.path.basename(clr.filename)}.{clr.binsize}.{chrom}.saddle_strength.pdf'))
    matplotlib.pyplot.close(f)


### print the P(s) curve
def draw_Ps(clr, hg19_arms, chr, output):
    os.makedirs(output, exist_ok=True)
  
    cvd = cooltools.expected_cis(clr=clr, view_df=hg19_arms[hg19_arms['chrom']==chr].reset_index(drop=True), intra_only=False, ignore_diags=0, nproc=12)
    regions = hg19_arms['name'].loc[hg19_arms['chrom'] == chr].values
    f, axs = matplotlib.pyplot.subplots(figsize=(7,9), nrows=2, gridspec_kw={'height_ratios':[2,1]}, sharex=True)
    for regionlist in [[regions[0]], [regions[1]], regions]:
        regBool = (cvd['region1']==regionlist[0]) & (cvd['region2']==regionlist[-1])
        dist = cvd['dist'].loc[regBool] * clr.binsize
        avg = cvd['balanced.avg.smoothed'].loc[regBool]
        der = numpy.gradient(numpy.log(avg), numpy.log(dist))
        
        axs[0].loglog(dist ,avg, alpha=0.5, label=regionlist)
        axs[1].semilogx(dist, der, alpha=0.5, label=regionlist)
        axs[0].set(ylabel='IC contact frequency')
        axs[1].set(xlabel='distance', ylabel='slope')
        
    axs[0].legend()
    axs[1].legend()
    f.tight_layout()
    f.savefig(os.path.join(output, f'{os.path.basename(clr.filename)}.{clr.binsize}.{chr}.P(s).pdf'))
    matplotlib.pyplot.close(f)


### Averaging interaction frequencies in blocks
def average_trans_frequency(clr, output, *args, **kwargs):
    os.makedirs(output, exist_ok=True)
    
    ac = cooltools.expected_trans(clr, view_df = None)
    acp = (ac.pivot_table(values="balanced.avg", index="region1", columns="region2", observed=True).reindex(index=clr.chromnames, columns=clr.chromnames))
    f, ax = matplotlib.pyplot.subplots(figsize=(7,6))
    im = ax.matshow(acp, extent=(0, len(clr.chromnames), len(clr.chromnames), 0), *args, **kwargs)
    f.colorbar(im, ax=ax, label='average normalized contact frequency')
    ax.set(xticks=[x+0.5 for x in range(len(clr.chromnames))], yticks=[x+0.5 for x in range(len(clr.chromnames))], xticklabels=clr.chromnames, yticklabels=clr.chromnames, title="average inter-interactions")
    for pos in range(1,len(clr.chromnames)):
        ax.axhline(pos, 0, len(clr.chromnames), color='k')
        ax.axvline(pos, 0, len(clr.chromnames), color='k')
    ax.tick_params(axis='x', rotation=90)
    f.tight_layout()
    f.savefig(os.path.join(output, f'{os.path.basename(clr.filename)}.{clr.binsize}.avg_trans.pdf'))
    matplotlib.pyplot.close(f)


### insulation boundaries
def tads_strength(clr, window, output):
    os.makedirs(output, exist_ok=True)

    insulation_table = cooltools.insulation(clr, [window], ignore_diags=2)

    f, ax = matplotlib.pyplot.subplots(sharex=True, figsize=(8,8), constrained_layout=True)
    ax.hist(insulation_table[f'boundary_strength_{window}'], bins=10**numpy.linspace(-4,1,200), histtype='step', lw=1, label=f'Window {window}')
    threshold_li = skimage.filters.threshold_li(insulation_table[f'boundary_strength_{window}'].dropna().values)
    threshold_otsu = skimage.filters.threshold_otsu(insulation_table[f'boundary_strength_{window}'].dropna().values)
    n_boundaries_li = (insulation_table[f'boundary_strength_{window}'].dropna()>=threshold_li).sum()
    n_boundaries_otsu = (insulation_table[f'boundary_strength_{window}'].dropna()>=threshold_otsu).sum()
    ax.axvline(threshold_li, c='green', label=f'{n_boundaries_li} boundaries (Li)')
    ax.axvline(threshold_otsu, c='magenta', label=f'{n_boundaries_otsu} boundaries (Otsu)')
    ax.set(xscale='log', ylabel='# boundaries')
    ax.legend()
    ax.set(xlabel='Boundary strength')

    f.tight_layout()
    f.savefig(os.path.join(output, f'{os.path.basename(clr.filename)}.{clr.binsize}.{window}.tad_strength.pdf'))
    matplotlib.pyplot.close(f)
    return insulation_table, threshold_li, threshold_otsu

def boundary_signal(clr, window, bwfile, insulation_table, threshold_li, threshold_otsu, output, flank=1000, pileup_flank=50000, nbins=100):
    os.makedirs(output, exist_ok=True)

    # ctcf strength
    ctcf_chip_signal = bbi.stackup(bwfile, insulation_table.chrom, insulation_table.start-flank, insulation_table.end+flank, bins=1).flatten()
    f, ax = matplotlib.pyplot.subplots()
    ax.loglog(insulation_table[f'boundary_strength_{window}'], ctcf_chip_signal, 'o', markersize=1, alpha=1)
    ax.set(xlabel='Boundary strength', ylabel='CTCF enrichment over input')
    ax.axvline(threshold_li, ls='--', color='green', label='boundaries (Li)')
    ax.axvline(threshold_otsu, ls='--', color='magenta', label='boundaries (Otsu)')
    ax.legend()

    f.tight_layout()
    f.savefig(os.path.join(output, f'{os.path.basename(clr.filename)}.{clr.binsize}.{window}.{os.path.basename(bwfile)}.ctcf_strength.pdf'))
    matplotlib.pyplot.close(f)

    # pileup ctcf
    top_boundaries = insulation_table[insulation_table[f'boundary_strength_{window}']>=threshold_otsu]
    mid = (top_boundaries.start + top_boundaries.end)//2
    stackup = bbi.stackup(bwfile, top_boundaries.chrom, mid-pileup_flank, mid+pileup_flank, bins=nbins)

    f, ax = matplotlib.pyplot.subplots(figsize=[7,5])
    ax.plot(numpy.linspace(0, 2*pileup_flank, nbins+1)[:-1], numpy.nanmean(stackup, axis=0))
    ax.set(xlabel='Distance from boundary, kbp', ylabel='CTCF ChIP-Seq mean fold change over input')
    ax.xaxis.set(major_formatter=matplotlib.ticker.EngFormatter('b'))

    f.tight_layout()
    f.savefig(os.path.join(output, f'{os.path.basename(clr.filename)}.{clr.binsize}.{window}.{os.path.basename(bwfile)}.ctcf_pileup.pdf'))
    matplotlib.pyplot.close(f)

def draw_region_tads(clr, region, window, insulation_table, output, *args, **kwargs):
    os.makedirs(output, exist_ok=True)

    start_pos_vector = numpy.arange(region[1], region[2]+1, clr.binsize)
    matrix_c = clr.matrix().fetch(region)
    matrix_a = numpy.dot(numpy.array([(i[1], i[0]) for i in itertools.product(start_pos_vector[::-1], start_pos_vector)]), numpy.array([[1, 0.5], [-1, 0.5]]))
    x = matrix_a[:, 1].reshape(len(start_pos_vector), len(start_pos_vector))
    y = matrix_a[:, 0].reshape(len(start_pos_vector), len(start_pos_vector))

    f, ax = matplotlib.pyplot.subplots(figsize=(18, 6))
    im = ax.pcolormesh(x, y, numpy.flipud(matrix_c), *args, **kwargs)
    im.set_rasterized(True)
    ax.set(aspect=0.5, xlim=[region[1], region[2]], ylim=[0, 10*window])
    ax.yaxis.set(major_formatter=matplotlib.ticker.EngFormatter('b'))
    ax.xaxis.set_visible(False)

    divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)
    cax = divider.append_axes("right", size="1%", pad=0.1)
    f.colorbar(im, cax=cax)

    insul_region = bioframe.select(insulation_table, region)
    # weak_boundaries = insul_region[~insul_region[f'is_boundary_{window}']]
    strong_boundaries = insul_region[insul_region[f'is_boundary_{window}']]

    cax = divider.append_axes("bottom", size="50%", pad=0., sharex=ax)
    # cax.scatter(weak_boundaries[['start', 'end']].mean(axis=1), weak_boundaries[f'log2_insulation_score_{window}'], label=f'Weak {window} boundaries')
    strong_poses = strong_boundaries[['start', 'end']].mean(axis=1).values
    strong_scores = strong_boundaries[f'log2_insulation_score_{window}'].values
    cax.scatter(strong_poses, strong_scores, label=f'Strong {window} boundaries')
    for sp, sc in zip(strong_poses, strong_scores):
        ax.text(sp, sc, f'{sc:.2f}', ha="center", va="bottom") 
    cax.plot(insul_region[['start', 'end']].mean(axis=1), insul_region[f'log2_insulation_score_{window}'], label=f'Window {window} bp')
    cax.set(ylabel=f'log2_insulation_score_{window}')
    cax.xaxis.set(major_formatter=matplotlib.ticker.EngFormatter('b'))
    cax.legend()

    mean_str = strong_boundaries[['start', 'end']].mean(axis=1).values
    for i in range(len(mean_str)):
        now = mean_str[i]
        if i>0:
            pre = mean_str[i-1]
        else:
            pre = region[1]
        if i<len(mean_str)-1:
            post = mean_str[i+1]
        else:
            post = region[2]
        ax.plot([now, (now+pre)/2], [0, now-pre], c='k', lw=1)
        ax.plot([now, (now+post)/2], [0, post-now], c='k', lw=1)
    
    f.tight_layout()
    f.savefig(os.path.join(output, f'{os.path.basename(clr.filename)}.{clr.binsize}.{region}.{window}.tads.pdf'))
    matplotlib.pyplot.close(f)


### draw_kernels
def draw_kernel(clr, output, cmap='viridis'):
    os.makedirs(output, exist_ok=True)

    kernels = cooltools.api.dotfinder.recommend_kernels(clr.binsize)
    f, axs = matplotlib.pyplot.subplots(ncols=len(kernels), figsize=(3*len(kernels),2.5))
    for ax, (ktype, kernel) in zip(axs, kernels.items()):
        ax.imshow(kernel[::-1,::-1], alpha=0.85, cmap=cmap, interpolation='nearest')
        # draw a square around the target pixel:
        x0 = kernel.shape[0] // 2 - 0.5
        y0 = kernel.shape[1] // 2 - 0.5
        rect = matplotlib.patches.Rectangle((x0, y0), 1, 1, lw=1, ec='r', fc='r')
        ax.add_patch(rect)
        # clean axis:
        ax.xaxis.set(visible=False)
        ax.yaxis.set(visible=False)
        ax.set(title=f'{ktype} kernel'.format(ktype))
        # add a checkerboard to highlight pixels:
        checkerboard = numpy.add.outer(range(kernel.shape[0]), range(kernel.shape[1])) % 2
        ax.imshow(checkerboard, cmap='gray', interpolation='nearest', alpha=0.3)
    f.tight_layout()
    f.savefig(os.path.join(output, f'{os.path.basename(clr.filename)}.{clr.binsize}.kernels.pdf'))
    matplotlib.pyplot.close(f)  


### loop-calling
def rectangles_around_dots(dots_df, region, loc="upper", lw=1, ec="cyan", fc="none"):
    # select dots from the region:
    df_reg = bioframe.select(bioframe.select(dots_df, region, cols=("chrom1","start1","end1")), region, cols=("chrom2","start2","end2"))
    # draw rectangular "boxes" around pixels called as dots in the "region":
    for s1, s2, e1, e2 in df_reg[["start1", "start2", "end1", "end2"]].itertuples(index=False):
        if loc == "upper":
            yield matplotlib.patches.Rectangle((s2, s1), e2 - s2, e1 - s1, lw=lw, ec=ec, fc=fc)
        elif loc == "lower":
            yield matplotlib.patches.Rectangle((s1, s2), e1 - s1, e2 - s2, lw=lw, ec=ec, fc=fc)
        else:
            raise ValueError("loc has to be upper or lower")

def loop_calling(clr, dots_df, region, output, *args, **kwargs):
    os.makedirs(output, exist_ok=True)

    bins = clr.bins().fetch(region)
    matrix = clr.matrix().fetch(region)
    for diag in [-1,0,1]:
        matrix = cooltools.lib.numutils.fill_diag(matrix, numpy.nan, i=diag)
    try:    
        f, ax, _ = binheat(matrix, bins, *args, **kwargs)
    except Exception as err:
        print(repr(err))
        return
    
    dot_df_reg = bioframe.select(bioframe.select(dots_df, region, cols=("chrom1","start1","end1")), region, cols=("chrom2","start2","end2"))
    for s1, s2, e1, e2 in dot_df_reg[["start1", "start2", "end1", "end2"]].itertuples(index=False):
        # upper box
        box = matplotlib.patches.Rectangle((s2, s1), e2 - s2, e1 - s1, lw=1, ec='blue', fc='none')
        # lower box
        # box = matplotlib.patches.Rectangle((s1, s2), e1 - s1, e2 - s2, lw=1, ec='blue', fc='none')
        ax.add_patch(box)

    f.tight_layout()
    f.savefig(os.path.join(output, f'{os.path.basename(clr.filename)}.{clr.binsize}.{region}.loop.pdf'))
    matplotlib.pyplot.close(f)

### pile up
def pile_up_clips(clr, sites, cvd, hg19_arms, name, output, flank=100_000, *args, **kwargs): # norm=matplotlib.colors.LogNorm()
    os.makedirs(output, exist_ok=True)

    stack = cooltools.pileup(clr, sites, view_df=hg19_arms, expected_df=cvd, flank=flank, nproc=12)
    mtx = numpy.nanmean(stack, axis=2)

    f, ax = matplotlib.pyplot.subplots()
    im = ax.matshow(mtx, extent=(-flank, flank, flank, -flank), *args, **kwargs)
    ax.text(0, -flank, f'{mtx[mtx.shape[0]//2,mtx.shape[1]//2]}', ha='center', va='bottom')
    ax.text(0, flank, f'{stack.shape[2]}', ha='center', va='top')
    f.colorbar(im, ax=ax, label ='normalized mean obs/exp')
    ax.set(xticks=[-flank, flank], xticklabels=[-flank, flank], xlabel='relative position', yticks=[-flank, flank], yticklabels=[-flank, flank], ylabel='relative position')
    f.savefig(os.path.join(output, f'{os.path.basename(clr.filename)}.{clr.binsize}.{name}.pileup.pdf'))
    matplotlib.pyplot.close(f)

    return stack



### standard peak assign operation, assign peak2 (columns in includes) to peak1, if several peak2 overlap peak1, select that with the min/max value in select_col
def standard_assign(peak1, peak2, peak1_name, peak2_name, select_col='distance', includes=['name', 'strand', 'distance'], select_mode='min', cols1=('chrom', 'start', 'end'), cols2=('chrom', 'start', 'end')):
    peak1 = peak1.copy()
    peak1_peak2 = bioframe.overlap(peak1, peak2, suffixes=('_1','_2'), cols1=cols1, cols2=cols2, return_index=True)
    peak1_peak2['distance_2'] = (peak1_peak2[f'{cols2[1]}_2'] + peak1_peak2[f'{cols2[2]}_2'] - peak1_peak2[f'{cols1[1]}_1'] - peak1_peak2[f'{cols1[2]}_1']) // 2
    peak1[f'{peak2_name}/{peak1_name}'] = peak1_peak2.groupby(["index_1"])["index_2"].count().values
    if select_col in ['pvalue']:
        peak1_peak2['tmp'] = peak1_peak2[f'{select_col}_2'].fillna(-1.0)
    elif select_col in ['distance']:
        peak1_peak2['tmp'] = peak1_peak2[f'{select_col}_2'].fillna(numpy.iinfo(peak1_peak2.dtypes[f'{select_col}_2'].type).max).abs()
    if select_mode=='min':
        idbest = peak1_peak2.groupby(["index_1"])["tmp"].idxmin().values
    elif select_mode=='max':
        idbest = peak1_peak2.groupby(["index_1"])["tmp"].idxmax().values
    for include in includes:
        peak1[(f'{peak1_name}_{peak2_name}_{include}').lstrip('_')] = peak1_peak2.loc[idbest,f'{include}_2'].values
    for chr in chromsizes:
        mask1chr = peak1[cols1[0]]==chr
        mask2chr = peak2[cols2[0]]==chr
        mid1 = pandas.DataFrame({'mid' : (peak1.loc[mask1chr,cols1[1]].values + peak1.loc[mask1chr,cols1[2]].values) // 2}).sort_values(by='mid')
        mid2 = pandas.DataFrame({'mid' : (peak2.loc[mask2chr,cols2[1]].values + peak2.loc[mask2chr,cols2[2]].values) // 2}).sort_values(by='mid')
        mid2['pos'] = mid2['mid']
        tmp = pandas.merge_asof(mid1, mid2, on='mid', direction='nearest')
        peak1.loc[mask1chr, f'{peak1_name}_{peak2_name}_closest'.lstrip('_')] = tmp.pos.values - tmp.mid.values

    return peak1

### get tss information from gtf
def get_tss(gtffile):
    genes = bioframe.read_table(gtffile, schema='gtf')
    tss_start = [row.start if row.strand=='+' else row.end for row in genes[genes.feature=='transcript'].itertuples()]
    gene_id = [row.attributes.split(';')[0].split()[1].strip('"') for row in genes[genes.feature=='transcript'].itertuples()]
    transcript_id = [row.attributes.split(';')[1].split()[1].strip('"') for row in genes[genes.feature=='transcript'].itertuples()]
    tss.loc[tss.chrom=='chrMT', 'chrom'] = 'chrM'
    tss = pandas.DataFrame({'chrom' : genes[genes.feature=='transcript'].chrom, 'start' : tss_start, 'end' : numpy.array(tss_start)+1, 'name' : transcript_id, 'score' : '.', 'strand' : genes[genes.feature=='transcript'].strand, 'gene_id' : gene_id})

    return tss


### draw deeptools-like plots
def deeptools_like(mat, extent, outfile):
    mat = mat[numpy.argsort(-mat.sum(axis=1)),:]
    f, ax = matplotlib.pyplot.subplots(figsize=(5,10))
    im = ax.matshow(mat, norm='linear', extent=extent, cmap='Reds', aspect='auto')
    ax.set(yticks=[0, mat.shape[0]])
    divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)
    cax = divider.append_axes("right", size="10%", pad=0.2)
    f.colorbar(im, cax=cax, label='normalized counts')
    cax = divider.append_axes("top", size="25%", pad=0.2)
    edges = numpy.linspace(extent[0], extent[1], mat.shape[1]+1)
    cax.plot((edges[:-1]+edges[1:])/2, numpy.nanmean(mat, axis=0))
    f.tight_layout() 
    f.savefig(outfile)

### get overlapping loops, must specify names for loops
def get_overlap_loops(loops1, loops2):
    loops1 = loops1.copy()
    loops2 = loops2.copy()
    G = networkx.Graph()
    G.add_nodes_from(range(len(loops1)+len(loops2)))
    loops1['over'] = 'no'
    loops2['over'] = 'no'

    anchor1_overlap = bioframe.overlap(loops1, loops2, suffixes=('_1','_2'), cols1=('chrom1', 'start1', 'end1'), cols2=('chrom1', 'start1', 'end1'), return_index=True)
    anchor2_overlap = bioframe.overlap(loops1, loops2, suffixes=('_1','_2'), cols1=('chrom2', 'start2', 'end2'), cols2=('chrom2', 'start2', 'end2'), return_index=True)
    anchor1_overlap_group = anchor1_overlap.groupby(['index_1'])
    anchor2_overlap_group = anchor2_overlap.groupby(['index_1'])
    anchor1_mask = anchor1_overlap_group['index_2'].count().values > 0
    anchor2_mask = anchor2_overlap_group['index_2'].count().values > 0
    for i in range(len(loops1)):
        if anchor1_mask[i] or anchor2_mask[i]:
            loops1.loc[i,'over'] = 'half'
            anchor1_index_2 = anchor1_overlap_group.get_group(i)['index_2'].values.astype('int64') if anchor1_mask[i] else []
            anchor2_index_2 = anchor2_overlap_group.get_group(i)['index_2'].values.astype('int64') if anchor2_mask[i] else []
            for j in itertools.chain(anchor1_index_2, anchor2_index_2):
                if loops2.loc[j,'over'] == 'no':
                    loops2.loc[j,'over'] = 'half'
            common_index_2 = numpy.intersect1d(anchor1_index_2, anchor2_index_2)
            if len(common_index_2) > 0:
                loops1.loc[i,'over'] = 'yes'
                for j in common_index_2:
                    loops2.loc[j,'over'] = 'yes'
                    G.add_edge(i, len(loops1)+j)            
    
    loops1['comp'] = -1
    loops2['comp'] = -1
    comps = networkx.connected_components(G)
    for cid, comp in enumerate(comps):
        comp = numpy.array(list(comp))
        loops1.loc[comp[comp<len(loops1)],'comp'] = cid
        loops2.loc[comp[comp>=len(loops1)]-len(loops1),'comp'] = cid
    return loops1, loops2

### merge neighbour bins
def merge_neighbour_bins(dige, nbdige):
    groups = []
    for _, group in dige.copy().groupby('chrom'):
        group['start'] = numpy.concatenate([numpy.full(nbdige,group['start'].values[0]), group['start'].values[numpy.arange(len(group)-nbdige)]])
        group['end'] = numpy.concatenate([group['end'].values[numpy.arange(nbdige,len(group))], numpy.full(nbdige,group['end'].values[-1])])
        groups.append(group)
    return pandas.concat(groups).reset_index(drop=True)


### align to peak
def align_to(peaks1, peaks2, cols1=('chrom', 'start', 'end'), cols2=('chrom', 'start', 'end')):
    peaks1 = peaks1.copy()
    peaks1 = standard_assign(peaks1, peaks2, '', 'align', select_col='distance', includes=cols2, select_mode='min', cols1=cols1, cols2=cols2)
    for col1, col2 in zip(cols1, cols2):
        peaks1[col1] = peaks1[f'align_{col2}']
    return peaks1.drop(columns=['align/','align_closest']+[f"align_{col2}" for col2 in cols2])

### call motifs, save motif instances to motif_file.motifs
def call_motifs(peak_file, pwm_files):
    peaks = bioframe.read_table(peak_file, schema='bed').query(f'chrom in {list(chromsizes.keys())}').reset_index(drop=True)
    sequences = [pymemesuite.common.Sequence(genome[record.chrom].ff.fetch(record.chrom, record.start, record.end), name=record.name.encode()) for record in peaks.itertuples()]
    fimo = pymemesuite.fimo.FIMO(both_strands=True)

    motif_files = []
    for pwm_file in pwm_files:
        with pymemesuite.common.MotifFile(pwm_file) as motif_handle:
            motif = motif_handle.read()
            pattern = fimo.score_motif(motif, sequences, motif_handle.background)
        names, starts, stops, scores, strands, pvalues, qvalues = [], [], [], [], [], [], []
        for i, m in enumerate(pattern.matched_elements):
            names.append(m.source.accession.decode())
            starts.append(m.start)
            stops.append(m.stop)
            scores.append(m.score)
            strands.append(m.strand)
            pvalues.append(m.pvalue)
            qvalues.append(m.qvalue)

        pattern_df = pandas.DataFrame({'name': names, 'start_m': starts, 'end_m': stops, 'score_m' : scores, 'strand_m' : strands, 'pvalue' : pvalues, 'qvalue' : qvalues})
        result = pandas.merge(peaks, pattern_df, on='name')
        start_m = result.start + [min(a,b)-1 for a, b in zip(result.start_m, result.end_m)]
        result.end_m = result.start + [max(a,b) for a, b in zip(result.start_m, result.end_m)]
        result.start_m = start_m
        result['name_m'] = result.index.astype('str')
        result.name_m = [f'{os.path.basename(pwm_file)}_{id}' for id in result.name_m]
        
        motif_files.append(f'{pwm_file}.motifs')
        result[['chrom', 'start_m', 'end_m', 'name_m', 'score_m', 'strand_m', 'pvalue', 'qvalue']].to_csv(path_or_buf=motif_files[-1], sep='\t', header=False, index=False)
    
    return motif_files











