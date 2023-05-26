import os, subprocess, pysam, shutil
try:
    import deeptools.bamCoverage
except ValueError:
    import matplotlib
    matplotlib.colormaps.unregister("rocket")
    matplotlib.colormaps.unregister("rocket_r")
    matplotlib.colormaps.unregister("mako")
    matplotlib.colormaps.unregister("mako_r")
    matplotlib.colormaps.unregister("vlag")
    matplotlib.colormaps.unregister("vlag_r")
    matplotlib.colormaps.unregister("icefire")
    matplotlib.colormaps.unregister("icefire_r")
    import deeptools.bamCoverage

### run hicpro
def run_HiC_Pro(filepairs, inpath, outpath, config):
    for filepair in filepairs:
        sample_dir = os.path.commonprefix([os.path.basename(filepair[0]),os.path.basename(filepair[1])])[:-2]
        os.makedirs(os.path.join(inpath,sample_dir), exist_ok=True)
        for file in filepair:
            file = os.path.join(os.getcwd(),file)
            subprocess.check_output(f'ln -s -f {file} {os.path.join(inpath,sample_dir,os.path.basename(file))}', shell=True)
    if os.path.exists(outpath):
        os.rmdir(outpath)
    subprocess.check_output(f'HiC-Pro -i {inpath} -o {outpath} -c {config}', shell=True)


### run macs2
def run_macs2(rep, outpath, adapter, minlen=20, utrim=0, dtrim=0, threads=6, genome_index='/home/ljw/hg19_with_bowtie2_index/hg19', qvalue=0.001):
    if dtrim>0:
        raise Exception("dtrim must be negative (trim reversely from 3' end)")
    for file in rep:
        subprocess.check_output(f'cutadapt -u {utrim} -u {dtrim} -a {adapter} -m {minlen} -j {threads} {file} -o {f"{file}.adp"}', shell=True)
        subprocess.check_output(f'bowtie2 -p {threads} -x {genome_index} -U {f"{file}.adp"} | samtools view -@ {threads} -b | samtools sort -@ {threads} > {f"{file}.adp.bam"}', shell=True)
        pysam.index("-b", "-@", f"{threads}", f"{file}.adp.bam")
        deeptools.bamCoverage.main(["-b", f"{file}.adp.bam", "-o", f"{file}.adp.bam.bw"])
    out_name = os.path.join(os.path.dirname(rep[0]), os.path.commonprefix([os.path.basename(rep[0]),os.path.basename(rep[1])])[:-4])
    pysam.merge("-@", f"{threads}", "-f", "-o", f"{out_name}.adp.bam", f"{rep[0]}.adp.bam", f"{rep[1]}.adp.bam")
    pysam.sort("-@", f"{threads}", "-o", f"{out_name}.adp.sort.bam", f"{out_name}.adp.bam")
    subprocess.check_output(f"mv {out_name}.adp.sort.bam {out_name}.adp.bam", shell=True)
    pysam.index("-b", "-@", f"{threads}", f"{out_name}.adp.bam")
    deeptools.bamCoverage.main(["--numberOfProcessors", f"{threads}", "-b", f"{out_name}.adp.bam", "-o", f"{out_name}.adp.bam.bw"])
    for inputs, on in zip([f"{rep[0]}.adp.bam {rep[1]}.adp.bam", f"{rep[0]}.adp.bam", f"{rep[1]}.adp.bam"],[out_name,f"{out_name}_rep1",f"{out_name}_rep2"]):
        subprocess.check_output(f'macs2 callpeak -t {inputs} -g hs -n {on} -B -q {qvalue} -f BAM --out {outpath} -m 1 50', shell=True)

### run hichipper
def run_hichipper(hichipper_outpath, hicpro_outpath, macs2_outpath, resfrags_file):
    shutil.rmtree(os.path.join(hichipper_outpath,"hicpro_div"))
    commons = os.listdir(os.path.join(hicpro_outpath,"bowtie_results","bwt2"))
    for common in commons:
        os.makedirs(os.path.join(hichipper_outpath,"hicpro_div",common,"bowtie_results","bwt2"), exist_ok=True)
        os.makedirs(os.path.join(hichipper_outpath,"hicpro_div",common,"hic_results","data"), exist_ok=True)
        os.makedirs(os.path.join(hichipper_outpath,"hicpro_div","yamls"), exist_ok=True)
        with open(os.path.join(hichipper_outpath,"hicpro_div","yamls",f"{common}.yaml"), "w") as f:
            f.write(f'peaks:\n - {os.path.join(os.getcwd(),macs2_outpath,common.replace("HiChIP", "ChIP-nexus"))}_peaks.narrowPeak\nresfrags:\n - {os.path.join(os.getcwd(),resfrags_file)}\nhicpro_output:\n - {os.path.join(os.getcwd(),hichipper_outpath,"hicpro_div",common)}\n')
        subprocess.check_output(f'ln -f -s {os.path.join(os.getcwd(),hicpro_outpath,"bowtie_results","bwt2",common)} {os.path.join(hichipper_outpath,"hicpro_div",common,"bowtie_results","bwt2",common)}', shell=True)
        subprocess.check_output(f'ln -f -s {os.path.join(os.getcwd(),hicpro_outpath,"hic_results","data",common)} {os.path.join(hichipper_outpath,"hicpro_div",common,"hic_results","data",common)}', shell=True)

    for common in commons:
        if os.path.exists(os.path.join(hichipper_outpath, common)):
            os.rmdir(os.path.join(hichipper_outpath, common))

    current_path = os.getcwd()
    os.chdir(hichipper_outpath)
    for common in commons:
        subprocess.check_output(f'hichipper --out {common} {os.path.join("hicpro_div","yamls",f"{common}.yaml")}', shell=True)
    os.chdir(current_path)