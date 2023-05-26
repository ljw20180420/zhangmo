import ljwlib.hic_module, pandas, requests, pathlib

def trackData2DataFrame(trackData):
    print(type(trackData))
    if isinstance(trackData[trackData['track']], dict):
        print(trackData.keys())
    if trackData[trackData['track']]==[]:
        return pandas.DataFrame({'chrom' : [], 'start' : [], 'end' : [], 'name' : [], 'score' : [], 'strand' : []})
    if isinstance(trackData[trackData['track']], dict):
        fragments = trackData[trackData['track']][trackData['chrom']]
    else:
        fragments =  trackData[trackData['track']]
    for schema in ['chrom', 'chromStart', 'chromEnd']:
        if schema not in fragments[0].keys():
            return None
    chroms, starts, ends, names, scores, strands = [], [], [], [], [], []
    for fragment in fragments:
        chroms.append(fragment.get('chrom'))
        starts.append(fragment.get('chromStart'))
        ends.append(fragment.get('chromEnd'))
        names.append('.' if fragment.get('name') is None else  fragment.get('name'))
        scores.append('.' if fragment.get('score') is None else fragment.get('score'))
        strands.append('.' if fragment.get('strand') is None else fragment.get('strand'))
    return pandas.DataFrame({'chrom' : chroms, 'start' : starts, 'end' : ends, 'name' : names, 'score' : scores, 'strand' : strands})

def iter_trackData(in_dict, func, loghandle, errhandle, donekeys, failkeys):
    params = {'genome': 'hg19', 'track': None, 'maxItemsOutput': None, 'chrom': None, 'chromStart': None, 'chromEnd': None, 'hubUrl': None}
    records = ['shortLabel', 'type', 'longLabel', 'parent', 'parentParent', 'origAssembly', 'subGroups']
    subGroups = [f'subGroup{i}' for i in range(10)]
    skipkeys = ['compositeContainer', 'sortOrder', 'controlledVocabulary', 'group', 'compositeTrack', 'dimensions', 'dragAndDrop', 'fileSortOrder', 'filterBy', 'itemRgb', 'priority', 'superTrack', 'noScoreFilter', 'visibility', 'pennantIcon', 'color', 'html'] + subGroups + records
    for key in in_dict.keys():
        if key in skipkeys or key in donekeys:
            continue
        if key not in failkeys:
            params['track'] = key
            failkey = False
            df_segs = []
            tlen = 0
            for chrom in ljwlib.hic_module.chromsizes.keys():
                params['chrom'] = chrom
                try:
                    trackData = requests.get('https://genome-asia.ucsc.edu/cgi-bin/hubApi/getData/track', params).json()
                except Exception as err:
                    errhandle.write(f'{key}\t{repr(err)}\n')
                    failkey = True
                    break
                if 'error' in trackData.keys():
                    errhandle.write(f'{key}\t{trackData["error"]}\n')
                    failkey = True
                    break
                df_seg = trackData2DataFrame(trackData)
                if not isinstance(df_seg, pandas.DataFrame):
                    errhandle.write(f'{key}\t{"no_chrom_start_end"}\n')
                    failkey = True
                    break
                df_segs.append(df_seg)
                tlen += len(df_seg)
                if tlen>1e6:
                    errhandle.write(f'{key}\t{"data_too_large"}\n')
                    failkey = True
                    break
        if key in failkeys or failkey:
            if isinstance(in_dict[key], dict):
                iter_trackData(in_dict[key], func, loghandle, errhandle, donekeys, failkeys)
        else:
            loghandle.write(f'{key}\t')
            for record in records:
                loghandle.write(f'{in_dict[key].get(record)}\t')
            df = pandas.concat(df_segs).reset_index(drop=True)
            df[['start','end']] = df[['start','end']].astype("int64")
            loghandle.write(f'{func(df)}\n')

def scrab_ucsc(logfile, errfile, func):
    ### create logfile and errfile if not exist
    pathlib.Path(logfile).touch(exist_ok=True)
    pathlib.Path(errfile).touch(exist_ok=True)
    ### clean the tail bad data
    donekeys = []
    failkeys = []
    for file in [logfile, errfile]:
        with open(file, 'r') as handle:
            lines = handle.readlines()
        if not lines:
            continue
        if len(lines[0].split('\t'))>len(lines[-1].split('\t')) or not lines[-1].split('\t')[-1]:
            lines.pop(-1)
            with open(file, 'w') as handle:
                handle.writelines(lines)
        if file==logfile:
            donekeys = set([line.split('\t')[0] for line in lines])
        else:
            failkeys = set([line.split('\t')[0] for line in lines])
    ### iterate over tracks
    track_dict = requests.get('https://genome-asia.ucsc.edu/cgi-bin/hubApi/list/tracks', {'genome': 'hg19'}).json()['hg19']
    with open(logfile, 'a', buffering=1) as loghandle, open(errfile, 'a', buffering=1) as errhandle:
        iter_trackData(track_dict, func, loghandle, errhandle, donekeys, failkeys)



### ucsc functional
def colocalize_count(bed1, bed2):
    bed1_bed2 = ljwlib.hic_module.standard_assign(bed1, bed2, '', 'assign', select_col='distance', includes=['name'], select_mode='min')
    bed2_bed1 = ljwlib.hic_module.standard_assign(bed2, bed1, '', 'assign', select_col='distance', includes=['name'], select_mode='min')
    # colocals = ['bed1_only', 'bed1_over', 'bed2_only', 'bed2_over']
    bed1_only_count = sum(bed1_bed2['assign_name'].isna())
    bed2_only_count = sum(bed2_bed1['assign_name'].isna())
    counts = [bed1_only_count, len(bed1)-bed1_only_count, bed2_only_count, len(bed2)-bed2_only_count]
    # return pandas.DataFrame({'colocal' : colocals, 'count' : counts})
    return counts