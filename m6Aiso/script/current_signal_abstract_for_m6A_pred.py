import os
import re
import sys
import gzip
import argparse
import pandas as pd
from collections import defaultdict as dd 

# python NanopolishSplit_Normalized.py -i nanopolish.tsv.gz -n 2 -o ./outdir/
def splitnanopolish_candidatecurrent(nanopolishfile, number, outpath):
    #path = os.getcwd()

    ##### split the nanopolish file
    with gzip.open(nanopolishfile, 'r') as inputfile:
        i = 0
        contiglist = {}
        contigread = []

        for line in inputfile:
            line = line.decode()
            if line.startswith('contig'):
                header = line
                continue
            line = line.rstrip()
            line_x = line.split('\t')
            if len(line_x) != 15:  #in order to delete some error line
                continue
            contigread.append(line)   
            if len(contiglist.keys()) < int(number):               
                #print(len(contiglist.keys()))
                if line_x[0] not in contiglist:
                    contiglist[line_x[0]] = []
                else:
                    continue
            else:
                tmpline1 = contigread[-2]
                tmpline2 = contigread[-1]
                with open(outpath + '/temp_' + str(i) + '.tmp', 'w') as f:
                    f.write(header)
                    for row in contigread[:-2]:
                        f.write(row + '\n')

                i = i + 1
                contiglist = {}
                contigread = []
                contigread.append(tmpline1)
                contigread.append(tmpline2)

        with open(outpath + '/temp_' + str(i) + '.tmp', 'w') as f:
            f.write(header)
            for row in contigread:
                f.write(row + '\n')               
        contiglist = {}
        contigread = []
    
    ## normalized the nanopolish current intensity
    for num in range(i+1):
        batchfile = outpath + '/temp_' + str(num) + '.tmp'
        print(batchfile)
        normlize_nanopolish(batchfile)

        batch = batchfile + '_' + 'tmp.normalized.tsv'
        candidatepos, contigline = get_center_position(batch)
        centerflank = candidate_position_flanking(candidatepos)
        centerregion = get_candidate_region(centerflank, contigline)
        candidate_position_filtering(centerregion, outpath)



def normlize_nanopolish(batchfile):   
    #(filepath, tempfilename) = os.path.split(batchfile)
    #(filename, _) = os.path.splitext(tempfilename)

    inputfile = open(batchfile, 'rt')
    eventalign_result = pd.read_csv(inputfile, delimiter='\t',names=['contig','position','reference_kmer','read_index','strand','event_index','event_level_mean','event_stdv','event_length','model_kmer','model_mean','model_stdv','standardized_level','start_idx','end_idx'],low_memory=False)
    cond_successfully_eventaligned = eventalign_result['reference_kmer'] == eventalign_result['model_kmer']
    if cond_successfully_eventaligned.sum() != 0:

        eventalign_result = eventalign_result[cond_successfully_eventaligned]

        keys = ['read_index','contig','position','reference_kmer'] # for groupby
        eventalign_result.loc[:, 'length'] = pd.to_numeric(eventalign_result['end_idx'])-pd.to_numeric(eventalign_result['start_idx'])
        eventalign_result.loc[:, 'sum_norm_mean'] = pd.to_numeric(eventalign_result['event_level_mean']) * eventalign_result['length']
        eventalign_result.loc[:, 'sum_norm_std'] = pd.to_numeric(eventalign_result['event_stdv']) * eventalign_result['length']
        eventalign_result.loc[:, 'sum_dwell_time'] = pd.to_numeric(eventalign_result['event_length']) * eventalign_result['length']
            
        eventalign_result = eventalign_result.groupby(keys)  
        sum_norm_mean = eventalign_result['sum_norm_mean'].sum() 
        sum_norm_std = eventalign_result["sum_norm_std"].sum()
        sum_dwell_time = eventalign_result["sum_dwell_time"].sum()

        start_idx = eventalign_result['start_idx'].min()
        end_idx = eventalign_result['end_idx'].max()
        total_length = eventalign_result['length'].sum()

        eventalign_result = pd.concat([start_idx,end_idx],axis=1)
        eventalign_result['norm_mean'] = (sum_norm_mean/total_length).round(1)
        eventalign_result["norm_std"] = sum_norm_std / total_length
        eventalign_result["dwell_time"] = sum_dwell_time / total_length
        eventalign_result.reset_index(inplace=True)

        eventalign_result['transcript_id'] = [contig.split('.')[0] for contig in eventalign_result['contig']]    #### CHANGE MADE ####

        eventalign_result['transcriptomic_position'] = pd.to_numeric(eventalign_result['position']) + 2 # the middle position of 5-mers.
        features = ['transcript_id','read_index','transcriptomic_position','reference_kmer','norm_mean','norm_std','dwell_time']
        df_events = eventalign_result[features]
        
        
        outname = open(batchfile + '_' + 'tmp.normalized.tsv', 'w')
        
        #if not os.path.getsize('CurrentIntensity_normalized.tsv'):
        df_events.to_csv(outname, sep='\t', index=False, header=True)
        # header: 'transcript_id','read_index','transcriptomic_position','reference_kmer','norm_mean','norm_std','dwell_time'
        #else:
        # df_events.to_csv(outname, sep='\t', index=False, header=False)
        os.remove(batchfile)

################################################################################################################################################################
################################################################################################################################################################

# get the DRACH motif position
def get_center_position(batchfile):
    candidatepos = {}
    contigline = {}
    
    with open(batchfile, 'r') as inputfile:
        for line in inputfile:
            line = line.strip('\n')
            line_x = line.split('\t')
            if line.startswith('transcript_id'):
                idx = {i:j for j,i in enumerate(line_x)}
                #header = line
                continue
            contig = line_x[idx['transcript_id']]
            read = line_x[idx['read_index']]
            position = line_x[idx['transcriptomic_position']].split('.')[0]
            kmer = line_x[idx['reference_kmer']]
            mean = line_x[idx['norm_mean']]
            std = line_x[idx['norm_std']]
            dwell = line_x[idx['dwell_time']]
            
            contigline[(contig, read, position)] = [mean, std, dwell, kmer]
            #key = '{}_{}_{}'.format(contig, read, position)
            if len(re.findall('[AGT][AG]AC[ACT]', kmer)) != 0:
                candidatepos[(contig, read, position)] = kmer
            else:
                continue
    os.remove(batchfile)
    return candidatepos, contigline


# flank the DRACH motifs
def candidate_position_flanking(candidatepos):
    centerflank = {}
    for key in candidatepos:
        contig, read, position = key[:]
        kmer = candidatepos[key]
        for i in range(int(position)-1, int(position)+2, 1):
            distance = int(position) - i
            if (contig, read, str(i)) not in centerflank:
                centerflank[(contig, read, str(i))] = [kmer, str(distance)]
            else:
                continue
    return centerflank

# get the flanking regions
def get_candidate_region(centerflank, contigline):
    centerregion = dd(lambda: dd(lambda: ''))
    for key in centerflank:
        contig, read, position = key[:]
        kmer, distance = centerflank[key][:]
        centerpos = str(int(position) + int(distance))
        if key in contigline:
            centerregion[(contig, read, centerpos)][distance] = contigline[key]
            centerregion[(contig, read, centerpos)]['kmer'] = kmer
        else:
            continue
    return centerregion

# process the candidate position
def candidate_position_filtering(centerregion, outpath):
    header = ['kmer', 'kmer_contig_readindex_tranpos', 'P1_mean', 'P1_std', 'P1_length', 'P0_mean', 
              'P0_std', 'P0_length', 'N1_mean', 'N1_std', 'N1_length', 'baseflank']

    candidatefilter = []

    relative = [-1, 0, 1]
    relative = [str(i) for i in relative]
    for key in centerregion:
        distance = []
        for kk in centerregion[key]:
            distance.append(kk)
        #print(distance)

        #idname = '{}_{}_{}'.format(key[0], key[1], key[2])
        if set(relative)<=set(distance):
            kmer = centerregion[key]['kmer']
            idname = '{}_{}_{}_{}'.format(kmer, key[0], key[1], key[2])

            normvalues =  centerregion[key]['-1'][:3] + centerregion[key]['0'][:3] + centerregion[key]['1'][:3]
            flankkmer = centerregion[key]['1'][-1][0] + centerregion[key]['0'][-1] + centerregion[key]['-1'][-1][-1]

            candidatefilter = '\t'.join(['{}'.format(i) for i in ([kmer, idname] + normvalues)])
            

            with open(outpath + '/Candidatecurrent.tsv', 'a') as f:
            #with open('Candidatecurrent.tsv', 'a') as f:
                if not os.path.getsize(outpath + '/Candidatecurrent.tsv'):
                #if not os.path.getsize('candidatecurrent.tsv'):
                    f.write('\t'.join(['{}'.format(i) for i in header]) + '\n')
                
                f.write(candidatefilter + '\t' + flankkmer + '\n')

################################################################################################################################################################
################################################################################################################################################################

def args_make(parser):
    parser.add_argument("--nanopolish_result",
                        help="nanopolish eventalign generated current result.",
                        metavar='\b',
                        required=True)
    
    parser.add_argument('--number', 
                        default=2,
                        type=int,
                        help='The number of each batchs (transcript)',
                        metavar='\b',
                        required=True)
    
    parser.add_argument('--out_dir',
                        help='The output directory',
                        metavar='\b',
                        required=True)
    return parser

def main(args):
    nanopolishfile = args.input_nanopolish
    number = args.number
    outpath = args.out_dir
    splitnanopolish_candidatecurrent(nanopolishfile, number, outpath)

