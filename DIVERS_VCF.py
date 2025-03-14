# python3.8
__author__ =    "Peng Zhang"
__copyright__ = "Copyright 2023, " \
                "Laboratory of Human Genetics of Infectious Diseases, " \
                "The Rockefeller University"
__license__ =   "CC BY-NC-ND 4.0"
__version__ =   "ver-1, 03-10-2025"

import os
import re
import time
import argparse

print('******************************************')
print(' ####   #####  #   #  #####  ####   ##### ')
print(' #   #    #    #   #  #      #   #  ##    ')
print(' #   #    #    #   #  ###    ####     #   ')
print(' #   #    #     # #   #      #  #      ## ')
print(' ####   #####    #    #####  #   #  ##### ')
print('          Deep Intronic Variants          ')
print('     with Effect on Recursive Splicing    ')
print('*******************************************\n')

###
# input parameters
###
timestamp = str(time.strftime('%Y%m%d-%H%M%S', time.localtime()))

parser = argparse.ArgumentParser(description="DIVERS")
parser.add_argument("-i", "--input", help="input variants, in VCF format")

args = parser.parse_args()
filename_var = args.input
filename_bed = filename_var+'.bed'
filename_map = filename_var+'.mapping'
filename_out = filename_var+'.DIVERS_'+timestamp+'.csv'

base_pairing = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}
def score_5ss(motif):
    if (motif[0] == 'G') and (motif[1] == 'G') and (motif[2] == 'T') and \
       (motif[3] in ['A','G']) and (motif[4] in ['A','G']) and (motif[5] == 'G'):
        return 3
    elif (motif[0] == 'G') and (motif[1] == 'G') and (motif[2] == 'T') and \
         (motif[3] in ['A','G']) and (motif[4] in ['A','G']) and (motif[5] in ['A','G']):
        return 2
    elif (motif[1] == 'G') and (motif[2] == 'T') and \
         (motif[3] in ['A','G']) and (motif[4] in ['A','G']) and (motif[5] in ['A','G']):
        return 1
    else:
        return 0

var_input_set = set()
var_output_set = set()
var_output_annot = 0
var_count_AGGT = 0
var_count_BP = 0
var_count_BP2 = 0
var_count_AGAIN = 0
var_count_DW5SS = 0
var_count_CRYPRS = 0

try:
    file_var = open(filename_var, 'r')
    file_bed = open(filename_bed, 'w')
    var_header = file_var.readline().strip()
    var_header_out = ','.join(var_header.split('\t')[5:])

    for eachline in file_var:
        if not eachline.startswith('#'):
            item = eachline.strip().split('\t')
            chrom = item[0]
            pos = item[1]
            var_id = item[2]
            ref = item[3]
            alt = item[4]
            if 'chr' not in chrom:
                chrom = 'chr' + chrom
            var_name = chrom+'*'+pos+'*'+var_id+'*'+ref+'*'+alt
            var_info = ','.join(item[5:])
            if not var_info:
                var_info = '.'
            var_start = var_end = '.'
            if (ref != '.') and (alt != '.') and (len(ref) == 1) and (len(alt) == 1) and (var_name not in var_input_set):
                var_start = str(int(pos)-1)
                var_end = pos
                var_input_set.add(var_name)
                file_bed.write(chrom+'\t'+var_start+'\t'+var_end+'\t'+var_name+'\t.\t+\t'+var_info+'\n')
                file_bed.write(chrom+'\t'+var_start+'\t'+var_end+'\t'+var_name+'\t.\t-\t'+var_info+'\n')
    file_var.close()
    file_bed.close()

    # DIVERS mapping, detection and annotation
    os.system("bedtools intersect -wo -s -a " + filename_bed + " -b DIVERS_detection.bed > " + filename_map)
    file_map = open(filename_map, 'r')
    file_out = open(filename_out, 'w')
    file_out.write('CHR,POS,ID,REF,ALT,STRAND,GENE,TRANSCRIPT,IVS#,IVS_LENGTH,'
                   'RS_TOTAL,RS#,BP_POS,PPT_Y,RS_POS,CLIP,RS_CONSEQ,'+var_header_out+'\n')

    for eachline in file_map:
        item = eachline.strip().split('\t')
        chrom = item[0]
        var_start = int(item[1])
        var_end = int(item[2])
        var_name = item[3]
        strand = item[5]
        var_info = item[6]
        RS_element_start = int(item[8])
        RS_element_end = int(item[9])
        gene = item[13]
        transcript = item[14]
        ivs = item[15]
        ivs_length = item[16]
        RS_total = item[17]
        RS_rank = item[18]
        BP_pos = item[19]
        PPT = item[20]
        RS_pos = item[21]
        clip = item[22]
        RS_conseq = item[23]
        seq_wt = item[24]
        seq_mt = '.'

        var_name_list = var_name.split('*')
        pos = var_name_list[1]
        var_id = var_name_list[2]
        ref = var_name_list[3]
        alt = var_name_list[4]
        if var_info == '.':
            var_info = ''

        # seq analysis for AGAIN and DOWN only
        if RS_conseq in ['RS_AGAIN','RS_DW5SS']:
            nucl_list = list(seq_wt)
            if strand == '+':
                var_index = int(pos) - RS_element_start - 1
                nucl_list[var_index] = alt
            else:
                var_index = RS_element_end - int(pos)
                nucl_list[var_index] = base_pairing[alt]
            seq_mt = ''.join(nucl_list)

        DIVERS_flag = 0

        # RS_AGGT
        if RS_conseq == 'RS_AGGT':
            var_count_AGGT += 1
            DIVERS_flag = 1
            var_output_set.add(var_name)

        # RS_BP & RS_BP2
        if RS_conseq == 'RS_BP':
            var_count_BP += 1
            DIVERS_flag = 1
            var_output_set.add(var_name)

        if RS_conseq == 'RS_BP2':
            if (strand == '+' and alt in ['A','G']) or (strand == '-' and alt in ['C','T']):
                var_count_BP2 += 1
                DIVERS_flag = 1
                var_output_set.add(var_name)

        # RS_AGAIN
        if RS_conseq == 'RS_AGAIN':
            AG_hit_wt = set(hit.end() for hit in re.finditer('AG', seq_wt))
            AG_hit_mt = set(hit.end() for hit in re.finditer('AG', seq_mt))
            AG_hit_set = AG_hit_mt - AG_hit_wt
            if AG_hit_set:
                AG_hit_list = list(AG_hit_set)
                AG_hit_list.sort()
                AG_hit = AG_hit_list[-1]
                if seq_mt[AG_hit-3] in ['C','T']:
                    var_count_AGAIN += 1
                    DIVERS_flag = 1
                    var_output_set.add(var_name)

        # RS_DW5SS
        if RS_conseq == 'RS_DW5SS':
            for i in range(0, len(seq_wt)-5):
                motif_wt = seq_wt[i:i+6]
                motif_mt = seq_mt[i:i+6]
                if 'GT' in motif_mt:
                    score_wt = score_5ss(motif_wt)
                    score_mt = score_5ss(motif_mt)
                    if score_mt > score_wt:
                        RS_conseq += '_'+str(i+3)+'nt'
                        var_count_DW5SS += 1
                        DIVERS_flag = 1
                        var_output_set.add(var_name)
                        break

        # CRYP_RS
        if ('CRYPRS' in RS_conseq) and (var_name not in var_output_set):
            CRYPRS_type = seq_wt
            CRYPRS_flag = 0
            if CRYPRS_type == 'NGGT':
                if (strand == '+' and alt == 'A') or (strand == '-' and alt == 'T'):
                    CRYPRS_flag = 1
            if CRYPRS_type == 'ANGT':
                if (strand == '+' and alt == 'G') or (strand == '-' and alt == 'C'):
                    CRYPRS_flag = 1
            if CRYPRS_type == 'AGNT':
                if (strand == '+' and alt == 'G') or (strand == '-' and alt == 'C'):
                    CRYPRS_flag = 1
            if CRYPRS_type == 'AGGN':
                if (strand == '+' and alt == 'T') or (strand == '-' and alt == 'A'):
                    CRYPRS_flag = 1
            if CRYPRS_flag:
                var_count_CRYPRS += 1
                DIVERS_flag = 1
                var_output_set.add(var_name)

        # DIVERS Output
        if DIVERS_flag:
            var_output_annot += 1
            file_out.write(chrom+','+pos+','+var_id+','+ref+','+alt+','+strand+','+
                           gene+','+transcript+','+ivs+','+ivs_length+','+RS_total+','+RS_rank+','+
                           BP_pos+','+PPT+','+RS_pos+','+clip+','+RS_conseq+','+var_info+'\n')

    file_out.close()
    file_map.close()
    os.remove(filename_bed)
    os.remove(filename_map)

    print('Input file:', filename_var)
    print('Output file:', filename_out, '\n')
    print('# Input variants:', str(len(var_input_set)))
    print('# Output variants:', str(len(var_output_set)))
    print('# Output annotations:', str(var_output_annot), '\n')
    print('# RS_AGGT:\t', str(var_count_AGGT))
    print('# RS_BP:\t', str(var_count_BP))
    print('# RS_BP-2:\t', str(var_count_BP2))
    print('# RS_AGAIN:\t', str(var_count_AGAIN))
    print('# RS_DW5SS:\t', str(var_count_DW5SS))
    print('# CRYP_RS:\t', str(var_count_CRYPRS), '\n')

except:
    print('Error occured. Please check your input file, or contact pzhang@rockefeller.edu.\n')
