# python3.8
__author__ =    "Peng Zhang"
__copyright__ = "Copyright 2025, " \
                "Laboratory of Human Genetics of Infectious Diseases, The Rockefeller University"
__license__ =   "CC BY-NC-ND 4.0"
__version__ =   "ver-1, 05-01-2025"

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
def rev_seq(fwd_seq):
    out_seq = ''
    for fwd_pos in range(0, len(fwd_seq)):
        out_seq += base_pairing[fwd_seq[len(fwd_seq) - fwd_pos - 1]]
    return out_seq

var_input_set = set()
var_output_set = set()
var_output_annot = 0
var_count_RS_AGGT = 0
var_count_RS_BP = 0
var_count_RS_BP2 = 0
var_count_RS_AGAIN = 0
var_count_RS_DW5SS = 0
var_count_CRYPRS = 0

###
# DIVERS running
###
try:
    file_var = open(filename_var, 'r')
    file_bed = open(filename_bed, 'w')
    var_info_header = ','.join(file_var.readline().strip().split('\t')[5:])

    for eachline in file_var:
        if not eachline.startswith('#'):
            eachline = eachline.replace(',', ';')
            eachline = eachline.replace('*', ';')
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

            var_start = var_end = var_type = '.'
            if ref.isalpha() and alt.isalpha():
                if len(ref) == 1 and len(alt) == 1:
                    var_type = 'snv'
                    var_start = str(int(pos)-1)
                    var_end = pos
                elif 10 >= len(ref) > 1 == len(alt):
                    var_type = 'del-'+str(len(ref)-1)+'nt'
                    var_start = pos
                    var_end = str(int(pos)+len(ref)-1)
                elif 10 >= len(alt) > 1 == len(ref):
                    var_type = 'ins-'+str(len(alt)-1)+'nt'
                    var_start = pos
                    var_end = str(int(pos)+len(alt)-1)

            if (var_type != '.') and (var_name not in var_input_set):
                var_input_set.add(var_name)
                file_bed.write(chrom+'\t'+var_start+'\t'+var_end+'\t'+var_name+'\t.\t+\t'+var_type+'\t'+var_info+'\n')
                file_bed.write(chrom+'\t'+var_start+'\t'+var_end+'\t'+var_name+'\t.\t-\t'+var_type+'\t'+var_info+'\n')

    file_var.close()
    file_bed.close()

    # DIVERS mapping, detection and annotation
    os.system("bedtools intersect -wo -s -a " + filename_bed + " -b DIVERS_Detection.bed > " + filename_map)
    file_map = open(filename_map, 'r')
    file_out = open(filename_out, 'w')
    file_out.write('CHR,POS,ID,REF,ALT,STRAND,VAR_TYPE,GENE,TRANSCRIPT,IVS#,IVS_SIZE,'
                   'RS#,RS_CONSEQ,RS_SCORE,RS_POS,BP_POS,PPT_Y,CLIP,RARE,CONSERV,RNALM,'+var_info_header+'\n')

    for eachline in file_map:
        item = eachline.strip().split('\t')
        chrom = item[0]
        var_start = int(item[1])
        var_end = int(item[2])
        var_name = item[3]
        strand = item[5]
        var_type = item[6]
        var_info = item[7]
        RS_region_start = int(item[9])
        RS_region_end = int(item[10])
        gene = item[14]
        transcript = item[15]
        ivs = item[16]
        ivs_size = item[17]
        RS = item[18]
        BP_pos = item[19]
        PPT = item[20]
        RS_pos = item[21]
        clip = item[22]
        rare = item[23]
        conserv = item[24]
        rnalm = item[25]
        RS_score = item[26]
        RS_conseq = item[27]
        seq_wt = item[28]
        seq_mt = '.'

        var_name_item = var_name.split('*')
        var_pos = var_name_item[1]
        var_id = var_name_item[2]
        var_ref = var_name_item[3]
        var_alt = var_name_item[4]
        if var_info == '.':
            var_info = ''

        DIVERS_flag = 0

        # RS-AGGT
        if RS_conseq == 'RS-AGGT':
            DIVERS_flag = 1
            var_count_RS_AGGT += 1
            var_output_set.add(var_name)

        # RS-BP/BP2
        if RS_conseq == 'RS-BP':
            DIVERS_flag = 1
            var_count_RS_BP += 1
            var_output_set.add(var_name)
        if RS_conseq == 'RS-BP2':
            if var_type == 'snv':
                if ((strand == '+' and var_alt in ['A','G']) or
                    (strand == '-' and var_alt in ['C','T'])):
                    DIVERS_flag = 1
                    var_count_RS_BP2 += 1
                    var_output_set.add(var_name)
            else:
                DIVERS_flag = 1
                var_count_RS_BP2 += 1
                var_output_set.add(var_name)

        # RS-AGAIN
        if RS_conseq == 'RS-AGAIN':
            nucl_list = list(seq_wt)
            if var_type == 'snv':
                if strand == '+':
                    var_index = int(var_pos)-RS_region_start-1
                    nucl_list[var_index] = var_alt
                else:
                    var_index = RS_region_end-int(var_pos)
                    nucl_list[var_index] = base_pairing[var_alt]
            elif 'ins' in var_type:
                if strand == '+':
                    var_index = int(var_pos)-RS_region_start-1
                    nucl_list[var_index] = var_alt
                else:
                    var_index = RS_region_end-int(var_pos)
                    nucl_list[var_index] = rev_seq(var_alt)
            elif ('del' in var_type) and (RS_region_end-len(var_ref) > int(var_pos) > RS_region_start):
                if strand == '+':
                    var_index = int(var_pos)-RS_region_start-1
                    for temp_index in range(var_index+1, var_index+len(var_ref)):
                        nucl_list[temp_index] = ''
                else:
                    var_index = RS_region_end-int(var_pos)
                    for temp_index in range(var_index-len(var_ref)+1, var_index):
                        nucl_list[temp_index] = ''

            seq_mt = ''.join(nucl_list)
            AG_hit_wt = set(hit.end() for hit in re.finditer('AG', seq_wt))
            AG_hit_mt = set(hit.end() for hit in re.finditer('AG', seq_mt))
            AG_hit_set = AG_hit_mt - AG_hit_wt

            if 'del' in var_type:
                AG_hit_set_temp = set()
                for each_pos in AG_hit_set:
                    if (each_pos+len(var_ref)-1) not in AG_hit_wt:
                        AG_hit_set_temp.add(each_pos)
                AG_hit_set = AG_hit_set_temp
            if 'ins' in var_type:
                AG_hit_set_temp = set()
                for each_pos in AG_hit_set:
                    if (each_pos-len(var_alt)+1) not in AG_hit_wt:
                        AG_hit_set_temp.add(each_pos)
                AG_hit_set = AG_hit_set_temp
            if AG_hit_set:
                AG_hit_list = list(AG_hit_set)
                AG_hit_list.sort()
                AG_hit = AG_hit_list[-1]
                if seq_mt[AG_hit-3] in ['C','T']:
                    DIVERS_flag = 1
                    var_count_RS_AGAIN += 1
                    var_output_set.add(var_name)

        # RS-DW5SS
        if RS_conseq.startswith('RS-DW5SS'):
            seq_wt_item = seq_wt.split('|')
            check_ref = seq_wt_item[0]
            check_alt = seq_wt_item[1]
            if var_ref == check_ref and var_alt == check_alt:
                DIVERS_flag = 1
                var_count_RS_DW5SS += 1
                var_output_set.add(var_name)

        # CRYPRS
        if ('CRYPRS' in RS_conseq) and (var_name not in var_output_set):
            CRYPRS_type = seq_wt
            CRYPRS_flag = 0
            if CRYPRS_type == 'NGGT':
                if (strand == '+' and var_alt == 'A') or (strand == '-' and var_alt == 'T'):
                    CRYPRS_flag = 1
            if CRYPRS_type == 'ANGT':
                if (strand == '+' and var_alt == 'G') or (strand == '-' and var_alt == 'C'):
                    CRYPRS_flag = 1
            if CRYPRS_type == 'AGNT':
                if (strand == '+' and var_alt == 'G') or (strand == '-' and var_alt == 'C'):
                    CRYPRS_flag = 1
            if CRYPRS_type == 'AGGN':
                if (strand == '+' and var_alt == 'T') or (strand == '-' and var_alt == 'A'):
                    CRYPRS_flag = 1
            if CRYPRS_flag:
                DIVERS_flag = 1
                var_count_CRYPRS += 1
                var_output_set.add(var_name)

        # DIVERS Output
        if DIVERS_flag:
            var_output_annot += 1
            file_out.write(chrom+','+var_pos+','+var_id+','+var_ref+','+var_alt+','+strand+','+var_type+','+
                           gene+','+transcript+','+ivs+','+ivs_size+','+RS+','+RS_conseq+','+RS_score+','+
                           RS_pos+','+BP_pos+','+PPT+','+clip+','+rare+','+conserv+','+rnalm+','+var_info+'\n')

    file_out.close()
    file_map.close()
    os.remove(filename_bed)
    os.remove(filename_map)

    print('Input file:', filename_var)
    print('Output file:', filename_out, '\n')
    print('# Input variants:', str(len(var_input_set)))
    print('# Output variants:', str(len(var_output_set)))
    print('# Output annotations:', str(var_output_annot), '\n')
    print('# RS-AGGT:\t', str(var_count_RS_AGGT))
    print('# RS-BP:\t', str(var_count_RS_BP))
    print('# RS-BP2:\t', str(var_count_RS_BP2))
    print('# RS-AGAIN:\t', str(var_count_RS_AGAIN))
    print('# RS-DW5SS:\t', str(var_count_RS_DW5SS))
    print('# CRYPRS:\t', str(var_count_CRYPRS), '\n')

except Exception as error_message:
    print(f"An error occurred: {error_message}")
    print("Please check your input file, or contact pzhang@rockefeller.edu.\n")
