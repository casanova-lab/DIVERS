# python3.8
__author__ =    "Peng Zhang"
__copyright__ = "Copyright 2026, " \
                "Laboratory of Human Genetics of Infectious Diseases"
__license__ =   "CC BY-NC-ND 4.0"
__version__ =   "02-28-2026"

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
timestamp = str(time.strftime('%Y%m%d_%H%M%S', time.localtime()))

parser = argparse.ArgumentParser(description="DIVERS - batch input")
parser.add_argument("-d", "--dir", help="directory of VCF files")
parser.add_argument("-s", "--sample", help="sample list")
parser.add_argument("-o", "--output", help="output filename, in CSV format")

args = parser.parse_args()
directory = args.dir
filename_sample = args.sample
filename_out = args.output

if filename_out.endswith('.csv'):
    file_out = open(filename_out, 'w')
else:
    file_out = open(filename_out+'.csv', 'w')

base_pairing = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}
def rev_seq(fwd_seq):
    out_seq = ''
    for fwd_pos in range(0, len(fwd_seq)):
        out_seq += base_pairing[fwd_seq[len(fwd_seq) - fwd_pos - 1]]
    return out_seq

file_info = open('DIVERS_Detection_INFO.txt', 'r')
RS_info_dict = dict()
for eachinfo in file_info:
    info_item = eachinfo.strip().split('\t')
    RS_name = info_item[0]
    RS_info = ','.join(info_item[1:])
    RS_info_dict[RS_name] = RS_info

file_sample = open(filename_sample, 'r')
first_read = True
for eachsample in file_sample:
    sample = eachsample.strip()
    filename_var = sample+'.vcf'
    filename_bed = filename_var+'.bed'
    filename_map_RS = filename_var+'.map1'
    filename_map_RSCRYP = filename_var+'.map2'

    var_input_set = set()
    var_output_set = set()
    var_output_annot = 0
    var_count_RS_AGGT = 0
    var_count_RS_BP = 0
    var_count_RS_3SS = 0
    var_count_RS_5SS = 0
    var_count_RS_CRYP = 0

    ###
    # DIVERS running
    ###
    try:
        file_var = open(directory + filename_var, 'r')
        file_bed = open(filename_bed, 'w')

        if first_read:
            var_info_header = ','.join(file_var.readline().strip().split('\t')[5:]) or 'OTHERS'
            file_out.write('SAMPLE,CHR,POS,ID,REF,ALT,STRAND,VAR_TYPE,GENE,TRANSCRIPT,IVS#,RS#,RS_CONSEQ,'
                           'SCORE,IVS_LEN,BP_POS,PPT_Y,RS_START,RS_END,CLIP,RNASEQ,RNALM,PHYLOP,PHASTCONS,RARE,'+var_info_header+'\n')
            first_read = False
        else:
            file_var.readline()

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
                chrom = chrom if 'chr' in chrom else 'chr' + chrom
                var_name = chrom+'*'+pos+'*'+var_id+'*'+ref+'*'+alt
                var_info = ','.join(item[5:]) or '.'
                var_start = var_end = var_type = '.'
                if ref.isalpha() and alt.isalpha():
                    if len(ref) == 1 and len(alt) == 1:
                        var_type = 'snv'
                        var_start = str(int(pos) - 1)
                        var_end = pos
                    elif 10 >= len(ref) > 1 == len(alt):
                        var_type = 'del-' + str(len(ref) - 1) + 'nt'
                        var_start = pos
                        var_end = str(int(pos) + len(ref) - 1)
                    elif 10 >= len(alt) > 1 == len(ref):
                        var_type = 'ins-' + str(len(alt) - 1) + 'nt'
                        var_start = pos
                        var_end = str(int(pos) + len(alt) - 1)
                if (var_type != '.') and (var_name not in var_input_set):
                    var_input_set.add(var_name)
                    file_bed.write(chrom+'\t'+var_start+'\t'+var_end+'\t'+var_name+'\t.\t+\t'+var_type+'\t'+var_info+'\n')
                    file_bed.write(chrom+'\t'+var_start+'\t'+var_end+'\t'+var_name+'\t.\t-\t'+var_type+'\t'+var_info+'\n')
        file_var.close()
        file_bed.close()

        # DIVERS mapping, detection and annotation
        os.system("bedtools intersect -wo -s -a "+filename_bed+" -b DIVERS_Detection_RS.bed > "+filename_map_RS)
        os.system("bedtools intersect -wo -s -a "+filename_bed+" -b DIVERS_Detection_RSCRYP.bed > "+filename_map_RSCRYP)

        file_map_RS = open(filename_map_RS, 'r')
        for eachline in file_map_RS:
            item = eachline.strip().split('\t')
            chrom = item[0]
            var_name = item[3]
            strand = item[5]
            var_type = item[6]
            var_info = item[7]
            region_start = int(item[9])
            region_end = int(item[10])
            RS_name = item[11]
            RS_conseq = item[14]
            seq_wt = item[15]

            _, var_pos, var_id, ref, alt = var_name.split('*')
            var_info_out = var_info.replace('*', ',')
            RS_name_out = RS_name.replace('|', ',')
            RS_info_out = RS_info_dict[RS_name]
            DIVERS_flag = 0
            DIVERS_conseq = ''

            # RS_AGGT
            if RS_conseq == 'RS_AGGT':
                DIVERS_conseq = 'RS_AGGT_loss'
                DIVERS_flag = 1
                var_count_RS_AGGT += 1
                var_output_set.add(var_name)

            # RS_BP/BP2
            if RS_conseq == 'RS_BP':
                if (var_type == 'snv') or ('del' in var_type):
                    DIVERS_conseq = 'RS_BP_loss'
                    DIVERS_flag = 1
                    var_count_RS_BP += 1
                    var_output_set.add(var_name)
            if RS_conseq == 'RS_BP2':
                if var_type == 'snv':
                    if (strand == '+' and alt in 'AG') or (strand == '-' and alt in 'CT'):
                        DIVERS_conseq = 'RS_BP_loss'
                        DIVERS_flag = 1
                        var_count_RS_BP += 1
                        var_output_set.add(var_name)
                elif 'del' in var_type:
                    DIVERS_conseq = 'RS_BP_loss'
                    DIVERS_flag = 1
                    var_count_RS_BP += 1
                    var_output_set.add(var_name)

            # RS_PPT
            if RS_conseq == 'RS_PPT':
                nucl_list = list(seq_wt)
                pos = int(var_pos)
                var_index = (pos - region_start - 1) if strand == '+' else (region_end - pos)
                if var_type == 'snv':
                    nucl_list[var_index] = alt if strand == '+' else base_pairing[alt]
                elif 'ins' in var_type:
                    nucl_list[var_index] = alt if strand == '+' else rev_seq(alt)
                elif ('del' in var_type) and (region_start < pos < region_end - len(ref)):
                    if strand == '+':
                        for temp_index in range(var_index + 1, var_index + len(ref)):
                            nucl_list[temp_index] = ''
                    else:
                        for temp_index in range(var_index - len(ref) + 1, var_index):
                            nucl_list[temp_index] = ''
                seq_mt = ''.join(nucl_list)

                AG_hit_wt = set(hit.end() for hit in re.finditer('AG', seq_wt))
                AG_hit_mt = set(hit.end() for hit in re.finditer('AG', seq_mt))
                AG_hit_set = AG_hit_mt - AG_hit_wt
                if 'del' in var_type:
                    shift = len(ref) - 1
                    AG_hit_set = {hit for hit in AG_hit_set if (hit + shift) not in AG_hit_wt}
                if 'ins' in var_type:
                    shift = -len(alt) + 1
                    AG_hit_set = {hit for hit in AG_hit_set if (hit + shift) not in AG_hit_wt}
                if AG_hit_set:
                    AG_hit = max(AG_hit_set)
                    if AG_hit >= 3 and seq_mt[AG_hit - 3] in 'CT':
                        DIVERS_conseq = 'RS_3SS_gain:'+str(AG_hit)+'nt*'
                        DIVERS_flag = 1
                        var_count_RS_3SS += 1
                        var_output_set.add(var_name)

            # RS_5SS_gain
            if RS_conseq.startswith('RS_5SS'):
                check_ref, check_alt = seq_wt.split('|')
                if ref == check_ref and alt == check_alt:
                    DIVERS_conseq = RS_conseq
                    DIVERS_flag = 1
                    var_count_RS_5SS += 1
                    var_output_set.add(var_name)

            # DIVERS_RS Output
            if DIVERS_flag:
                var_output_annot += 1
                file_out.write(sample+','+chrom+','+var_pos+','+var_id+','+ref+','+alt+','+strand+','+var_type+','+
                               RS_name_out+','+DIVERS_conseq+','+RS_info_out+','+var_info_out+'\n')
        file_map_RS.close()


        file_map_RSCRYP = open(filename_map_RSCRYP, 'r')
        for eachline in file_map_RSCRYP:
            item = eachline.strip().split('\t')
            chrom = item[0]
            var_name = item[3]
            strand = item[5]
            var_type = item[6]
            var_info = item[7]
            RSCRYP_name = item[11]
            RSCRYP_conseq = item[14]
            RSCRYP_type = item[15]

            _, var_pos, var_id, ref, alt = var_name.split('*')
            var_info_out = var_info.replace('*', ',')
            RSCRYP_name_out = RSCRYP_name.replace('|', ',')
            RSCRYP_info_out = '1,.,.,.,.,.,.,.,.,.,.,.'
            DIVERS_flag = 0
            DIVERS_conseq = ''

            # RSCRYP
            if var_name not in var_output_set:
                if RSCRYP_type=='NGGT':
                    if (strand=='+' and alt=='A') or (strand=='-' and alt=='T'):
                        DIVERS_flag = 1
                elif RSCRYP_type=='ANGT':
                    if (strand=='+' and alt=='G') or (strand=='-' and alt=='C'):
                        DIVERS_flag = 1
                elif RSCRYP_type=='AGNT':
                    if (strand=='+' and alt=='G') or (strand=='-' and alt=='C'):
                        DIVERS_flag = 1
                elif RSCRYP_type=='AGGN':
                    if (strand=='+' and alt=='T') or (strand=='-' and alt=='A'):
                        DIVERS_flag = 1

            # DIVERS_RSCRYP Output
            if DIVERS_flag:
                DIVERS_conseq = RSCRYP_conseq
                var_count_RS_CRYP += 1
                var_output_annot += 1
                var_output_set.add(var_name)
                file_out.write(sample+','+chrom+','+var_pos+','+var_id+','+ref+','+alt+','+strand+','+var_type+','+
                               RSCRYP_name_out+',.,'+DIVERS_conseq+','+RSCRYP_info_out+','+var_info_out+'\n')
        file_map_RSCRYP.close()

        os.system('rm ' + filename_bed)
        os.system('rm ' + filename_map_RS)
        os.system('rm ' + filename_map_RSCRYP)

        print('Sample:', sample)
        print('Input file:', filename_var)
        print('Output file:', filename_out, '\n')
        print('# Input variants:', str(len(var_input_set)))
        print('# Output variants:', str(len(var_output_set)))
        print('# Output annotations:', str(var_output_annot), '\n')
        print('# RS_AGGT_loss:\t', str(var_count_RS_AGGT))
        print('# RS_BP_loss:\t', str(var_count_RS_BP))
        print('# RS_3SS_gain:\t', str(var_count_RS_3SS))
        print('# RS_5SS_gain:\t', str(var_count_RS_5SS))
        print('# RS_CRYP:\t', str(var_count_RS_CRYP), '\n')

    except Exception as error_message:
        print(f"An error occurred: {error_message}")
        print("Please check your input files, or contact pzhang@rockefeller.edu.\n")

file_out.close()
