import vcf
import pandas as pd
import os

def read_vcf_data(vcf_file):
    '''
    Extracts dictionary of variants and samples from vcf file
    :param vcf_file: vcf file name (with full or relative path)
    '''
    vcf_reader = vcf.Reader(filename=vcf_file)
    vcf_dict = {}
    for record in vcf_reader:
        key = f"{record.CHROM}:{record.POS}"
        if key not in vcf_dict:
            vcf_dict[key]=record
        else:
            if isinstance(vcf_dict[key], list):
                vcf_dict[key].append(record)
            else:
                rec_list = [vcf_dict[key], record]
                vcf_dict[key]=rec_list

    return vcf_dict, vcf_reader.samples

def read_vcf_meta(vcf_file):
    '''
    Extracts meta data (info fields and filters) from vcf file
    :param vcf_file: vcf file name (with full or relative path)
    '''
    vcf_reader = vcf.Reader(filename=vcf_file)
    ids = [vcf_reader.infos[r].id for r in vcf_reader.infos.keys()]
    desc = [vcf_reader.infos[r].desc for r in vcf_reader.infos.keys()]
    dict_info = dict(zip(ids, desc))
    
    return dict_info, dict(vcf_reader.filters)


def extract_attr(var_dict, info_dict, attr, colname, sample):
    vcf_dict = var_dict[sample]
    vcf_info = info_dict[sample]
    df = pd.DataFrame({colname:[vcf_dict[r].INFO[attr][0] for r in vcf_dict.keys()]})
    df['Sample'] = sample
    return df

