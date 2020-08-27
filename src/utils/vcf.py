#!/usr/bin/env python
# -*- coding: utf-8 -*-

import vcf
import pandas as pd
import numpy as np
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


def extract_attr(vcf_dict, vcf_info, attr, colname, sample=None, as_list = False):
    if type(sample) != type(None):
        vcf_dict = vcf_dict[sample]
        vcf_info = vcf_info[sample]
    
    if as_list:
        df = pd.DataFrame({colname:[vcf_dict[r].INFO[attr][0] for r in vcf_dict.keys()]})
    else:
        df = pd.DataFrame({colname:[vcf_dict[r].INFO[attr] for r in vcf_dict.keys() if attr in vcf_dict[r].INFO]})
    if type(sample) != type(None):
        df['Sample'] = sample
    
    return df


def accuracy_metrics(sample_vcf, ref_vcf, var_stat, common_stat, compare_to='Target'):
    '''
    Get accuracy metrics for calls
    :param sample: vcf filename for sample
    :param ref: vcf filename for reference
    :param compare_to: str to define if compare within Target or Reference regions
    '''
    if compare_to == 'Target':
        label = 'On_target'
    else:
        label = 'On_ref'

    ref_tot = var_stat.loc[var_stat.File.isin([ref_vcf]), label].tolist()[0]
    sample_tot = var_stat.loc[var_stat.File.isin([sample_vcf]), label].tolist()[0]
    TP = common_stat.loc[common_stat.Region.isin([compare_to]) & common_stat.Files.isin([f"{sample_vcf}, reference"]), 'Common_variants'].tolist()[0]
    FP = sample_tot-TP
    FN = ref_tot-TP
    
    precision=np.around((TP/(TP+FP))*100, 2)
    sensitivity=np.around((TP/(TP+FN))*100, 2)
    f1=2*((sensitivity*precision)/(sensitivity+precision))
    
    return precision, sensitivity, np.around(f1, 2)
