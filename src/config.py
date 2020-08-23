#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 2020.08.23
Author: Oksana Riba Grognuz

Contains project-level configurations
"""
import os

# Paths
basedir = os.path.abspath(os.path.dirname(__file__))
DATA_PATH = os.path.join(os.path.dirname(basedir), "data")
DATA_RAW_PATH = os.path.join(DATA_PATH, "raw")
DATA_INTERIM_PATH = os.path.join(DATA_PATH, "interim")
DATA_PROCESSED_PATH = os.path.join(DATA_PATH, "processed")

# STAT and METRICS FILES
SAM_STAT = "samtools.stats"
BAM_ALN = "bam_alignment.stats"
BAM_MET = "bam_metric.stats"
VARS_COMMON = "common_variants.stats"
VARS_REF = "reference_variant.stats"
VARS_TAR = "target_variant.stats" 
IGV = "igv_viewer.html"
SAMPLOTS = "samplots"
ANNOT_REG = "regions_annot.bed"
ANNOTR_TARGET = "targets_annot.bed"

# PRETTY NAMES
PRETTY_NAMES = {'SG001_1_GT.vcf.gz':'SG001_1',
                'SG001_1_GT_target.vcf':'SG001_1',
                'SG001_2_GT.vcf.gz':'SG001_2',
                'SG001_2_GT_target.vcf':'SG001_2',
                'SG001_ref_target.vcf': 'Reference'}

SNS_BARPLOT_STYLE = {'axes.facecolor': 'white',
                     'axes.edgecolor': '.8',
                     'axes.grid': False,
                     'axes.axisbelow': True,
                     'axes.labelcolor': '.15',
                     'figure.facecolor': 'white',
                     'grid.color': '.0',
                     'grid.linestyle': '-',
                     'text.color': '.15',
                     'xtick.color': '.15',
                     'ytick.color': '.15',
                     'xtick.direction': 'out',
                     'ytick.direction': 'out',
                     'lines.solid_capstyle': 'round',
                     'patch.edgecolor': 'w',
                     'patch.force_edgecolor': False,
                     'image.cmap': 'rocket',
                     'font.family': ['sans-serif'],
                     'font.sans-serif': ['Arial', 'DejaVu Sans', 'Liberation Sans',
                                         'Bitstream Vera Sans', 'sans-serif'],
                     'xtick.bottom': False,
                     'xtick.top': False,
                     'ytick.left': False,
                     'ytick.right': False,
                     'axes.spines.left': False,
                     'axes.spines.bottom': False,
                     'axes.spines.right': False,
                     'axes.spines.top': False}

SNS_BOXPLOT_STYLE = {'axes.facecolor': 'white',
                     'axes.edgecolor': '.8',
                     'axes.grid': False,
                     'axes.axisbelow': True,
                     'axes.labelcolor': '.15',
                     'figure.facecolor': 'white',
                     'grid.color': '.0',
                     'grid.linestyle': '-',
                     'text.color': '.15',
                     'xtick.color': '.3',
                     'ytick.color': '.4',
                     'xtick.direction': 'out',
                     'ytick.direction': 'out',
                     'lines.solid_capstyle': 'round',
                     'patch.edgecolor': 'w',
                     'patch.force_edgecolor': False,
                     'image.cmap': 'rocket',
                     'font.family': ['sans-serif'],
                     'font.sans-serif': ['Arial', 'DejaVu Sans', 'Liberation Sans',
                                         'Bitstream Vera Sans', 'sans-serif'],
                     'xtick.bottom': False,
                     'xtick.top': False,
                     'ytick.left': True,
                     'ytick.right': False,
                     'axes.spines.left': True,
                     'axes.spines.bottom': False,
                     'axes.spines.right': False,
                     'axes.spines.top': False}
