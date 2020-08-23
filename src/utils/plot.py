import os
import sys
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.image as mpimg
from matplotlib_venn import venn2, venn2_circles
from matplotlib.cbook import boxplot_stats
import numpy as np
import seaborn as sns
import pandas as pd
import warnings
from src import config
from src.config import SNS_BARPLOT_STYLE, SNS_BOXPLOT_STYLE

sns.set_style('whitegrid', {'axes.grid': False})
sns.set(style='ticks')
pd.set_option('max_colwidth', 100)
warnings.filterwarnings("ignore")


def show_samplots_image(png_name):
    image = mpimg.imread(os.path.join(config.DATA_PROCESSED_PATH, config.SAMPLOTS, png_name))
    plt.figure(figsize=(14,10))
    plt.imshow(image)
    plt.axis('off')
    plt.show()


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n}, {a:.2f}, {b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


def plot_var_venn(filename, var_stat, common_stat, compare_to='Target'):
    '''
    Plot venn diagram for sam ple vs reference
    :param filename: vcf filename for sample
    :param compare_to: str to define if compare within Target or Reference regions
    '''
    if compare_to == 'Target':
        label = 'On_target'
    else:
        label = 'On_ref'
    
    cmap = plt.get_cmap('Blues')
    cmap_tr = truncate_colormap(cmap, 0.3, 1)

    ref_tot = var_stat.loc[var_stat.File.isin(['SG001_ref.vcf.gz']), label].tolist()[0]
    sample_tot = var_stat.loc[var_stat.File.isin([filename]), label].tolist()[0]
    common = common_stat.loc[common_stat.Region.isin([compare_to]) & common_stat.Files.isin([filename]), 'Common_variants'].tolist()[0]
    subsets = (sample_tot-common, ref_tot-common, common)

    venn2(subsets=subsets, set_labels = (config.PRETTY_NAMES[filename], 'Reference'), set_colors = (cmap_tr(10), cmap_tr(30), cmap_tr(50)), alpha=0.75)
    venn2_circles(subsets=subsets, color='black', alpha=0.7, linewidth=0.4)
    plt.title(f"{compare_to} regions")
    plt.show()
    plt.close()


def show_values_on_bars(axs, h_v="v", space=0.4, fontsize=12, color='black', symbol='', fontweight='bold'):
    def _show_on_single_plot(ax):
        if h_v == "v":
            for p in ax.patches:
                _x = p.get_x() + p.get_width() / 2
                _y = p.get_y() + p.get_height() + float(space)
                value = f"{p.get_height():.1f}{symbol}"
                ax.text(_x, _y, value, ha="center", fontsize=fontsize, color=color, fontweight=fontweight)
        elif h_v == "h":
            for p in ax.patches:
                _x = p.get_x() + p.get_width() + float(space)
                _y = p.get_y() + p.get_height()
                value = f"{p.get_width():.1f}{symbol}"
                ax.text(_x, _y, value, ha="left", fontsize=fontsize, color=color, fontweight=fontweight)

    if isinstance(axs, np.ndarray):
        for idx, ax in np.ndenumerate(axs):
            _show_on_single_plot(ax)
    else:
        _show_on_single_plot(axs)


def plot_nice_barplot(df, x_col, y_col, hue_col, title_text, x_lab, y_lab, hue_lab,
                      space=-9, fontsize=29, symbol='', figsize=(10, 6)):
    """
    Plot barplot
    :param df:
    :param x_col:
    :param y_col:
    :param hue_col:
    :param title_text:
    :param x_lab:
    :param y_lab:
    :param hue_lab:
    :param space:
    :param fontsize:
    :param symbol:
    :param figsize:
    :return: barplot
    """
    sns.set_style(SNS_BARPLOT_STYLE)
    plt.figure(figsize=figsize)
    ax = sns.barplot(x=x_col, y=y_col, hue=hue_col, data=df, palette="Blues_r")
    ax.set_title(title_text, fontsize=20, color='#505050')
    ax.set_ylabel(y_lab, fontsize='15', color='#A9A9A9', fontweight='bold')
    ax.set_xlabel(x_lab, fontsize='15', color='#A9A9A9', fontweight='bold')
    show_values_on_bars(ax, space=space, fontsize=fontsize, color='white', symbol=symbol)
    # ax.axes.get_yaxis().set_visible(False)
    ylabels = plt.yticks()
    ax.set_yticklabels(ylabels, color='white', size=1)
    _, xlabels = plt.xticks()
    ax.set_xticklabels(xlabels, size=14, color='#505050')
    legend = ax.get_legend()
    legend.set_title(hue_lab)
    legend.get_frame().set_linewidth(0.0)
    plt.setp(legend.get_texts(), fontsize='14', color='#505050')
    plt.setp(legend.get_title(), fontsize='15', color='#A9A9A9', fontweight='bold')
    # plt.legend(fontsize='large', title_fontsize='20')
    plt.show()
    plt.close()


def make_boxplot(df_data, x_col, y_cols, main_title='', x_lab='', y_labs='', refline=None):
    """
    Produces boxplots
    :param df_data: dataframe with data to plot
    :param x_col: (string) column in df_data to plot on x axis
    :param y_cols: (list of strings) columns in df_data to plot on y axis
    :param refline: (list of floats) reference line on y axis
    :param main_title: main plot title (str)
    :param x_lab: x axis label (str)
    :param y_labs: list if strings with y axis labels
    """
    sns.set_style(SNS_BOXPLOT_STYLE)
    n_subplots = len(y_cols)
    fig_sizes = {1: (6, 4), 2: (11, 4), 3: (16, 4)}
    f, axes = plt.subplots(1, n_subplots, figsize=fig_sizes[n_subplots])
    f.suptitle(main_title, fontsize=14, color='#505050', fontweight='bold')
    f.subplots_adjust(top=0.85)
    meanprops = {"marker": "s", "markerfacecolor": "white", "markeredgecolor": "gray"}

    for i in range(0, n_subplots):
        col = y_cols[i]
        ax = sns.boxplot(x=x_col, y=col, data=df_data, showmeans=True, meanprops=meanprops,
                         palette="Blues_r", width=0.3, linewidth=0.5, orient='v', ax=axes[i])
        if len(x_lab) > 0:
            ax.set_xlabel(x_lab, fontsize=12, color='#A9A9A9', fontweight='bold')
        if len(y_labs) > 0:
            ax.set_ylabel(y_labs[i], fontsize=12, color='#A9A9A9', fontweight='bold')
        if refline is not None:
            ax.axhline(refline[i], ls='--', color='#5e6f80')
            max_ref = max(df_data[col]) - refline[i]
            ax.text(0.23, refline[i] + max_ref * 0.05, f"Global: {refline[i]:.2f}", color='#5e6f80')
    plt.show()
    plt.close()

