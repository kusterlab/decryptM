# pipeline.py
#
# This script is build for simple MQ evidence based decryptM data-analysis.
#
# Run this script as follows:
# $python pipeline.py -q -p <./config.toml>
#
# -q : activate quantification script
# -p : activate plotting script
# <./config.toml> toml file with all relevant parameters to run the script
#
# Florian P. Bayer - Aug. 2022
#
# version 1.0
#

# ####### #
# Imports #
# ####### #
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import sys
import os
import argparse
import toml

import numpy as np
import pandas as pd
from scipy import optimize, interpolate

# If on linux server use agg backend
import matplotlib
if sys.platform.startswith('linux'):
    matplotlib.use('Agg')
from matplotlib import pyplot as plt
plt.rcParams['text.usetex'] = False
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

import tqdm
from tqdm.autonotebook import tqdm as autonotebook
autonotebook.pandas()

import filters as ddPTM_filters


# ######### #
# Functions #
# ######### #
def check_path(path, create_folder=False):
    """
    Check the given path whether it exists and is accessible.

    if create_folder is True and folder dose not exist yet it will be created
    """
    if not os.path.exists(path):
        if create_folder:
            os.mkdir(path)
        else:
            raise FileNotFoundError(f'Path "{path}" does not exist.')
    if not os.access(path, os.R_OK):
        raise PermissionError(f'File at "{path}" cannot be opened. Try to close it elsewhere.')


def check_toml(config):
    """
    Check the given script input on validity and report values that are not meeting the specs to the user

    Parameters
    ----------
    config : toml object
        The config file with the parameters imported (dict)
    """
    # Report Meta data
    print('Data set: ', config['Meta Data']['dataset_name'])
    print('Description: ', config['Meta Data']['experiment_description'])

    # Check TMT plexing
    if int(config['TMT']['plex']) != len(config['TMT']['channels']):
        warnings.warn('Warning: plex number and channel number dose not equal.')

    # Check experiment type specific input
    if config['Meta Data']['experiment_type'] == 'dd':
        if len(config['TMT']['channels']) != len(config['TMT']['doses']):
            raise ValueError("Channels and Doses do no correspond")
    elif config['Meta Data']['experiment_type'] == 'td':
        if len(config['TMT']['channels']) != len(config['TMT']['time_points']):
            raise ValueError("Channels and time points do no correspond")
    else:
        raise ValueError(f"Unknown experiment type {config['Meta Data']['experiment_type']}.")

    # Report Processing
    print('Control Channel: ', config['Processing']['control_channel'])
    print('Normalization based on: ', config['Processing']['subproteome_normalization'])


def clean_modified_sequence(mod_seq, remove_mox=True):
    """
    clean_modified_sequence(mod_seq, remove_mox=True)

    Cleans up the modified sequence. Columns-wise operation is 10x faster than apply.
    Try to add different MaxQuant nomenclatures.

    Can clean away mox if wanted.

    Parameters
    ----------
    mod_seq : pd.Series(<seqs>)
    remove_mox : boolean, default True

    Returns
    -------
    mod_seq : pd.Series(<seqs>)
    """
    mod_seq = mod_seq.str.replace('_', '')
    mod_seq = mod_seq.str.replace('\(Acetyl \(Protein N-term\)\)', '(ac)')
    mod_seq = mod_seq.str.replace('\(Oxidation \(M\)\)', '(ox)')

    mod_seq = mod_seq.str.replace('pS', 'S(ph)').str.replace('pT', 'T(ph)').str.replace('pY', 'Y(ph)')
    mod_seq = mod_seq.str.replace('\(Phospho \(STY\)\)', '(ph)')
    mod_seq = mod_seq.str.replace('\(Acetyl \(K\)\)', '(ac)')
    mod_seq = mod_seq.str.replace('\(GG \(K\)\)', '(ub)')

    if remove_mox:
        mod_seq = mod_seq.str.replace('M\(ox\)', 'M')

    return mod_seq


def load_cleaned_evidence(path, sum_cols=[], first_cols=[], max_cols=[], min_cols=[], concat_cols=[],
                          merge_duplicates=True, remove_contaminants=True, remove_reverse=True):
    """
    load_cleaned_evidence(path, sum_cols=[], first_cols=[], max_cols=[], min_cols=[], concat_cols=[],
                          merge_duplicates=True, remove_contaminants=True, remove_reverse=True)

    This function loads the evidence.txt from MaxQuant and cleans it up as specified.
    The modified sequence and the experiment name are used for unique grouping.

    Parameters
    ----------
    path : str
        path to the file
    sum_cols : array-like
        column names where duplicated entries will be summed up (cols of dtype numeric)
    first_cols : array-like, optional
        column names where the first row will be used for duplicated entries (cols of dtype categorical)
    max_cols : array-like, optional
        column names where the max value will be used for duplicated entries (cols of dtype numeric)
    min_cols : array-like, optional
        column names where the min value will be used for duplicated entries (cols of dtype numeric)
    concat_cols : array-like, optional
        column names where all values will be concatenated for duplicated entries (cols of dtype object)
    merge_duplicates : bool, optional
        If duplicated entries should be merged to a single entry. The merge action depends on the column groups above.
        Default value is True.
    remove_contaminants : bool, optional
        If contaminant peptides should be removed the analysis. Contaminant entries come from MaxQuants contaminant list
        and are indicated by a '+'. Default value is True.
    remove_reverse : bool, optional
        If reversed peptides should be removed the analysis. Reversed entries come from MaxQuants FDR decoy approach
        and are indicated by a '+'. Default value is True.

    Returns
    -------
    evd : pd.DataFrame -> the cleaned MaxQuant output
    """
    # Load
    evd = pd.read_csv(path, sep='\t', low_memory=False)

    # Different MQ version use different Mod. sequence notations... Rewrite the PTM positions consistently
    evd['Modified sequence'] = clean_modified_sequence(evd['Modified sequence'])

    # Remove contaminants & decoys
    if remove_contaminants:
        evd = evd[evd['Potential contaminant'].replace(np.nan, '') != '+']
    if remove_reverse:
        evd = evd[evd['Reverse'].replace(np.nan, '') != '+']

    # Remove peptide duplicates and solve duplicated entries issue by sum, first, min, max
    if merge_duplicates:
        if not sum_cols:
            raise ValueError("Specify a list of columns to sum when using the sum_duplicates option")

        grouped_evd = evd.groupby(['Modified sequence', 'Experiment'])

        # Count & Sum the duplicates, then merge results to new_evd
        evd_count = grouped_evd.size().to_frame().rename({0: 'N duplicates'}, axis=1)
        evd_sum = grouped_evd[sum_cols].sum()
        new_evd = pd.merge(left=evd_count, right=evd_sum, left_index=True, right_index=True)

        # do further aggregations on other columns and merge to new_evd
        if first_cols:
            evd_first = grouped_evd[first_cols].first()
            new_evd = pd.merge(left=new_evd, right=evd_first, left_index=True, right_index=True)

        if max_cols:
            evd_max = grouped_evd[max_cols].max()
            new_evd = pd.merge(left=new_evd, right=evd_max, left_index=True, right_index=True)
            new_max_cols = list(map(lambda s: 'Max ' + s, max_cols))
            new_evd.rename({old: new for old, new in zip(max_cols, new_max_cols)}, axis=1, inplace=True)
            max_cols = new_max_cols

        if min_cols:
            evd_min = grouped_evd[min_cols].min()
            new_evd = pd.merge(left=new_evd, right=evd_min, left_index=True, right_index=True)
            new_min_cols = list(map(lambda s: 'Min ' + s, min_cols))
            new_evd.rename({old: new for old, new in zip(min_cols, new_min_cols)}, axis=1, inplace=True)
            min_cols = new_min_cols

        if concat_cols:
            evd[concat_cols] = evd[concat_cols].replace(np.nan, '')

            for c_col in concat_cols:
                evd_concat = grouped_evd[c_col].apply(';'.join)
                new_evd = pd.merge(left=new_evd, right=evd_concat, left_index=True, right_index=True)

            new_concat_cols = list(map(lambda s: 'All ' + s, concat_cols))
            new_evd.rename({old: new for old, new in zip(concat_cols, new_concat_cols)}, axis=1, inplace=True)
            concat_cols = new_concat_cols

        # Resort the columns
        evd = new_evd[['N duplicates'] + first_cols + concat_cols + max_cols + min_cols + sum_cols]

    evd.reset_index(inplace=True)
    return evd


def build_reporter_col_names(col_name, id_list):
    """
    Build the names for the reporter columns
    """
    return [f'{col_name} {int(i)}' for i in id_list]


def add_ptm_mask_column(df, col_name, str_symbol):
    """
    Adds a mask column to the input df, that checks if a modification is present is present in the 'Modified sequence'.
    """
    df[col_name] = df['Modified sequence'].str.count(str_symbol) > 0
    return df


def add_fullproteome_mask_column(df, ptms):
    """
    Add a full proteome mask. Everything that is not a ptm is considered a full proteome peptide.
    """
    ptm_cols = [col_name for (col_name, ptm_name, ptm_symbol) in ptms.values()]
    df['Fullproteome'] = (df[ptm_cols].sum(axis=1) == 0)
    return df


def get_expected_ptms(config):
    """
    Parse the expected ptms from the config file and get the respective lookup table.

    Parameters
    ----------
    config : toml object
        The config file with the parameters

    Returns
    -------
    ptms : dict(<PTM>: tuple(<ptm_proteome, mq_ptm_column_name, ptm_str_symbol>))
        A dictionary containing the ptm as key and a tuple with ptm relevant information.
        pos0: Name of the sub-proteome, pos1: MQ evidence column name, pos2: modified sequence symbol
    """
    ptms = {}
    if config['PTMs'].get('Phosphorylation'):
        ptms['Phosphorylation'] = ('Phosphoproteome', 'Phospho (STY)', '\(ph\)')
    if config['PTMs'].get('Ubiquitinylation'):
        ptms['Ubiquitinylation'] = ('Ubiproteome', 'GG (K)', 'K\(ub\)')
    if config['PTMs'].get('Acetylation'):
        ptms['Acetylation'] = ('Acetylproteome', 'Acetyl (K)', 'K\(ac\)')
    return ptms


def build_drug_log_concentrations(steps, scale, dmso_offset=1000):
    """
    build_drug_concentrations(steps, scale)

    This function takes a sequence of drug steps and its corresponding scale and calculates the respective
    concentration numpy array. The DMSO (0) will be corrected to a n times smaller number than the smallest value.
    This is necessary for log-step transformations later on. Default shift is 1000.

    Parameters
    ----------
    steps : array-like
        The drug steps in real space. Order does not matter
    scale : numeric
        The drug scaling factor
    dmso_offset : numeric
        The dmso offset to display the dmso in log space. By default 1000 (=3 orders of magnitude)

    Returns
    -------
    log_steps : array-like
        The drug log steps in the original order
    """
    # Return empty array if empty input
    if len(steps) == 0:
        return np.array([])
    # Else calculate the log steps with defined dmso offset based on input
    steps = np.array(steps)
    steps[steps == 0] = min(steps[steps != 0]) / dmso_offset
    log_steps = np.log10(steps * scale)
    return log_steps


def add_median_normalization(df, raw_cols, norm_cols, ref_col=''):
    """
    add_median_normalization(df, raw_cols, norm_cols, ref_col)

    Median centric normalization of raw_cols using the rows of the reference only.
    The normalized values will be added to the df under the name of norm_cols.

    Parameters
    ----------
    df : pd.DataFrame
        A data frame containing at least raw_cols, ref_col
    raw_cols : array-like
        A array-like object containing the column labels of the raw data
    norm_cols : array-like
        A array-like object containing the column labels of the future normalized data
    ref_col : string
        A string indicating a column with booleans which is used as a reference for normalization
        By default, '', meaning that all rows are used to calculate the normalization factors.

    Returns
    -------
    df : pd.DataFrame
        Output data frame with the new normalized columns
    """
    mask = df[ref_col] if ref_col else df.index
    ref_medians = np.log2(df.loc[mask, raw_cols].replace(0, np.nan)).median()
    normalization_factors = ref_medians.mean() - ref_medians
    df[norm_cols] = 2 ** (np.log2(df[raw_cols].replace(0, np.nan)) + normalization_factors)
    df[norm_cols] = df[norm_cols].replace([np.nan, -np.inf, np.inf], 0)
    return df


def calculate_median_normalization_factors(df, cols, ref_col):
    """
    This function is for reporting purpose and equivalent to add_median_normalization
    """
    mask = df[ref_col] if ref_col else df.index
    ref_medians = np.log2(df.loc[mask, cols].replace(0, np.nan)).median()
    normalization_factors = ref_medians.mean() - ref_medians
    normalization_factors['Reference column'] = ref_col
    return normalization_factors


def add_ratios(df, cols, ratio_cols, ref_col):
    """
    add_ratios(df, cols, ratio_cols, ref_col)

    Calculate ratios of cols / ref_col
    The ratio values will be added to the df under the name of ratio_cols.

    Parameters
    ----------
    df : pd.DataFrame
        A data frame containing at least raw_cols, ref_col
    cols : array-like
        A array-like object containing the column labels of the data
    ratio_cols : array-like
        A array-like object containing the column labels of the future ratio data
    ref_col : string
        A string indicating a column used as a reference for ratio calculations.

    Return
    ------
    df : pd.DataFrame
        Output data frame with the new ratio columns
    """
    df[ratio_cols] = df[cols].div(df[ref_col], axis=0).replace([np.inf], np.nan)
    return df


def logisitc_model(x, ec50, slope, top, bottom):
    """
    Logistic model to fit the drug response data.

    Parameters
    ----------
    x: array-like
        Drug concentrations in log space
    ec50: float
        Inflection point in log space of the sigmoid value
    slope: float
        slope of the transition between top and bottom
    top: float
        top asymptote
    bottom: float
        bottom asymptote

    Returns
    -------
    y: array-like
        Drug response
    """
    # predict y with given parameters using the standard model
    return (top - bottom) / (1 + 10 ** (slope * (x - ec50))) + bottom


def fit_logistic_function(y_data, x_data, x_interpolated=None, curve_guess=None, curve_bounds=None, max_opt=None):
    """
    fit_logistic_function(y_data, x_data, x_interpolated=None, max_opt=None)

    This function is fitting a 4 parameter logistic function to the given y & x data.

    Parameters
    ----------
    y_data : pd.Series
        A Series object, containing the observed DMSO normalized ratios (including the DMSO ratio)
    x_data : array_like
        An array, containing the used drug concentrations in Log space
    x_interpolated : array_like, optional
        If None (default), the fit will be performed on the pure y&x data. If object is an array,
        containing drug concentrations, then the logistic fit is performed on those x_interpolated
        values. To obtain y_interpolated values that are corresponding to the x_interpolated, a
        linear interpolation is performed. The values of x_interpolated should be within within the
        x_data range to not extrapolate the data.
    curve_guess : array_like, optional
        If None (default), the fit will use no initial parameter guesses.
    curve_bounds : tuple of array_like, optional
        If None (default), the fit will use no bounds for the parameter values.
        Else provide the following structure ([min1, ..], [max1, ..]).
    max_opt : int, optional
        If None (default), the fit will use scipy default estimation of a maximum number of fit
        iterations before the fit is terminated. Otherwise, provide a integer to specify the maximum
        number yourself to save time or increase number of successful fits in longer time.

    Returns
    -------
    out : tuple
        An tuple object containing the fitted parameter, R2, and the parameter standard errors in
        the following order: ('Log EC50', 'Curve slope', 'Curve top', 'Curve bottom', 'R2', 'Curve RMSE',
        'Log EC50 error', 'Curve slope error', 'Curve top error', 'Curve bottom error')
    """
    # fast return
    if any(np.isnan(y_data)):
        return 10 * (np.nan,)

    # Fit logistic model
    try:
        y_data = y_data.values

        if x_interpolated is not None:
            f = interpolate.interp1d(x_data, y_data, kind='linear')
            popt, pcov = optimize.curve_fit(logisitc_model, x_interpolated, f(x_interpolated),
                                            p0=curve_guess, bounds=curve_bounds, max_nfev=max_opt)
            perr = np.sqrt(np.diag(pcov))

        else:
            popt, pcov = optimize.curve_fit(logisitc_model, x_data, y_data, p0=curve_guess,
                                            bounds=curve_bounds, max_nfev=max_opt)
            perr = np.sqrt(np.diag(pcov))

        # Calculate R2 and RMSE
        ss_res = np.sum((y_data - logisitc_model(x_data, *popt)) ** 2)
        ss_tot = np.sum((y_data - np.mean(y_data)) ** 2)
        r2 = 1 - ((ss_res + 1e-10) / (ss_tot + 1e-10))
        rmse = np.sqrt(ss_res/len(y_data))

        return (*popt, r2 if r2 >= 0 else 0, rmse, *perr)

    except RuntimeError:
        return 10 * (np.nan,)


def build_interpolation_points(x, exclude_low_n=1, exclude_top_n=0, interpolation_size=np.log10(10 / 3) / 2):
    """
    build_interpolation_points(x, exclude_low_n=1, exclude_top_n=0)

    Constructs an interpolated array, where the new elements lay with space of maximum interpolation_size between
    the x elements array

    Parameters
    ----------
    x : array-like
        Input array
    exclude_low_n : int
        The number of elements to exclude from the interpolation at the lower end. By default 1 (exclude DMSO)
    exclude_top_n : int
        The number of elements to exclude from the interpolation at the higher end. By default 0
    interpolation_size : float
        The maximum interpolation size before distance between x is halved.
        Default is np.log10(10/3)/2 for half log10 step interpolation

    Returns
    -------
    x_interpolated : array-like
        An array with added interpolation points that are equidistant smaller than interpolation_size
    """
    while_counter = 0
    while True:
        # sort x and calculate the x difference to
        x = sorted(x)
        x_diff = np.diff(x[exclude_low_n:len(x) - exclude_top_n])

        # if all smaller then full log step
        if all(x_diff <= interpolation_size):
            return x
        # else: half the distance where it is too large
        half_x = list((x[pos + exclude_low_n] + x[pos + exclude_low_n + 1]) / 2 for pos in np.where(x_diff >= interpolation_size)[0])
        x = np.append(x, half_x)

        # while loop trap protection
        while_counter += 1
        if while_counter > 10:
            raise ValueError('Wrong interpolation size triggered while-loop escape.')


def add_logistic_model(df, ratio_cols, x_data, interpolation=True, max_opt=None,
                       interpolation_kwargs={'exclude_low_n': 1, 'exclude_top_n': 0}):
    """
    add_logistic_model(ddf, ratio_cols, x_data, interpolation=True, max_opt=None,
                       interpolation_kwargs={'exclude_low_n': 1, 'exclude_top_n': 0})

    Fits a logistic model to the ratio columns. This is a wrapper function for "fit_logistic_function(...)".

    Parameters
    ----------
    df : pd.DataFrame
        A data frame containing at least ratio_cols
    ratio_cols : array_like
        A array-like object containing the column labels of the ratio data
    x_data : array_like
        A array-like object containing the drug concentrations in log space.
    interpolation: boolean, optional
        If True (default), the interpolated approach is used for the model fit. Else not.
    max_opt : int, optional
        If None (default), the fit will use scipy default estimation of a maximum number of fit
        iterations before the fit is terminated. Otherwise, provide a integer to specify the maximum
        number yourself to save time or increase number of successful fits in longer time.
    interpolation_kwargs : dict, optional
        parameters passed over to the interpolation point function. By default the dmso control (lowest dose) is
        excluded therefore exclude 1 from the lowest. All other doses should be interpolated therefore exclude 0 (none).

    Returns
    -------
    df : pd.DataFrame
        Output data frame with the fitted columns
    """
    # Calculate the interpolation schema if required. Interpolation only in between drug concentrations.
    x_interpolated = build_interpolation_points(x_data, **interpolation_kwargs) if interpolation else None

    # Fitting parameter [EC50, slope, top, bottom]
    curve_bounds = ([min(x_data) - 3, 0, 0, 0], [max(x_data) + 3, 25, 10, 200])
    curve_guess = [np.median(x_data), 1, 1, 0]

    # Fit the ratio data  to the logistic function
    fits = df[ratio_cols].progress_apply(fit_logistic_function,
                                         x_data=x_data,
                                         x_interpolated=x_interpolated,
                                         curve_bounds=curve_bounds,
                                         curve_guess=curve_guess,
                                         max_opt=max_opt,
                                         axis=1)

    # Transform and add data
    fit_cols = ['Log EC50', 'Curve slope', 'Curve top', 'Curve bottom', 'R2', 'Curve RMSE',
                'Log EC50 error', 'Curve slope error', 'Curve top error', 'Curve bottom error']
    fits = pd.DataFrame.from_records(data=fits.apply(np.array), columns=fit_cols)
    df[fit_cols] = fits

    # Add more curve columns
    df['EC50'] = 10 ** (df['Log EC50'])
    df['pEC50'] = -1 * df['Log EC50']
    df['Curve effect size'] = df['Curve bottom'] - df['Curve top']
    return df


def define_regulated_curves(df, drug_concs, col_names, min_signal, data_type='ms3'):
    """
    define_regulated_curves(df, first_dose_col, high_dose_col, data_type='ms3')

    Apply the curve filters and add a new column <Regulation> to the input df. This function defines what is considered
    a regulation and what not. The definitions for each quantification approach and direction are stored in the
    filters.py next to this script.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame
    drug_concs: array_like
        an array containing the drug concentrations in real space
    col_names: array_like
        an array containing the column names with the ratio channels
    min_signal: float
        the minimum DMSO signal required for a curve to have good quantification (machine-specific).
    data_type : str, optional
        the measurement data type. Default <ms3>, alternatively <ms2>
    Return
    ------
    df : pd.DataFrame
    """
    # Select the right filter
    f_up = ddPTM_filters.ms2_up if data_type == 'ms2' else ddPTM_filters.ms3_up
    f_down = ddPTM_filters.ms2_down if data_type == 'ms2' else ddPTM_filters.ms3_down

    # Extract the name of the lowest and highest dose column
    drug_concs = np.array(drug_concs)
    first_dose = np.nanmin(drug_concs[np.nonzero(drug_concs)])
    first_dose_idx = np.argwhere(drug_concs == first_dose)[0, 0]
    high_dose_idx = np.nanargmax(drug_concs)
    first_dose_col = col_names[first_dose_idx]
    high_dose_col = col_names[high_dose_idx]

    # Apply the up filter
    curve_filter_up = (
          (df['R2'] >= f_up['R2 MINIMUM'])
        & (df['Curve RMSE'] <= (f_up['RMSE MAXIMUM'](df[col_names].mean(axis=1))))
        & (df['Curve top'] <= (f_up['TOP MAXIMUM'](df['Curve effect size'])))
        & (df['Curve bottom'] >= f_up['BOTTOM MINIMUM'])
        & (df['Curve effect size'] >= f_up['EFFECT SIZE MINIMUM'])
        & (df['Curve slope'] <= f_up['SLOPE MAXIMUM'])
        & (df[high_dose_col] >= f_up['HIGH DOSE MINIMUM'])
        & (df[first_dose_col] <= (f_up['FIRST DOSE MAXIMUM'](df['Curve effect size'])))
    )

    # Apply the down filter
    curve_filter_down = (
          (df['R2'] >= f_down['R2 MINIMUM'])
        & (df['Curve RMSE'] <= f_down['RMSE MAXIMUM'])
        & (df['Curve top'].between(*f_down['TOP RANGE']))
        & (df['Curve bottom'] <= f_down['BOTTOM MAXIMUM'])
        & (df['Curve effect size'] <= f_down['EFFECT SIZE MAXIMUM'])
        & (df['Curve slope'] <= f_down['SLOPE MAXIMUM'])
        & (df[high_dose_col] <= f_down['HIGH DOSE MAXIMUM'])
        & (df[first_dose_col].between(*f_down['FIRST DOSE RANGE']))
    )

    # Add labels
    df['Regulation'] = np.nan
    df.loc[curve_filter_down, 'Regulation'] = 'down'
    df.loc[curve_filter_up, 'Regulation'] = 'up'

    # Apply the min_signal filter by removing potential regulations that are below the min signal threshold.
    df.loc[df['Curve signal'] < min_signal, 'Regulation'] = np.nan

    return df


def get_sorting_channel(col_names, experiment_type, drug_cons, time_points):
    """
    Get the chanel name that should be used for sorting the responses.
    This will be the max dose or time depending on the experiment.
    """
    if experiment_type == 'dd':
        return col_names[np.nanargmax(drug_cons)]
    if experiment_type == 'td':
        return col_names[np.nanargmax(time_points)]
    raise ValueError('Experiment type is unknown. Set it to td (time) or dd (dose)')


def sort_data(df, channel):
    """
    Sort the data frame according to the given column channel.
    The column will be log2 transformed and the sign gets removed. Large values to top.
    """
    min_value = df[channel].replace(0, np.nan).min()
    mask = np.abs(np.log2(df[channel].replace(0, min_value))).sort_values(ascending=False).index
    df = df.loc[mask].reset_index(drop=True)
    return df


def sort_td_data(df, cols):
    """
    Sort the time-dependent data frame according to the given intensity columns.
    The max absolute difference will be sorted to the top
    """
    max_vals = df[cols].replace(0, np.nan).apply(np.log2).apply(np.abs).max(axis=1)
    mask = max_vals.sort_values(ascending=False).index
    df = df.loc[mask].reset_index(drop=True)
    return df


def add_intensity_dist(ax, dist_values, cmap, colorrange, signal):
    """
    Add an intensity distribution to the figure axis.
    """
    color = cmap(colorrange(signal))
    sns.distplot(dist_values, vertical=True, ax=ax, hist_kws={'color': 'gray'}, kde_kws={'color': 'black'})

    ax.axhline(signal, color=color)
    ax.set_ylabel('Control Channel Log2 Intensity Distribution')
    ax.set_xticks([])
    ax.set_xticklabels([])


def add_plot_title(ax, gene_name, mod_seq):
    """
    Add the plot title to the figure axis.
    """
    fig_title = '{} @ {}'
    gene_name, mod_seq = str(gene_name), str(mod_seq)

    gene_name = gene_name if len(gene_name) < 30 else (gene_name[:28] + '..')
    mod_seq = mod_seq if len(mod_seq) < 55 else (mod_seq[:53] + '..')

    ax.set_title(fig_title.format(mod_seq, gene_name), fontsize=10, fontweight='bold')


def add_dd_plot(ax, curve, cmap, colorrange, tmt_cols_ratio, drug_log_concs, drug_unit, curve_x_extrapolate,
                plot_helper_line, plot_fit=True):
    """
    add_dd_plot(ax, curve, cmap, colorrange, tmt_cols_ratio, drug_log_concs, drug_unit, curve_x_extrapolate,
                plot_helper_line, plot_fit=True

    Add the logistic dose-dependent plot to the given figure axis.

    Parameters
    ----------
    ax : matplotlib axis object
    curve : pd.Series with the data, curve fit values, and meta data
    cmap : matplotlib color map object
    colorrange : matplotlib color range object
    tmt_cols_ratio : array-like, with the tmt column names
    drug_log_concs : array-like, with the drug concentrations in log10 space
    drug_unit : str, the unit of the x-axis.
    curve_x_extrapolate : array-like, with left and right extrapolation log distances of the fitted curved.
    plot_helper_line : array-like, with the y values of the lower and upper boundary. Will be drawn as dotted lines.
    plot_fit : bool, optional, default True, flag if the fitted curve should be plotted.

    Return
    ------
    None
    """
    # Signal Noise based coloring
    color_sn = cmap(colorrange(curve['Curve signal']))

    # experimental data
    x_data = drug_log_concs
    y_data = curve[tmt_cols_ratio]

    # logistic fit
    if plot_fit:
        ax.scatter(10 ** x_data, y_data, color=color_sn, edgecolors='none', linewidth=0)
        params = curve[['Log EC50', 'Curve slope', 'Curve top', 'Curve bottom']]
        x_curve = np.arange(min(x_data) + curve_x_extrapolate[0], max(x_data) + curve_x_extrapolate[1], 0.01)
        ax.plot(10 ** x_curve, logisitc_model(x_curve, *params), color=color_sn)

        # Parameter info box
        infobox = 'LogSignal = {:.1f}\nEC50 = {:.2e}\nSlope = {:.2f}\nTop = {:.2f}\nBottom = {:.2f}\nEffekt = {:.2f}\nR² = {:.3f}'
        infobox = infobox.format(curve['Curve signal'], curve['EC50'], *params.values[1:],
                                 curve['Curve effect size'], curve['R2'])
        infobox_props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(0.01, 0.98, infobox, transform=ax.transAxes,
                fontsize=11, verticalalignment='top', bbox=infobox_props)
    else:
        ax.plot(10**x_data, y_data, '-o', color=color_sn)

    # add the title
    add_plot_title(ax, curve['Gene names'], curve['Modified sequence'])

    # axes make the plot pretty
    ax.set_xscale('log')
    ax.set_xlim((10 ** (min(x_data) + curve_x_extrapolate[0] - 0.5),
                 10 ** (max(x_data) + curve_x_extrapolate[1] + 0.5)))
    ax.set_xlabel(f'Drug Concentration [{drug_unit}]')
    ax.set_ylabel('Relative Response')
    ax.set_ylim((0, None if max(y_data) > 1.9 else 2))
    ax.axhline(y=plot_helper_line[0], ls=':', color='black')
    ax.axhline(y=plot_helper_line[1], ls=':', color='black')


def add_td_plot(ax, curve, cmap, colorrange, tmt_cols_ratio, time_points, time_unit, plot_helper_line):
    """
    add_td_plot(ax, curve, cmap, colorrange, tmt_cols_ratio, time_points, plot_helper_line)

    Add the time-dependent plot to the given figure axis. This is just a linear connection between the observed
    data points.

    Parameters
    ----------
    ax : matplotlib axis object
    curve : pd.Series with the data, curve fit values, and meta data
    cmap : matplotlib color map object
    colorrange : matplotlib color range object
    tmt_cols_ratio : array-like, with the tmt column names
    time_points : array-like, with the drug concentrations in log10 space
    time_unit : str, the unit of the x-axis
    plot_helper_line : array-like, with the y values of the lower and upper boundary. Will be drawn as dotted lines.

    Return
    ------
    None
    """
    # Signal Noise based coloring
    color_sn = cmap(colorrange(curve['Curve signal']))

    # experimental data
    x_data = time_points
    y_data = curve[tmt_cols_ratio]
    ax.plot(x_data, y_data, '-o', color=color_sn)

    # add the title
    add_plot_title(ax, curve['Gene names'], curve['Modified sequence'])

    # axes make the plot pretty
    ax.set_xscale('log')
    ax.set_xlabel(f'Time [{time_unit}]')
    ax.set_ylabel('Relative Response')
    ax.set_ylim((0, None if max(y_data) > 1.9 else 2))
    ax.axhline(y=plot_helper_line[0], ls=':', color='black')
    ax.axhline(y=plot_helper_line[1], ls=':', color='black')


def write_curves_to_pdf(pdf_file_path, df, intensity_distribution, cmap, colorrange, tmt_cols_ratio, experiment_type,
                        drug_log_concs, time_points, curve_x_extrapolate, plot_helper_line, drug_unit, time_unit,
                        plot_logisitc_flag):
    """
    Write the curves to pdf.

    This is a wrapper function that calls the actual plotting function for each curve present in the DataFrame.
    Depending on the data type (dd or td) the respective plotting function is used.
    """
    with PdfPages(pdf_file_path) as pdf:
        # For each curve, create a figure, plot the data, and add to the opened pdf
        for i, row in tqdm.tqdm(df.iterrows(), total=len(df)):

            fig, (plot_ax, cbar_ax) = plt.subplots(nrows=1, ncols=2, figsize=(8, 4),
                                                   gridspec_kw={'width_ratios': [14, 1]}, dpi=100)

            add_intensity_dist(cbar_ax, intensity_distribution, cmap, colorrange, row['Curve signal'])

            if experiment_type == 'dd':
                add_dd_plot(plot_ax, row, cmap, colorrange, tmt_cols_ratio, drug_log_concs, drug_unit,
                            curve_x_extrapolate, plot_helper_line,  plot_logisitc_flag)

            elif experiment_type == 'td':
                add_td_plot(plot_ax, row, cmap, colorrange, tmt_cols_ratio, time_points, time_unit, plot_helper_line)

            else:
                raise ValueError('Experiment type is unknown. Set it to td (time) or dd (dose)')

            fig.tight_layout()
            pdf.savefig(fig)
            plt.close()


def run_quantification(config):
    """
    Run the quantification pipeline on evidence level using MaxQuant output.

    Parameters
    ----------
    config : toml object
        The config file with the parameters imported (dict)

    Return
    -----
    None
    """
    print('Hello from the quantification engine', end='\n\n')

    # Extract the relevant info from the config file
    in_file_path = os.path.join(config['MaxQuant'].get('txt_foler'), config['MaxQuant'].get('input_file'))
    out_file_path = os.path.join(config['Paths'].get('table_out_folder'), config['Paths'].get('output_file'))
    quantification_type = config['TMT'].get('quantification')
    drug_concs = config['TMT'].get('doses')
    drug_log_concs = build_drug_log_concentrations(drug_concs, float(config['TMT'].get('dose_scale')))
    reporter_ids = config['TMT'].get('channels')
    control_channel = config['Processing'].get('control_channel')
    subproteome_normalization = config['Processing'].get('subproteome_normalization')
    ptms = get_expected_ptms(config)
    min_dmso_signal = float(config['CurveFilter'].get('signal_min'))

    # build the new column names
    tmt_cols_raw = build_reporter_col_names('Reporter intensity corrected', reporter_ids)
    tmt_cols_normal = build_reporter_col_names('TMT Channel Normal', reporter_ids)
    tmt_cols_ratio = build_reporter_col_names('TMT Channel Ratio', reporter_ids)

    # specify the information cols
    intensity_cols = ['Intensity'] + tmt_cols_raw
    info_cols = config['MaxQuant'].get('evd_info_cols').copy() + list(v[1] for v in ptms.values())

    # Load evidence file
    print('Reading evidence file...', end=' ')
    evd = load_cleaned_evidence(in_file_path,
                                sum_cols=intensity_cols,
                                first_cols=info_cols,
                                max_cols=config['MaxQuant'].get('evd_max_cols'),
                                min_cols=config['MaxQuant'].get('evd_min_cols'),
                                concat_cols=config['MaxQuant'].get('evd_concat_cols'))

    # Add ptm / fullproteome mask columns
    for ptm_proteome, ptm_col, ptm_symbol in ptms.values():
        evd = add_ptm_mask_column(evd, ptm_proteome, ptm_symbol)
    evd = add_fullproteome_mask_column(evd, ptms)
    print('Done!')

    # Median Normalization grouped on experiment
    evd = evd.groupby('Experiment').apply(add_median_normalization, raw_cols=tmt_cols_raw, norm_cols=tmt_cols_normal,
                                          ref_col=subproteome_normalization)
    normalization_factors = evd.groupby('Experiment').apply(calculate_median_normalization_factors, cols=tmt_cols_raw,
                                                            ref_col=subproteome_normalization)

    # ratios to relative to DMSO control
    evd = add_ratios(evd, tmt_cols_normal, tmt_cols_ratio, ref_col=f'TMT Channel Normal {int(control_channel)}')

    # Define 'Curve signal' as the Log2 intensity of the control channel
    evd['Curve signal'] = np.log2(evd[f'TMT Channel Normal {int(control_channel)}'].replace(0, 1))

    # Fit the logistic model
    if config['Processing'].get('fit_logistic_model'):
        print('Fitting Models:')
        flag = (config['Processing'].get('simple_fit') is False)
        i_kwargs = {'exclude_low_n': int(config['Processing'].get('interpolation_low_exclude')),
                    'exclude_top_n': int(config['Processing'].get('interpolation_high_exclude'))}
        evd = add_logistic_model(evd, tmt_cols_ratio, x_data=drug_log_concs, interpolation=flag, max_opt=150,
                                 interpolation_kwargs=i_kwargs)

        # Add a regulation column
        evd = define_regulated_curves(evd, drug_concs, tmt_cols_ratio, min_dmso_signal, quantification_type)

    # Save the data
    print('Writing quantification results to file...', end=' ')
    evd.to_csv(out_file_path, sep='\t', index=False)
    normalization_factors.to_csv(os.path.join(os.path.dirname(out_file_path), 'normalization_factors.txt'), sep='\t')
    print(f'Done!', end='\n\n')


def run_plotting(config):
    """
    Run the plotting pipeline based on the curves file.

    Parameters
    ----------
    config : toml object
        The config file with the parameters imported (dict)

    Return
    -----
    None
    """
    print('Hello from the plotting engine')

    # Extract the relevant info from the config file
    experiment_type = config['Meta Data'].get('experiment_type')
    drug_concs = np.array(config['TMT'].get('doses'))
    time_points = np.array(config['TMT'].get('time_points'))
    reporter_ids = np.array(config['TMT'].get('channels'))
    curve_x_extrapolate = config['Plot'].get('curve_x_extrapolate')
    plot_helper_line = config['Plot'].get('plot_helper_line')
    plot_logisitc_flag = config['Processing'].get('fit_logistic_model')
    drug_unit = str(config['TMT'].get('dose_label'))
    time_unit = str(config['TMT'].get('time_label'))

    # Process toml input
    mask = len(reporter_ids) * [True]
    drug_log_concs = np.array([])

    if experiment_type == 'dd':
        mask = (drug_concs != 0)
        drug_log_concs = build_drug_log_concentrations(drug_concs[mask], float(config['TMT'].get('dose_scale')))
    elif experiment_type == 'td':
        mask = (time_points != 0)
        time_points = time_points[mask]

    tmt_cols_ratio = build_reporter_col_names('TMT Channel Ratio', reporter_ids[mask])
    sort_channel = get_sorting_channel(tmt_cols_ratio, experiment_type, drug_log_concs, time_points)

    # Build the expected ptms
    ptms = get_expected_ptms(config)
    if config['PTMs'].get('Fullproteome'):
        ptms['Fullproteome'] = ('Fullproteome', '', '')

    # Load the curve file
    curve_file_path = os.path.join(config['Paths'].get('table_out_folder'), config['Paths'].get('output_file'))
    curve_df = pd.read_csv(curve_file_path, sep='\t')

    # Apply the min signal filter
    curve_df = curve_df[curve_df['Curve signal'] > float(config['CurveFilter'].get('signal_min'))].reset_index(drop=True)

    # Plot the curves to pdf. Each experiment and subproteome in an individual file.
    print('\nStart with experiment export 1/{}'.format(curve_df['Experiment'].nunique()))
    for exp_name, exp_df in curve_df.groupby('Experiment'):
        for ptm_proteome, ptm_col, ptm_symbol in ptms.values():

            # Signal coloring and Curve signal distributions for plots:
            cmap = matplotlib.cm.get_cmap(config['Plot'].get('colormap_color'))
            colorrange = matplotlib.colors.Normalize(vmin=config['Plot'].get('colormap_lower'),
                                                     vmax=config['Plot'].get('colormap_upper'))

            intensity_distribution = exp_df.loc[exp_df[ptm_proteome], 'Curve signal'].copy().replace(0, np.nan).dropna().values

            # Plot all curves
            pdf_file = f'{exp_name}_{ptm_proteome}_all.pdf'
            pdf_file_path = os.path.join(config['Paths'].get('pdf_out_folder'), pdf_file)

            # Sort the data
            exp_df_sorted = sort_data(exp_df[exp_df[ptm_proteome]], sort_channel)
            if experiment_type == 'td':
                exp_df_sorted = sort_td_data(exp_df[exp_df[ptm_proteome]], tmt_cols_ratio)

            # Write plots to file
            print(pdf_file)
            write_curves_to_pdf(pdf_file_path, exp_df_sorted, intensity_distribution, cmap, colorrange, tmt_cols_ratio,
                                experiment_type, drug_log_concs, time_points, curve_x_extrapolate, plot_helper_line,
                                drug_unit, time_unit, plot_logisitc_flag)

            # Plot the up and down regulated curves if specified in the toml file
            if not config['Plot'].get('up_down_pdfs'):
                continue

            # Plot the up / down going curves into individual PDFs
            for r_type, r_df in exp_df[exp_df['Regulation'].notna()].groupby('Regulation', as_index=False):
                # Paths
                pdf_file = f'{exp_name}_{ptm_proteome}_{r_type}.pdf'
                pdf_file_path = os.path.join(config['Paths'].get('pdf_out_folder'), pdf_file)

                # Sort the data
                r_df_sorted = sort_data(r_df[r_df[ptm_proteome]], sort_channel)
                if experiment_type == 'td':
                    r_df_sorted = sort_td_data(r_df[r_df[ptm_proteome]], tmt_cols_ratio)

                    # Write plots to file
                print(pdf_file)
                write_curves_to_pdf(pdf_file_path, r_df_sorted, intensity_distribution, cmap, colorrange, tmt_cols_ratio,
                                    experiment_type, drug_log_concs, time_points, curve_x_extrapolate, plot_helper_line,
                                    drug_unit, plot_logisitc_flag)


# ###### #
# Script #
# ###### #

if __name__ == '__main__':
    # Add a command line parser for the config file
    parser = argparse.ArgumentParser(
        description='Run ddPTM data analysis based on the evidence',
        epilog="It will analyze the evidence file corresponding to IDs and paths specified in the toml file.\n\nFPB-2022"
    )
    parser.add_argument(
        dest="config_path",
        metavar="config.toml",
        type=str,
        help="Path to the config.toml file to run the script",
    )
    parser.add_argument(
        '--quantification', '-q',
        action='store_true',
        help='Run the quantification',
    )
    parser.add_argument(
        '--plotting', '-p',
        action='store_true',
        help='Run the plotting',
    )

    print('\n\nEVIDENCE BASED PIPELINE')
    print(80*'~')
    print('Script runs on: ', sys.platform)

    # Parse the config file
    args = parser.parse_args()
    check_path(args.config_path)
    config_file = toml.load(args.config_path)
    check_toml(config_file)
    abs_folder_path = os.path.abspath(os.path.dirname(args.config_path))

    print(80*'~', end='\n\n')

    # Overwrite the relative paths with absolute paths in the config file
    config_file['MaxQuant']['txt_foler'] = os.path.join(abs_folder_path, config_file['MaxQuant']['txt_foler'])
    config_file['Paths']['table_out_folder'] = os.path.join(abs_folder_path, config_file['Paths']['table_out_folder'])
    config_file['Paths']['pdf_out_folder'] = os.path.join(abs_folder_path, config_file['Paths']['pdf_out_folder'])

    # linux paths use forward slashes instead of backslashes
    if sys.platform.startswith('linux'):
      config_file['MaxQuant']['txt_foler'] = config_file['MaxQuant']['txt_foler'].replace("\\", "/")
      config_file['Paths']['table_out_folder'] = config_file['Paths']['table_out_folder'].replace("\\", "/")
      config_file['Paths']['pdf_out_folder'] = config_file['Paths']['pdf_out_folder'].replace("\\", "/")
    
    # Extract the folder paths from the config file and check them
    check_path(config_file['MaxQuant'].get('txt_foler'))
    check_path(config_file['Paths'].get('table_out_folder'))
    check_path(config_file['Paths'].get('pdf_out_folder'), create_folder=True)

    # Run Quantification
    if args.quantification:
        run_quantification(config_file)

    # Run Plotting
    if args.plotting:
        run_plotting(config_file)

    print('\nSCRIPT DONE !!', end='\n\n\n')
