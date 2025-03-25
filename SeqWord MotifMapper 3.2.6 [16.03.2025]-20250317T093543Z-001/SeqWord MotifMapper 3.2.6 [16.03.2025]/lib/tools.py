from functools import reduce
import os, string, random, math, copy, re
import numpy as np

nuc_table = {
    "A":"A",
    "T":"T",
    "G":"G",
    "C":"C",
    "U":"U",
    "R":"[A,G]",
    "Y":"[C,T]",
    "S":"[G,C]",
    "W":"[A,T]",
    "K":"[G,T]",
    "M":"[A,C]",
    "B":"[C,G,T]",
    "V":"[A,C,G]",
    "D":"[A,G,T]",
    "H":"[A,C,T]",
    "N":"[A,C,G,T]",
    }
nucleotides = [["A","T"],["T","A"],["G","C"],["C","G"],["R","Y"],["Y","R"],["S","S"],
        ["W","W"],["K","M"],["M","K"],["B","V"],["V","B"],["D","H"],["H","D"],["N","N"]]

rev_translation = [["A","T"],["T","A"],["G","C"],["C","G"],["R","Y"],["Y","R"],
                    ["K","M"],["M","K"],["B","V"],["V","B"],["D","H"],["H","D"]]

def msg(msg):
    print(f"\n{msg}\n")

def alert(msg):
    print(f"\n{msg}\n")

def format_number(num,dig,zoom=0):
    return str(int((10**(dig+zoom))*num)/float(10**dig))

def format_numeric_string(num):
    outnum = str(num)
    if len(outnum) <= 3:
        return outnum
    for i in range(len(outnum)-3,0,-3):
        outnum = outnum[:i]+","+outnum[i:]
    return outnum

def format_string(text,L=30,flg_fill=False):
    if len(text) < L:
        if flg_fill:
            text += (" "*(L-len(text)))
        return text
    else:
        return text[:L-3]+"..."
    
def get_random_filename(path='',length=10):
    letters = string.ascii_lowercase
    fname = "?"
    while fname:
        fname = "~"+"".join(random.choice(letters) for i in range(length))
        if not os.path.exists(os.path.join(path,fname)):
            break
    return fname

def dereplicate(ls):
    if len(ls) < 2:
        return ls
    ls = copy.deepcopy(ls)
    ls.sort()
    for i in range(len(ls)-1,0,-1):
        if ls[i] == ls[i-1]:
            del ls[i]
    return ls

def dereplicate_and_count(ls, flg_absolute_falues=False):  # ls = [[value, count],...]
    
    def _is_equal(v1, v2, flg_absolute_falues):
        if flg_absolute_falues:
            try:
                v1 = float(v1)
                v2 = float(v2)
                return abs(v1) == abs(v2)
            except:
                v1 = str(v1).upper()
                v2 = str(v2).upper()
                return v1 == v2
        else:
            return v1 == v2
            
    if len(ls) < 2:
        return ls
    ls = copy.deepcopy(ls)
    if flg_absolute_falues:
        ls.sort(key=lambda ls: abs(ls[0]) if (isinstance(ls[0],int) or isinstance(ls[0],float)) else str(ls[0]).upper())
    else:
        ls.sort(key=lambda ls: ls[0])
    for i in range(len(ls) - 1, 0, -1):
        if _is_equal(ls[i][0], ls[i-1][0], flg_absolute_falues):
            ls[i-1][1] += ls[i][1]
            del ls[i]
    return ls

def get_modbase_description(ft):
    if 'description' not in ft.qualifiers:
        return {}
    data = list(map(lambda item: item.split("="), ft.qualifiers['description'][0].replace(" ","").split(";")))
    return dict(zip(list(map(lambda item: item[0], data)),list(map(lambda item: item[1], data))))

def compile_motif(word=""):
    if not word:
        return ""
    letters = list(word.upper())
    return "".join(list(map(lambda l: nuc_table[l], letters)))

def inlist_motifs(motif,flg_reverse_complement = False):
    ls = []
    p = 0
    wlength = len(motif)
    while p < wlength:
        if motif[p] == "[":
            s = motif.find("]",p+1)
            if not ls:
                ls = motif[p+1:s].split(",")
            else:
                ls = reduce(lambda ls1,ls2: ls1+ls2, list(map(lambda l: list(map(lambda w: w+l, ls)), motif[p+1:s].split(","))))
            p = s+1
        else:
            if not ls:
                ls = [motif[p]]
            else:
                ls = list(map(lambda w: w+motif[p], ls))
            p += 1
    if flg_reverse_complement:
        ls = list(map(lambda w: reverse_complement(w), ls))
    return ls
    
def generate_motif(ls):
    if not ls:
        return ""
    if len(ls)==1:
        return ls[0]
    items = nuc_table.items()
    table = dict(zip(list(map(lambda item: item[1], items)),list(map(lambda item: item[0], items))))
    table = {
        "[A]":"A",
        "[T]":"T",
        "[G]":"G",
        "[C]":"C",
        "[U]":"U",
        "[A,G]":"R",
        "[C,T]":"Y",
        "[C,G]":"S",
        "[A,T]":"W",
        "[G,T]":"K",
        "[A,C]":"M",
        "[C,G,T]":"B",
        "[A,G,T]":"D",
        "[A,C,T]":"H",
        "[A,C,G]":"V",
        "[A,C,G,T]":"N"
        }
    motif = []
    for i in range(min([len(item) for item in ls])):
        motif.append(table[str(dereplicate([item[i] for item in ls])).replace("'","").replace(" ","")])
    return "".join(motif)

def reverse_complement(w, reverse=True):
    transtable = "".maketrans("".join(list(map(lambda ls: ls[0], rev_translation))),"".join(list(map(lambda ls: ls[1], rev_translation))))
    w = list(w.upper().translate(transtable))
    if reverse:
        w.reverse()
    return "".join(w)

def is_dna_sequence(seq,ambiguous=True):
    letters = list(map(lambda ls: ls[0], nucleotides))
    if not ambiguous:
        letters = letters[:4]
    return len(list(filter(lambda val: val==False, list(map(lambda l: l.upper() in letters, list(seq))))))==0
    
def get_site_location_relative_to_start_codon(genes,site,strand,modbase_location,reverse_modbase_location):
    if not genes or genes==["non-coding"]:
        return "NA"
    locations = []
    start,end = list(map(lambda v: int(v), site.split("..")))
    for cds in genes:
        if cds.strand==strand:
            p = modbase_location
        elif reverse_modbase_location:
            p = end-reverse_modbase_location+1
        else:
            continue
        if strand == 1:
            locations.append(start-int(cds.location.start)+p-2)
        else:
            locations.append(int(cds.location.end)-end+p-2)
    if len(locations) > 1:
        locations.sort()
        if all(list(map(lambda v: v > 0, locations))):
            locations = [min(locations)]
        elif any(list(map(lambda v: v > 0, locations))):
            locations = [max(locations)]
    return ",".join(list(map(lambda v: str(v), locations)))
    
def get_locus_by_coordinates(seq, start, end, flg_complement=False, flg_reverse=False):    # start and end are counted starting form 1
    locus = seq[int(start)-1:int(end)]
    if flg_complement:
        return reverse_complement(w=locus.upper(), reverse=flg_reverse)
    
    if flg_reverse:
        return "".join(reversed(list(locus.upper())))
    
    return locus.upper()
    
def ask_filename(file_description='file'):
    filename = ""
    while filename == "":
        filename = input(f"Enter name to save {file_description} or * to abort: ").strip()
        if filename == "*":
            return False
        filename = check_file_name_symbols(fname)
    return filename
    
def check_file_name_symbols(fname, replacement="_"):
    for symbol in ("/","\\","\"","|","*","<",">",":","?"):
        fname = fname.replace(symbol,replacement)
    return fname
    
def random_id(length=1):
    id = ""
    for i in range(length):
        id += str(random.randint(0,9))
    return id

def kendal_correlation(values_1, values_2, parameter="", bootstraps=1000, alpha=0.05, output='ci'):
    """
    Calculate the Kendall tau correlation, its statistical error, and a confidence interval.

    Args:
        values_1 (list or array): First set of values.
        values_2 (list or array): Second set of values.
        parameter (str): Optional parameter name for labeling.
        bootstraps (int): Number of bootstrap samples for confidence interval calculation.
        alpha (float): Significance level for confidence interval.
        output (str): 'ci' to return confidence interval, 'stderr' for standard error.

    Returns:
        dict: {'correlation': float, 'correlation_error': float or confidence interval, 'p-value': float}
    """
    if len(values_1) != len(values_2):
        raise ValueError("To calculate Kendall correlation, the two arrays must be of the same length!")

    try:
        from scipy.stats import kendalltau
    except ImportError:
        alert("Module SCIPY is not installed! Kendall tau correlation cannot be performed")
        return None

    import numpy as np

    # Calculate Kendall tau correlation and p-value
    correlation, p_value = kendalltau(values_1, values_2)

    # Calculate the standard error for the correlation
    n = len(values_1)
    if n < 2:
        raise ValueError("Kendall correlation requires at least two data points.")
    corr_error = np.sqrt(2 * (2 * n + 5) / (9 * n * (n - 1)))

    # Bootstrapping for confidence interval
    bootstrap_correlations = []
    for _ in range(bootstraps):
        # Resample data with replacement
        indices = np.random.choice(range(n), size=n, replace=True)
        resampled_1 = np.array(values_1)[indices]
        resampled_2 = np.array(values_2)[indices]
        # Calculate Kendall tau for the resampled data
        bootstrap_corr, _ = kendalltau(resampled_1, resampled_2)
        bootstrap_correlations.append(bootstrap_corr)

    # Calculate confidence interval
    lower_bound = np.percentile(bootstrap_correlations, 100 * alpha / 2)
    if lower_bound < -1:
        lower_bound = -1
    upper_bound = np.percentile(bootstrap_correlations, 100 * (1 - alpha / 2))
    if upper_bound > 1:
        upper_bound = 1
    confidence_interval = (lower_bound, upper_bound)

    # Return results
    return {
        'correlation': correlation,
        'correlation error': confidence_interval if output == 'ci' else corr_error,
        'p-value': p_value,
        'stderr': corr_error,
        'confidence_interval': confidence_interval,
    }
        
def spearman_correlation(values_1, values_2, parameter="", bootstraps=1000, alpha=0.05, output='ci'):
    """
    Calculate the Spearman correlation, its statistical error, and a confidence interval.

    Args:
        values_1 (list or array): First set of values.
        values_2 (list or array): Second set of values.
        parameter (str): Optional parameter name for labeling.
        bootstraps (int): Number of bootstrap samples for confidence interval calculation.
        alpha (float): Significance level for confidence interval.
        output (str): 'ci' to return confidence interval, 'stderr' for standard error.

    Returns:
        dict: {'correlation': float, 'correlation_error': float or confidence interval, 'p-value': float}
    """
    if len(values_1) != len(values_2):
        raise ValueError("To calculate Spearman correlation, the two arrays must be of the same length!")

    try:
        from scipy.stats import spearmanr
    except ImportError:
        alert("Module SCIPY is not installed! Spearman correlation cannot be performed")
        return None

    import numpy as np

    # Calculate Spearman correlation and p-value
    correlation, p_value = spearmanr(values_1, values_2)

    # Calculate the standard error for the correlation
    n = len(values_1)
    if n < 2:
        raise ValueError("Spearman correlation requires at least two data points.")
    corr_error = np.sqrt(2 * (2 * n + 5) / (9 * n * (n - 1)))

    # Bootstrapping for confidence interval
    bootstrap_correlations = []
    for _ in range(bootstraps):
        # Resample data with replacement
        indices = np.random.choice(range(n), size=n, replace=True)
        resampled_1 = np.array(values_1)[indices]
        resampled_2 = np.array(values_2)[indices]
        # Calculate Spearman correlation for the resampled data
        bootstrap_corr, _ = spearmanr(resampled_1, resampled_2)
        bootstrap_correlations.append(bootstrap_corr)

    # Calculate confidence interval
    lower_bound = np.percentile(bootstrap_correlations, 100 * alpha / 2)
    if lower_bound < -1:
        lower_bound = -1
    upper_bound = np.percentile(bootstrap_correlations, 100 * (1 - alpha / 2))
    if upper_bound > 1:
        upper_bound = 1
    confidence_interval = (lower_bound, upper_bound)

    # Return results
    return {
        'correlation': correlation,
        'correlation error': confidence_interval if output == 'ci' else corr_error,
        'p-value': p_value,
        'stderr': corr_error,
        'confidence_interval': confidence_interval,
    }
        
def u_correlation(values_1, values_2, parameter="", bootstraps=1000, alpha=0.05, output='ci'):
    """
    Calculate the Mann-Whitney U Test and its statistical error.

    Args:
        num_resamples (int): Number of resumples for statistical testing (default: 1000).
        values_1 (list or array): First set of values.
        values_2 (list or array): Second set of values.
        values_1 and values_2 can be of different lengths

    Returns:
        dict: {'correlation': U statistics from 0 to n1 * n2 where 0 indicates no overlaps between two distributions and n1 * n2 - complete overlap, 
        'correlation_error': (lower,upper) or float, 'p-value': float, 'stdev': float standard deviations}
    """
    if len(values_1) < 1 or len(values_2) < 1:
        raise ValueError("Both input arrays must contain at least one element.")

    try:
        from scipy.stats import mannwhitneyu
    except:
        alert("Module SCIPY is not installed! Mann-Whitney U Test cannot be performed")
        return None
        
    # Calculate Mann-Whitney U Test to assess relationships without assuming normality or linearity
    u_stat, p_value = mannwhitneyu(values_1, values_2)
    
    # Resampling (bootstrap)
    resampled_u_stats = []
    for _ in range(bootstraps):
        resampled_1 = np.random.choice(values_1, size=len(values_1), replace=True)
        resampled_2 = np.random.choice(values_2, size=len(values_2), replace=True)
        resampled_u_stat, _ = mannwhitneyu(resampled_1, resampled_2)
        resampled_u_stats.append(resampled_u_stat)

    # Calculate the standard error and confidence intervals
    u_stat_std = np.std(resampled_u_stats)
    u_stat_ci = [np.percentile(resampled_u_stats, 97.5), np.percentile(resampled_u_stats, 2.5)]
    if u_stat_ci[0] < 0:
        u_stat_ci[0] = 0

    # Return results
    return {
        'correlation': 1 - u_stat / len(values_1) / len(values_2),  # Convertion to 0..1 range; 0 - no overlap; 1 - complete overlap
        'correlation error': [1 - v  / len(values_1) / len(values_2) for v in u_stat_ci] if output == 'ci' else u_stat_std  / len(values_1) / len(values_2),
        'p-value': p_value,
        'stdev': u_stat_std  / len(values_1) / len(values_2),
        'stderr': None,
        'confidence_interval': u_stat_ci,
    }

def ks_correlation(values_1, values_2, parameter="", bootstraps=1000, alpha=0.05, output='ci'):
    """
    Calculate the Kolmogorov-Smirnov Test and its statistical error.

    Args:
        num_resamples (int): Number of resumples for statistical testing (default: 1000).
        values_1 (list or array): First set of values.
        values_2 (list or array): Second set of values.
        values_1 and values_2 can be of different lengths

    Returns:
        dict: {'correlation': KS statistics from 0 to 1 where 0 Indicates the two distributions are identical and 1 - they are different, 
        'correlation_error': (lower,upper) or float, 'p-value': float, 'stdev': float standard deviations}
    """
    if len(values_1) < 1 or len(values_2) < 1:
        raise ValueError("Both input arrays must contain at least one element.")

    try:
        from scipy.stats import ks_2samp
    except:
        alert("Module SCIPY is not installed! Kolmogorov-Smirnov Test cannot be performed")
        return None
        
    # Calculate Kolmogorov-Smirnov Test to assess relationships without assuming normality or linearity
    ks_stat, p_value = ks_2samp(values_1, values_2)
    
    # Resampling (bootstrap)
    resampled_ks_stats = []
    for _ in range(bootstraps):
        resampled_1 = np.random.choice(values_1, size=len(values_1), replace=True)
        resampled_2 = np.random.choice(values_2, size=len(values_2), replace=True)
        resampled_ks_stat, _ = ks_2samp(resampled_1, resampled_2)
        resampled_ks_stats.append(resampled_ks_stat)

    # Calculate the standard error and confidence intervals
    ks_stat_std = np.std(resampled_ks_stats)
    ks_stat_ci = [1 - np.percentile(resampled_ks_stats, 97.5), 1 - np.percentile(resampled_ks_stats, 2.5)]
    if ks_stat_ci[0] < 0:
        ks_stat_ci[0] = 0

    # Return results
    return {
        'correlation': 1 - ks_stat,     # Convertion to 0..1 range; 0 - no overlap; 1 - complete overlap
        'correlation error': ks_stat_ci if output == 'ci' else ks_stat_std,
        'p-value': p_value,
        'stdev': ks_stat_std,
        'stderr': None,
        'confidence_interval': ks_stat_ci,
    }

def calculate_mutual_info(values_1, values_2, parameter="", rank_number=4, bootstraps=100, alpha=0.05, output='ci'):
    """
    Calculate the Mutual Information Regression using non-parametric methods (e.g., kernel density estimation or k-nearest neighbors),
    including estimates for standard deviation, standard error, and confidence intervals.

    Args:
        values_1 (list or array): Feature variables.
        values_2 (list or array): Dependent target variables.
        bootstraps (int): Number of bootstrap samples for statistical testing.
        alpha (float): Significance level for confidence interval.
        output (str): 'ci' to return confidence interval, 'stderr' for standard error.

    Returns:
        dict: {'correlation': float, 'p-value': float, 'correlation_error': float or CI, 'stdev': float, 'stderr': float, 'confidence_interval': tuple}
    """
    try:
        from sklearn.feature_selection import mutual_info_regression
    except ImportError:
        alert("Module SKLEARN is not installed! Mutual Information Regression cannot be performed")
        return None

    if len(values_1) != len(values_2):
        raise ValueError("Input arrays must have the same length!")

    import numpy as np

    # Convert inputs to NumPy arrays
    X = np.array(values_1).reshape(-1, 1)  # Feature
    y = np.array(values_2)  # Target

    # Calculate observed mutual information
    observed_mi = mutual_info_regression(X, y, random_state=42)[0]

    # Bootstrapping for statistical estimates
    bootstrapped_mis = []
    for _ in range(bootstraps):
        resampled_indices = np.random.choice(len(y), len(y), replace=True)
        resampled_X = X[resampled_indices]
        resampled_y = y[resampled_indices]
        bootstrapped_mi = mutual_info_regression(resampled_X, resampled_y, random_state=42)[0]
        bootstrapped_mis.append(bootstrapped_mi)

    # Convert to numpy array for processing
    bootstrapped_mis = np.array(bootstrapped_mis)

    # Calculate standard deviation
    stdev = np.std(bootstrapped_mis)

    # Calculate standard error
    stderr = stdev / np.sqrt(len(bootstrapped_mis))

    # Calculate confidence interval
    lower_bound = np.percentile(bootstrapped_mis, 100 * alpha / 2)
    if lower_bound < 0:
        lower_bound = 0
    upper_bound = np.percentile(bootstrapped_mis, 100 * (1 - alpha / 2))
    confidence_interval = (lower_bound, upper_bound)

    # Validate bootstrapped values against observed_mi
    if not (np.min(confidence_interval) <= observed_mi <= np.max(confidence_interval)):
        observed_mi = sum(confidence_interval)/2

    # Permutation testing for p-value
    permuted_mis = []
    for _ in range(bootstraps):
        shuffled_y = np.random.permutation(y)
        permuted_mi = mutual_info_regression(X, shuffled_y, random_state=42)[0]
        permuted_mis.append(permuted_mi)

    p_value = np.mean([mi >= observed_mi for mi in permuted_mis])

    # Return results
    return {
        'correlation': observed_mi,
        'p-value': p_value,
        'correlation error': confidence_interval if output == 'ci' else stderr,
        'stdev': stdev,
        'stderr': stderr,
        'confidence_interval': confidence_interval,
    }

def calculate_rank_association(target_parameters, modified_nucleotides, parameter="", rank_number=4, bootstraps=1000, alpha=0.05, output='ci'):
    """
    Analyze association using Kendall Tau, K-S Test, and standard error, with confidence intervals.

    Args:
        modified_nucleotides (list or array): Modified nucleotides per sliding window.
        target_parameters (list or array): DNA parameters per sliding window.
        rank_number (int): Number of ranks to divide sliding windows into.
        alpha (float): Significance level for confidence interval.

    Returns:
        dict: Results including correlations, p-values, K-S test results, and confidence intervals.
        output argument can be 'ci' - confidence interval, or 'stderr' - standard error
    """
    from scipy.stats import kendalltau, ks_2samp
    import numpy as np

    if len(modified_nucleotides) != len(target_parameters):
        raise ValueError("Lengths of modified_nucleotides and target_parameters must be the same.")

    # Sort target parameters and assign ranks
    sorted_indices = np.argsort(target_parameters)
    sorted_params = np.array(target_parameters)[sorted_indices]
    sorted_modified = np.array(modified_nucleotides)[sorted_indices]
    
    # Define rank boundaries
    rank_boundaries = np.linspace(0, len(sorted_params), rank_number + 1, dtype=int)

    # Sum of modified nucleotides for each rank
    rank_sums = []
    for i in range(rank_number):
        start, end = rank_boundaries[i], rank_boundaries[i + 1]
        rank_sums.append(np.sum(sorted_modified[start:end]))
        
    # Normalize ranks to represent their central trend
    rank_midpoints = np.arange(1, rank_number + 1)
    
    #print("\ttools:536", parameter, rank_midpoints, rank_sums)

    # Kendall Tau correlation
    tau_value, tau_p_value = kendalltau(rank_midpoints, rank_sums)

    # Kolmogorov-Smirnov Test between first and last ranks
    ks_stat, ks_p_value = ks_2samp(
        sorted_modified[rank_boundaries[0]:rank_boundaries[1]],
        sorted_modified[rank_boundaries[-2]:rank_boundaries[-1]]
    )

    # Bootstrapping for Kendall Tau for stderr and confidence interval
    bootstrap_taus = []
    for _ in range(bootstraps):
        bootstrap_indices = np.random.choice(len(sorted_modified), len(sorted_modified), replace=True)
        bootstrap_sums = []
        for i in range(rank_number):
            start, end = rank_boundaries[i], rank_boundaries[i + 1]
            bootstrap_sums.append(np.sum(sorted_modified[bootstrap_indices[start:end]]))
        bootstrap_tau, _ = kendalltau(rank_midpoints, bootstrap_sums)
        bootstrap_taus.append(bootstrap_tau)

    stderr_tau = np.std(bootstrap_taus)
    # Calculate confidence interval
    lower_bound = np.percentile(bootstrap_taus, 100 * alpha / 2)
    if lower_bound < -1:
        lower_bound = -1
    upper_bound = np.percentile(bootstrap_taus, 100 * (1 - alpha / 2))
    if upper_bound > 1:
        upper_bound = 1
    confidence_interval = [lower_bound, upper_bound]

    '''
    return {
        'kendall_tau': {'value': tau_value, 'p-value': tau_p_value, 'confidence_interval': confidence_interval, 'standard error': stderr_tau},
        'ks_test': {'statistic': ks_stat, 'p-value': ks_p_value},
    }
    '''

    return {
        'correlation': tau_value,
        'p-value': max(tau_p_value, ks_p_value),
        'correlation error': confidence_interval if output == 'ci' else stderr_tau,
        'stderr': stderr_tau
    }

def calculate_distribution_statistics(seq_length, modified_nucleotides, regions_of_interest, mode="Z-score"):    # mode = LD | Z-score
    """
    Compare distributions of modified nucleotides against specific regions of interest.

    Args:
        seq_length (int): Length of the full sequence.
        modified_nucleotides (list): List of [start, end] for modified nucleotides.
        regions_of_interest (list): List of [start, end] for regions of interest.

    Returns:
        dict: {'linkage_disequilibrium': float, 'p_value': float}
    """
    # Initialize full sequence and region-specific markers
    full_sequence = np.zeros(seq_length)
    region_sequence = np.zeros(seq_length)

    # Mark modified nucleotides
    for start, end, strand in modified_nucleotides:
        full_sequence[start-1:end] = 1

    # Mark regions of interest
    for start, end in regions_of_interest:
        region_sequence[start-1:end] = 1
    
    # Total counts of modifications
    total_modified = np.sum(full_sequence)
    total_unmodified = seq_length - total_modified

    # Overlap of modified nucleotides with regions of interest
    modified_in_regions = int(np.sum(full_sequence * region_sequence))
    unmodified_in_regions = int(np.sum(region_sequence)) - modified_in_regions

    # Modified outside regions
    modified_outside_regions = total_modified - modified_in_regions
    unmodified_outside_regions = total_unmodified - unmodified_in_regions
    region_length = int(np.sum(region_sequence))

    if mode == "LD":
        normalized_LD, ld_error, p_value = contingency_table_statistics(seq_length, region_length, modified_in_regions, 
            modified_outside_regions, unmodified_in_regions, unmodified_outside_regions)
        # Return results
        return {
            'LD': normalized_LD,
            'LD error': ld_error,
            'p-value': p_value
        }
    elif mode == "Z-score":
        z_score, z_error, p_value, expected_count, observed_count = calculate_zsore(target_seqlength=region_length, total_seqlength=seq_length, 
            target_counts=modified_in_regions, total_counts=total_modified)
        # Return results
        return {
            'Z-score': z_score,
            'Z error': z_error,
            'p-value': p_value,
            'expected': expected_count,
            'observed': observed_count,
        }
    else:
        raise ValueError(f"Unknown statistics {mode}!")
    
def contingency_table_statistics(seq_length, region_length, modified_in_regions, modified_outside_regions, unmodified_in_regions, unmodified_outside_regions):
    try:
        from scipy.stats import chi2_contingency, fisher_exact
    except:
        alert("Distribution of modified nucleotides across regions cannot be evaluated. \nLibraris SCIPY or NUMPY are not installed on this computer!")
        return None
    
    # Construct contingency table for statistical testing
    contingency_table = np.array([
        [modified_in_regions, modified_outside_regions],
        [unmodified_in_regions, unmodified_outside_regions]
    ])
    # Statistical Test (Fisher's Exact Test if small counts, Chi-Square otherwise)
    if np.any(contingency_table < 5):
        _, p_value = fisher_exact(contingency_table)
    else:
        _, p_value, _, _ = chi2_contingency(contingency_table)

    # Linkage Disequilibrium normalization
    total_modified = modified_in_regions + modified_outside_regions
    observed_proportion = modified_in_regions / region_length
    expected_proportion = total_modified / seq_length
    max_deviation = max(observed_proportion, 1 - observed_proportion)
    
    # Normalized LD in range [-1, 1]
    normalized_LD = (observed_proportion - expected_proportion) / max_deviation

    # Calculate the error of LD
    ld_error = np.sqrt((observed_proportion * (1 - observed_proportion)) / region_length)
    
    return normalized_LD, ld_error, p_value

# Calculate Z-score of a the observed deviation from expectation normalized by sqrt(expectation)
def calculate_zsore(target_seqlength, total_seqlength, target_counts, total_counts):
    try:
        import scipy.stats as stats
    except:
        tools.alert("Calculation cannot be performed as the library scipy has not been installed on this computer!")
        return None, None, None, None, None
        
    # Expected proportions
    expected = (target_seqlength / total_seqlength) * total_counts
    # Z-values and p-values
    z_score = (target_counts - expected) / (expected ** 0.5)
    p_value = stats.norm.sf(abs(z_score)) * 2
    
    # Calculate z-score error (standard error of the z-score)
    if expected > 0:
        z_error = (expected ** 0.5) / (expected)
    else:
        z_error = 0  # Avoid division by zero
    
    return z_score, z_error, p_value, expected, target_counts
        
def compare_sites_across_contigs(contigs, modified_sites):
    try:
        import numpy as np
        from scipy.stats import chi2_contingency
    except:
        alert("Distribution of modified nucleotides across contigs cannot be evaluated \nas libraris SCIPY or NUMPY are not installed on this computer!")
        return 1
        
    # Calculate contig lengths and total sequence length
    contig_lengths = [end - start + 1 for start, end in contigs]
    total_length = sum(contig_lengths)

    # Count modified sites per contig
    modified_counts = np.zeros(len(contigs))
    for i, (c_start, c_end) in enumerate(contigs):
        for m_start, m_end in modified_sites:
            # Check if the modified site overlaps with the current contig
            overlap_start = max(c_start, m_start)
            overlap_end = min(c_end, m_end)
            if overlap_start <= overlap_end:
                modified_counts[i] += overlap_end - overlap_start + 1

    # Expected number of modifications per contig based on contig length
    total_modified = sum(modified_counts)
    expected_counts = [total_modified * (length / total_length) for length in contig_lengths]

    # Create observed and expected counts
    observed = modified_counts
    expected = expected_counts

    # Perform Chi-Square Test
    chi2_stat, p_value = chi2_contingency([observed, expected], correction=False)[:2]

    return p_value


def calculate_GC(seq):
    if not len(seq):
        return 0
    return 100.0*float(seq.upper().count("G")+seq.upper().count("C"))/len(seq)
    
def calculate_AT(seq):
    if not len(seq):
        return 0
    return 100.0*float(seq.upper().count("A")+seq.upper().count("T"))/len(seq)
    
def calculate_GC_skew(seq):
    GC = seq.upper().count("G")+seq.upper().count("C")
    if not GC:
        return 0
    return float(seq.upper().count("G")-seq.upper().count("C"))/GC
    
def calculate_AT_skew(seq):
    AT = seq.upper().count("A")+seq.upper().count("T")
    if not AT:
        return 0
    return float(seq.upper().count("A")-seq.upper().count("T"))/AT
    
def getMeanAndStDev(data_list):
    if not data_list:
        raise ValueError("The input list cannot be empty.")
    
    try:
        data_list = [float(v) for v in data_list]
    except ValueError:
        raise ValueError("All values must be numbers!")
    
    if len(data_list) < 2:
        return data_list[0], 0.0
    
    s = sum(data_list)
    s2 = sum(v * v for v in data_list)
    n = len(data_list)
    mean = s / n
    variance = (s2 - mean * mean * n) / (n - 1)
    std_dev = math.sqrt(variance)
    
    return mean, std_dev    
def get_MinMaxAvr(DataList):
    try:
        DataList = [float(v) for v in DataList]
    except:
        raise ValueError("All values must be numbers!")
    return min(DataList), max(DataList), sum(DataList)/len(DataList)
    
def joing_sliding_windows(ls):  # it is assumed that sliding windows overlap
    if len(ls) < 2:
        return ls
    ls.sort(key=lambda d: d["start"])
    for i in range(len(ls)-1,0,-1):
        # Check if the counts of two overlapping windows are significantly above the expectation 
        if ls[i]["confirmed"] and ls[i-1]["confirmed"]:
            # Joing the windows and their counts
            # Extend border
            ls[i-1]["end"] = ls[i]["end"]
            # Add counts
            ls[i-1]["count"] += ls[i]["count"]
            # Set maximum p-value
            ls[i-1]["p-value"] = max(ls[i]["p-value"],ls[i-1]["p-value"])
            # Remove previous sliding window
            del ls[i]
    return ls

#def density_distribution(seqlength, loci, window_length, window_step):
def density_distribution(frames, modified_sites, flg_join_confirmed_windows=False):
    sliding_windows = sorted(frames, key=lambda ls: int(ls[0]))
    N = len(modified_sites)
    sliding_windows, counts = count_loci_in_windows(frames=sliding_windows, modified_sites=modified_sites)
    # Seq length is the end of the last sliding window
    seqlength = sliding_windows[-1][-1]
    # Calculate window length
    window_length = sliding_windows[0][1] - sliding_windows[0][0]
    p_values, fdr_cutoff_p = analyze_modified_sites(seqlength, N, window_length, counts)
    try:
        sliding_windows = [{"start":int(sliding_windows[i][0]),
                            "end":int(sliding_windows[i][1]),
                            "count":int(counts[i]),
                            "p-value":float(p_values[i]),
                            "confirmed":p_values[i] <= fdr_cutoff_p}
                            for i in range(len(sliding_windows))]
        if flg_join_confirmed_windows:
            sliding_windows = joing_sliding_windows(sliding_windows)
        return counts, fdr_cutoff_p
    except:
        return counts, 0
    
def count_loci_in_windows(frames, modified_sites):
    counts = []
    # Count loci in each sliding window
    for window_start, window_end in frames:
        count = 0
        for entry in modified_sites:
            locus_start, locus_end = entry[:2]
            # Check if locus falls within the sliding window
            if (locus_start >= window_start and locus_start < window_end) or \
               (locus_end > window_start and locus_end <= window_end):
                count += 1
        counts.append(count)
    return frames, counts

def select_promoter_regions(genes, promoter_length, seqlength):
    def _get_length(start, strand_1, end, strand_2, plength):
        if end - start < plength and strand_1 == strand_2:
            return end - start
        else:
            return plength
        
    promoter_regions = [
        [
            genes[i].location.start - _get_length(
                    start=genes[i-1].location.end if i else 0,
                    strand_1 = genes[i-1].strand if i else genes[i].strand,
                    end=genes[i].location.start,
                    strand_2 = genes[i].strand,
                    plength=promoter_length) \
                if genes[i].strand in (1, '+', 'DIR') else genes[i].location.end + 1,
            genes[i].location.start if genes[i].strand in (1, '+', 'DIR') else \
                genes[i].location.end + _get_length(
                    start=genes[i].location.end,
                    strand_1 = genes[i].strand,
                    end=genes[i+1].location.end if i < len(genes)-1 else seqlength,
                    strand_2 = genes[i+1].strand if i < len(genes)-1 else genes[i].strand,
                    plength=promoter_length) + 1,
            genes[i].strand
        ] for i in range(len(genes))
    ]
    # Check and correct overlapped regions
    overlap = 0
    if len(promoter_regions) > 1:
        for i in range(len(promoter_regions) - 1):
            region_1 = promoter_regions[i]
            region_2 = promoter_regions[i+1]
            # Check if the end of region_1 is after the start of region_2 
            if region_1[1] > region_2[0]:
                overlap += region_1[1] - region_2[0]
    # Return promoter regions in like [[0, 61, 1], [1320, 1353, 1], [2963, 2960, 1], [start, stop, strand], ...
    return promoter_regions, overlap
    
def count_loci_in_genes(genes, modified_sites, loci_to_filter=[], echo=""):
    counts = []
    # Removing sites withing regions set to filter out
    if loci_to_filter:
        for i in range(len(modified_sites) - 1, -1, -1):
            locus_start, locus_end, locus_strand = modified_sites[i]
            for gene_start, gene_end, gene_strand in loci_to_filter:
                # Check if locus falls within filtered regions
                if ((locus_start >= gene_start and locus_start < gene_end) or \
                   (locus_end > gene_start and locus_end <= gene_end)) and gene_strand==locus_strand:
                       del modified_sites[i]
                       break
                       
    # Count loci in each sliding window
    for gene_start, gene_end, gene_strand in genes:
        count = 0
        for locus_start, locus_end, locus_strand in modified_sites:
            # Check if locus falls within the sliding window
            if ((locus_start >= gene_start and locus_start < gene_end) or \
               (locus_end > gene_start and locus_end <= gene_end)) and gene_strand==locus_strand:
                count += 1
                if echo:
                    print(f"\ttarget: {gene_start}..{gene_end}, locus: {locus_start}..{locus_end}, strands: {gene_strand}, {locus_strand}, " + 
                        f"{locus_start - gene_end if gene_strand == 1 else gene_start - locus_start}")
        counts.append(count)
    return sum(counts)

def analyze_modified_sites(sequence_length, N, sliding_window_length, counts):
    from scipy.stats import lognorm, rankdata    
    # Calculate the expected number of sites per window
    expectation = (N / sequence_length) * sliding_window_length

    # Filter counts to avoid invalid values for lognormal fit
    positive_counts = [count for count in counts if count > 0]
    if len(positive_counts) < 2:
        return 1, 0
        
    # Fit lognormal distribution to positive counts
    shape, loc, scale = lognorm.fit(positive_counts, floc=0)

    # Generate the p-values for counts
    p_values = []
    for count in counts:
        if count <= expectation:
            p_values.append(1.0)  # Assign p-value = 1 for counts <= expectation
        else:
            # Estimate p-value for counts > expectation assuming lognormal distribution
            p = 1 - lognorm.cdf(count, shape, loc, scale)
            p_values.append(p)

    # Apply Benjamini-Hochberg FDR estimation
    p_values_filtered = [p for p in p_values if p < 1]
    m = len(p_values_filtered)
    if m == 0:
        return p_values, 0  # No p-values < 1 for FDR estimation
    sorted_p_values = sorted(p_values_filtered)
    adjusted_p_values = [p * m / rank for rank, p in enumerate(sorted_p_values, start=1)]
    fdr_cutoff_p = next((p for p in adjusted_p_values if p <= 0.05), 0)
    # Return results
    return p_values, fdr_cutoff_p

# Exclude modified sites within given genomic regions
def filter_regions(sites, regions):
    """
    Filters sites to exclude those with start or end positions within specified regions.

    Args:
        sites (list of dict): List of site dictionaries, each containing 'start' and 'end' keys.
        regions (list of list): List of regions, where each region is a list [start, end].

    Returns:
        list of dict: Filtered list of sites excluding those within the regions.
    """
    def is_within_regions(start, end, region_bounds):
        for region_start, region_end in region_bounds:
            if region_start <= int(start) <= region_end or region_start <= int(end) <= region_end:
                return True
        return False

    # Convert region boundaries to a list of tuples for easier comparison
    region_bounds = [(region[0], region[1]) for region in regions]

    # Filter out sites that fall within the specified regions
    filtered_sites = [
        site for site in sites
        if not is_within_regions(site['start'], site['end'], region_bounds)
    ]

    return filtered_sites
    

def filtered_motif_refining(records, motifs):   # motifs in the format: [['CRGKGATC', 1, 6, -2], ...]
    '''
    Record example:
    {'genome': 'Escherichia', 'method': 'kinModCall', 'modtype': 'm4C', 'start': '14311', 'end': '14311', 'score': '102', 'strand': '+', 'para': '.', , 'nucleotide': 'C',
        'data': {'coverage': '135', 'context': 'GTGATAAACCTCAACGCTGGCGGGGATCAGCGATTTCATCC', 'IPDRatio': '3.01', 'frac': '0.842', 
        'fracLow': '0.553', 'fracUp': '0.989', 'identificationQv': '50'}}
    '''
    selected_motifs = []
    # Select records corresponding to the given motifs
    for motif_set in motifs:
        points = motif_set[1:]
        # Dereplication is needed when both points are 0
        points = dereplicate(points)
        motif = motif_set[0]
        length = len(motif)
        # Motifs are converted to list of words (reverse complement words for negative locations)
        # First motif is turned into this format: C[A,G]G[G,T]GATC, then the templet is converted into a list of all possible motif words
        reverse_complement_words = []
        direct_words = inlist_motifs(motif = compile_motif(motif), flg_reverse_complement = False)
        if any([p < 0 for p in points]):
            reverse_complement_words = inlist_motifs(motif = compile_motif(motif), flg_reverse_complement = True)
        for k in points:
            words = direct_words
            if k < 0:
                # Reverse complement locations should be used
                words = reverse_complement_words
            # Assumption that context sequence is 41 bp
            p = 20
            if k:
                p = 21 - k
            # Filtering records
            if words:
                selected_motifs += [entry for entry in records if entry["data"]["context"][p:p+length] in words]
    return selected_motifs
                
def filter_motifs(records, motifs):
    """
    Filters a list of records based on motifs and modified nucleotide positions.

    Args:
        ls (list): List of dictionaries, each containing a 'data' key with a 'context' field.
        motifs (list): List of motifs. Each motif is a list with the motif string and nucleotide positions.
        flg_include (bool): If True, include matching records. If False, exclude them.

    Returns:
        list: Filtered list of records.
    """
    # Selected records can be either included or excluded
    flg_positive_filter = False
    records_to_include = []
    records_to_exclude = []
    for motif in motifs:    # motif is like that: CTAG, 3, -3
        if not motif.strip().strip("-"):
            continue
        
        # Motifs with '-' at the beginning are excluded
        flg_include=True
        if motif.startswith("-"):
            flg_include = False
            motif = motif[1:]
        else:
            # At list one motif was selected to be included
            flg_positive_filter = True
        
        # Extract motif sequence (word) and modified site location
        motif = [s.strip() for s in motif.split(",")]
        word = motif[0]
        if len(motif) == 1:
            locations = [0]
        else:
            locations = [int(v) for v in motif[1:]]

        # Convert to motif template, for example, 'GDGCHC' to 'G[A,G,T]GC[A,C,T]C'
        template = compile_motif(word)
        
        # Create list of direct motif instances
        direct_instances = []
        if any([v >= 0 for v in locations]):
            direct_instances = inlist_motifs(template)
            
        # If there are modified bases on reverse strand (negative locations), reverse complement instances of motifs are created
        reverse_instances = []
        if any([v < 0 for v in locations]):
            reverse_instances = inlist_motifs(template, True)
        
        # Word length
        wlength = 0
        if direct_instances or reverse_instances:
            wlength = len(direct_instances[0]) if direct_instances else len(reverse_instances[0])

        # Cycle across modied locations
        for k in locations:
            # Select which word instances use for filtering, direct or reverse complement 
            searched_motifs = direct_instances
            reverse_complement = False
            if k < 0:
                searched_motifs = reverse_instances
                # Convert negative k (counted from the end of the word) to a positive number counted from the beginning of the word
                k = abs(k)
                reverse_complement = True
                
            # Location of the modified cite in the context sequence assuming that the context sequence is 41 bp
            p = [21 - i for i in range(1,wlength+1,1)]
            if k != 0:
                p = 21 - k
            
            # Select entries containing sought motifs at proper locations in context sequences
            if isinstance(p, int):
                selected_entries = [entry for entry in records if entry["data"]["context"][p:p+wlength] in searched_motifs]
            elif isinstance(p, list):
                # Select multiple centered words, if modified base location is not specified
                selected_entries = []
                for l in p:
                    selected_entries += [entry for entry in records if entry["data"]["context"][l:l+wlength] in searched_motifs]
                
            if flg_include:
                # Add selected records
                records_to_include += selected_entries
            else:
                # Remove selected records
                records_to_exclude += selected_entries
                
    # If no motifs were selected to include, all records are included
    if not flg_positive_filter:
        records_to_include = records
        
    # Filter excluded motifs
    if records_to_exclude:
        excluded_sites = [int(entry['start']) for entry in records_to_exclude]
        records_to_include = [entry for entry in records_to_include if int(entry['start']) not in excluded_sites]
        
    # Return filtered modified site entries
    return records_to_include
    

def match(rec, reg_motif, n, motif_template_length):
    """
    Checks if a record's context matches the motif aligned at position n.

    Args:
        rec (dict): A record containing a 'data' key with a 'context' field (41-character string).
        reg_motif (regex): Compiled regex motif to match against the context.
        n (int): Position in the motif to align with the central nucleotide.

    Returns:
        bool: True if the context matches the motif, False otherwise.
    """
    context = rec['data']['context']  # 41-character string
    central_nucleotide_index = 20    # Central nucleotide index in the context sequence

    # Calculate the start and end positions for matching the motif
    start = central_nucleotide_index - (n - 1)
    end = start + motif_template_length  # Use the regex pattern length for alignment

    # Ensure the alignment falls within the bounds of the context
    if start < 0 or end > len(context):
        return False

    # Extract the fragment of the context sequence
    context_fragment = context[start:end]

    # Match the fragment against the motif regex
    return bool(reg_motif.match(context_fragment))
                    
###############################################################################
if __name__ == "__main__":
    print(inlist_motifs("CTRGAWA"))
