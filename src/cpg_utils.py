import numpy as np
import pandas as pd

def predict_cpg(sequences: dict, readout: str, scale_requested, prediction_ranges):
    #If no specific scale is request then default to linear
    if scale_requested is None:
        scale_actual = "linear"
    else:
        scale_actual = scale_requested

    predictions = {}
    
    # NOTE: Do not modify the original sequences dictionary directly,
    # as it may be needed for multiple prediction tasks.
    # Instead, create a copy to work with.
    seqs = {k: v for k, v in sequences.items()}
    
    # prediction_ranges are in FLANKED coordinates (relative to the sequence after
    # any upstream/downstream flanks were applied in preprocess_data).
    #   point -> crop the INPUT so the mean covers exactly the range.
    #   track -> crop the OUTPUT so edge positions keep their full window context
    # If prediction ranges were sent from the Evaluator use them to crop the sequences
    window_size = 50
    half_window = window_size // 2

    if readout == "point":
        for seq_id, sequence in seqs.items():
            pr = prediction_ranges.get(seq_id) if prediction_ranges else None
            if pr:
                start, end = pr
                sequence = sequence[start:end]
                print(f"Sequence '{seq_id}' cropped to range [{start}, {end}] before point calc.")
            predictions[seq_id] = cpg_mean(sequence, scale_actual)
   
    if readout == "track":
        for seq_id, sequence in seqs.items():
            pr = prediction_ranges.get(seq_id) if prediction_ranges else None
            if pr:
                start, end = pr
                lo = max(0, start - half_window)
                hi = min(len(sequence), end + half_window)   # minimal sufficient right pad
                padded = calculate_cpg_per_base(sequence[lo:hi], scale_actual, window_size=window_size)
                left_trim = start - lo
                track = padded[left_trim:left_trim + (end - start)]
                print(f"Track for '{seq_id}' windowed on [{lo}, {hi}) and trimmed to range [{start}, {end}].")
            else:
                track = calculate_cpg_per_base(sequence, scale_actual, window_size=window_size)
            predictions[seq_id] = track

    return predictions, scale_actual


def cpg_mean(seq: str, scale_actual) -> int:
    """Calculate average CpG sites ('CG') in a DNA sequence."""
    s = seq.upper()
    total_cpg = sum(1 for i in range(len(s)-1) if s[i:i+2] == "CG")

    epsilon = 1e-9  # tiny constant for zero smoothing
    mean_val = (total_cpg+ epsilon)/len(s)

    if scale_actual == "linear":
        cpg_mean = mean_val
    if scale_actual == "log":
        cpg_mean = np.log2(mean_val)

    return cpg_mean

def calculate_cpg_per_base(sequence, scale_actual, window_size):
    """
    Calculate CpG density and counts at base pair level using a sliding window.
    
    Calculate CpG density at base pair level using a sliding window.
    
    Parameters:
    -----------
    sequence : str
        DNA sequence (ACGTN)
    scale_actual : str
        "linear" or "log"; log applies log2 with a small epsilon
    window_size : int
        Size of the sliding window centred on each position
    
    Returns:
    --------
    list of float
        One value per base: CpGs per 100 bp within the window around that
        position. Windows are clipped at the sequence edges, so the rate is
        normalised by the actual window length.
    """
    sequence = sequence.upper().replace(' ', '').replace('\n', '')
    seq_len = len(sequence)
    
    # Initialize arrays
    positions = np.arange(seq_len)
    bases = list(sequence)
    is_cpg_start = np.zeros(seq_len, dtype=bool)
    cpg_density = np.zeros(seq_len, dtype=float)
    cpg_count = np.zeros(seq_len, dtype=int)
    
    # Mark CpG dinucleotide start positions
    for i in range(seq_len - 1):
        if sequence[i:i+2] == 'CG':
            is_cpg_start[i] = True
    
    # Calculate sliding window CpG density for each position
    half_window = window_size // 2
    
    for i in range(seq_len):
        # Define window boundaries
        start = max(0, i - half_window)
        end = min(seq_len, i + half_window)
        
        # Count CpGs WHOLLY INSIDE the window. A CpG whose C is the window's last base
        # and whose G falls just past the edge is not counted
        # NOTE this differs from cpg_mean (point readout), which scans the whole
        # sequence and so has no window boundary to straddle.
        window_seq = sequence[start:end]
        cpgs_in_window = window_seq.count('CG')
        window_len = end - start
        
        cpg_count[i] = cpgs_in_window
        # Density as CpGs per 100 bp
        cpg_density[i] = (cpgs_in_window / window_len) * 100 if window_len > 0 else 0
    
    # Create DataFrame
    df = pd.DataFrame({
        'position': positions,
        'base': bases,
        'is_cpg_start': is_cpg_start,
        'cpg_density': cpg_density,
        'cpg_count': cpg_count
    })

    epsilon = 1e-9  # tiny constant for zero smoothing

    if scale_actual == "linear":
        cpg_bp = df['cpg_density'].tolist()
    if scale_actual == "log":
        cpg_bp = np.log2(df['cpg_density'].to_numpy() + epsilon).tolist()
    
    return cpg_bp
