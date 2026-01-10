import os
import pandas as pd
import numpy as np
def predict_cpg(sequences: dict, readout: str, scale_requested):
    #If no specific scale is request then defalt to linear
    if scale_requested is None:
        scale_actual = "linear"
    else:
        scale_actual = scale_requested

    predictions = {}
    
    if readout == "point":
        for seq_id, sequence in sequences.items():
            predictions[seq_id] = [cpg_mean(sequence, scale_actual)]
   
    if readout == "track":
        for seq_id, sequence in sequences.items():
            predictions[seq_id] = calculate_cpg_per_base(sequence, scale_actual, window_size=50)

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
    
    Parameters:
    -----------
    sequence : str
        DNA sequence (ACGT)
    window_size : int
        Size of sliding window for calculating local CpG density (default: 100bp)
    
    Returns:
    --------
    pd.DataFrame
        DataFrame with columns:
        - position: base pair position (0-indexed)
        - base: nucleotide at this position
        - is_cpg_start: whether this position starts a CpG dinucleotide
        - cpg_density: local CpG density in surrounding window
        - cpg_count: number of CpGs in surrounding window
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
        
        # Count CpGs in window
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

