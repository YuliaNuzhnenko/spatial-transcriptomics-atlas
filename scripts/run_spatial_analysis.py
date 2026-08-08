#!/usr/bin/env python3
"""
Spatial Transcriptomics Autocorrelation & Neighborhood Analyzer
Real Data Processor for 10x Genomics Visium Spatial Coordinates
Author: Yulia Nuzhnenko
"""
import os
import numpy as np
import pandas as pd

def load_visium_sample_data():
    """
    Loads real 10x Visium tissue spot coordinates and expression from CSV.
    """
    csv_path = os.path.join(os.path.dirname(__file__), "..", "examples", "data", "visium_brain_spatial.csv")
    df = pd.read_csv(csv_path)
    
    coords = df[["x", "y"]].values
    expression_matrix = df[["Mbp_expression", "Actb_expression"]].values
    gene_names = ["Mbp (Myelin Basic Protein)", "Actb (Beta-Actin)"]
    
    return coords, expression_matrix, gene_names

def compute_morans_i(coords, gene_expression, distance_cutoff=200.0):
    """
    Computes Moran's I spatial autocorrelation statistic.
    Formula: I = (N / S0) * sum_i sum_j w_ij (x_i - mean)(x_j - mean) / sum_i (x_i - mean)^2
    """
    N = len(coords)
    mean_x = np.mean(gene_expression)
    denom = np.sum((gene_expression - mean_x) ** 2)
    if denom == 0:
        return 0.0
        
    diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
    dist_matrix = np.sqrt(np.sum(diff ** 2, axis=-1))
    
    W = (dist_matrix < distance_cutoff).astype(float)
    np.fill_diagonal(W, 0)
    
    S0 = np.sum(W)
    if S0 == 0:
        return 0.0
        
    delta = gene_expression - mean_x
    num = np.sum(W * np.outer(delta, delta))
    
    return float((N / S0) * (num / denom))

def main():
    print("==================================================")
    print(" 10x Visium Spatial Autocorrelation Processor")
    print("==================================================")
    coords, expr, genes = load_visium_sample_data()
    print(f"Loaded 10x Visium Dataset: {len(coords)} spatial tissue spots.")
    print("Calculating Moran's I spatial autocorrelation per marker gene:\n")
    
    for i, gene in enumerate(genes):
        score = compute_morans_i(coords, expr[:, i])
        print(f"  * Gene: {gene:<30} Moran's I: {score:.4f}")

if __name__ == "__main__":
    main()
