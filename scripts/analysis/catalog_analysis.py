#!/usr/bin/env python3
"""
Catalog analysis tools for galaxy properties.

Based on 01_data_exploration.ipynb notebook analysis.
Tools for analyzing galaxy catalog data including:
- Selection criteria feasibility analysis
- NIRSpec FoV constraint validation
- Data quality assessment
- Key parameter distributions
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
import warnings
warnings.filterwarnings('ignore')

# Define cosmology (Planck 2018)
cosmo = FlatLambdaCDM(H0=67.4, Om0=0.315)


def load_selection_data(filepath):
    """
    Load the selection criteria data from CSV file.
    
    Parameters:
    -----------
    filepath : str
        Path to the selection data CSV file
        
    Returns:
    --------
    pandas.DataFrame
        Loaded selection data
    """
    try:
        selection_data = pd.read_csv(filepath)
        print(f"Selection data loaded successfully!")
        print(f"Shape: {selection_data.shape}")
        print(f"Total galaxies: {len(selection_data):,}")
        return selection_data
    except Exception as e:
        print(f"Error loading data: {e}")
        return None


def analyze_data_quality(catalog):
    """
    Check data quality and missing values.
    
    Parameters:
    -----------
    catalog : pandas.DataFrame
        Galaxy catalog data
        
    Returns:
    --------
    dict
        Dictionary with data quality metrics
    """
    print("=== DATA QUALITY ASSESSMENT ===")
    
    # Missing values analysis
    missing_values = catalog.isnull().sum()
    missing_percent = (missing_values / len(catalog)) * 100
    missing_df = pd.DataFrame({
        'Missing Count': missing_values,
        'Missing Percentage': missing_percent
    }).sort_values('Missing Count', ascending=False)
    
    print("Missing values per column:")
    missing_with_data = missing_df[missing_df['Missing Count'] > 0]
    if len(missing_with_data) > 0:
        print(missing_with_data)
    else:
        print("No missing values found!")
    
    # Key variables summary
    print("\n=== Key Variables Summary ===")
    key_vars = ['z', 'logm', 'sfr', 'Re_kpc', 'logSFR', 'DeltaMS']
    summary = {}
    
    for var in key_vars:
        if var in catalog.columns:
            valid_data = catalog[var].dropna()
            summary[var] = {
                'min': valid_data.min(),
                'max': valid_data.max(),
                'median': valid_data.median(),
                'count': len(valid_data)
            }
            print(f"{var}: min={valid_data.min():.3f}, max={valid_data.max():.3f}, median={valid_data.median():.3f}, count={len(valid_data)}")
    
    return {
        'missing_analysis': missing_df,
        'key_variables': summary,
        'total_objects': len(catalog)
    }


def analyze_selection_feasibility(catalog):
    """
    Analyze feasibility of selection criteria for JWST observations.
    
    Parameters:
    -----------
    catalog : pandas.DataFrame
        Galaxy catalog with selection criteria columns
        
    Returns:
    --------
    dict
        Dictionary with feasibility analysis results
    """
    print("=== SELECTION CRITERIA FEASIBILITY ===")
    
    results = {}
    
    # 1. Star-forming galaxies criterion (DeltaMS >= 0)
    sf_galaxies = catalog[catalog['DeltaMS'] >= 0]
    sf_fraction = len(sf_galaxies) / len(catalog)
    results['star_forming'] = {
        'count': len(sf_galaxies),
        'total': len(catalog),
        'fraction': sf_fraction
    }
    print(f"1. Star-forming galaxies (ΔMS ≥ 0): {len(sf_galaxies):,} / {len(catalog):,} ({sf_fraction*100:.1f}%)")
    
    # 2. Redshift range for Halpha observability with NIRSpec
    # NIRSpec wavelength range: 0.97-5.27 μm, Halpha rest: 0.6563 μm
    z_min_halpha = (0.97 / 0.6563) - 1  # ~0.48
    z_max_halpha = (5.27 / 0.6563) - 1  # ~7.0
    
    halpha_observable = catalog[(catalog['z'] >= z_min_halpha) & (catalog['z'] <= z_max_halpha)]
    halpha_fraction = len(halpha_observable) / len(catalog)
    results['halpha_observable'] = {
        'count': len(halpha_observable),
        'z_range': (z_min_halpha, z_max_halpha),
        'fraction': halpha_fraction
    }
    print(f"2. Halpha observable ({z_min_halpha:.2f} < z < {z_max_halpha:.1f}): {len(halpha_observable):,} / {len(catalog):,} ({halpha_fraction*100:.1f}%)")
    
    # 3. Combined star-forming + Halpha observable
    sf_and_halpha = sf_galaxies[(sf_galaxies['z'] >= z_min_halpha) & (sf_galaxies['z'] <= z_max_halpha)]
    combined_fraction = len(sf_and_halpha) / len(catalog)
    results['combined_criteria'] = {
        'count': len(sf_and_halpha),
        'fraction': combined_fraction
    }
    print(f"3. Star-forming AND Halpha observable: {len(sf_and_halpha):,} / {len(catalog):,} ({combined_fraction*100:.1f}%)")
    
    # 4. Size constraint analysis
    valid_size = sf_and_halpha.dropna(subset=['Re_kpc'])
    size_fraction = len(valid_size) / len(sf_and_halpha) if len(sf_and_halpha) > 0 else 0
    results['valid_size'] = {
        'count': len(valid_size),
        'fraction': size_fraction
    }
    print(f"4. With valid size measurements: {len(valid_size):,} / {len(sf_and_halpha):,} ({size_fraction*100:.1f}% of SF+Halpha sample)")
    
    # 5. Sample properties summary
    if len(valid_size) > 0:
        print(f"\n=== AVAILABLE SAMPLE PROPERTIES ===")
        results['sample_properties'] = {
            'mass_range': (valid_size['logm'].min(), valid_size['logm'].max()),
            'redshift_range': (valid_size['z'].min(), valid_size['z'].max()),
            'size_range': (valid_size['Re_kpc'].min(), valid_size['Re_kpc'].max()),
            'deltams_range': (valid_size['DeltaMS'].min(), valid_size['DeltaMS'].max())
        }
        print(f"Stellar mass range: {valid_size['logm'].min():.2f} - {valid_size['logm'].max():.2f} dex")
        print(f"Redshift range: {valid_size['z'].min():.3f} - {valid_size['z'].max():.3f}")
        print(f"Effective radius range: {valid_size['Re_kpc'].min():.2f} - {valid_size['Re_kpc'].max():.2f} kpc")
        print(f"ΔMS range: {valid_size['DeltaMS'].min():.2f} - {valid_size['DeltaMS'].max():.2f}")
    
    # 6. Pair selection feasibility
    target_pairs = 10
    target_control = 10
    total_needed = (target_pairs * 2) + target_control  # 30 total
    
    print(f"\n=== PAIR SELECTION FEASIBILITY ===")
    print(f"Available galaxies for pair selection: {len(valid_size):,}")
    print(f"Target: {target_pairs} pairs ({target_pairs * 2} galaxies) + {target_control} control = {total_needed} total galaxies")
    
    selection_ratio = total_needed / len(valid_size) if len(valid_size) > 0 else float('inf')
    results['selection_feasibility'] = {
        'available': len(valid_size),
        'needed': total_needed,
        'ratio': selection_ratio,
        'feasible': len(valid_size) >= total_needed
    }
    
    if len(valid_size) >= total_needed:
        print(f"✓ Selection is feasible! Ratio needed: {selection_ratio*100:.3f}%")
    else:
        print(f"⚠ May need to relax constraints. Available: {len(valid_size)}, Needed: {total_needed}")
    
    return results, valid_size


def check_nirspec_fov_constraints(catalog, nirspec_fov_arcsec=3.0):
    """
    Check NIRSpec IFU field of view constraints for galaxy sample.
    
    Parameters:
    -----------
    catalog : pandas.DataFrame
        Galaxy catalog with z and Re_kpc columns
    nirspec_fov_arcsec : float
        NIRSpec IFU field of view in arcseconds (default: 3.0)
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame with FoV analysis results
    """
    print("=== NIRSpec IFU FIELD OF VIEW CONSTRAINT ===")
    
    def check_fov_constraint(z, re_kpc, n_re=2):
        """Check if n_re effective radii fit within NIRSpec FoV"""
        # Convert effective radius from kpc to arcsec at redshift z
        angular_diameter_distance = cosmo.angular_diameter_distance(z).to(u.kpc).value
        re_arcsec = (re_kpc / angular_diameter_distance) * 206265  # conversion to arcsec
        
        # Check if 2*Re fits within FoV
        diameter_needed = n_re * re_arcsec
        return diameter_needed <= nirspec_fov_arcsec, diameter_needed, re_arcsec
    
    # Apply FoV constraint to catalog
    fov_results = []
    for idx, galaxy in catalog.iterrows():
        if pd.notna(galaxy['z']) and pd.notna(galaxy['Re_kpc']):
            fits_fov, diameter_needed, re_arcsec = check_fov_constraint(galaxy['z'], galaxy['Re_kpc'])
            fov_results.append({
                'id': galaxy.get('id', idx),
                'z': galaxy['z'], 
                'Re_kpc': galaxy['Re_kpc'],
                'Re_arcsec': re_arcsec,
                'diameter_needed': diameter_needed,
                'fits_fov': fits_fov
            })
    
    fov_df = pd.DataFrame(fov_results)
    
    if len(fov_df) > 0:
        # Summary statistics
        n_fits_fov = fov_df['fits_fov'].sum()
        success_rate = n_fits_fov / len(fov_df)
        
        print(f"Galaxies fitting NIRSpec FoV (2Re ≤ {nirspec_fov_arcsec}''): {n_fits_fov} / {len(fov_df)} ({success_rate*100:.1f}%)")
        
        return fov_df, success_rate
    else:
        print("No valid galaxies found for FoV analysis")
        return pd.DataFrame(), 0.0


def plot_key_distributions(catalog, save_path=None):
    """
    Create comprehensive distribution plots of key galaxy properties.
    
    Parameters:
    -----------
    catalog : pandas.DataFrame
        Galaxy catalog data
    save_path : str, optional
        Path to save the plot
    """
    # Set up plotting style
    plt.style.use('default')
    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['legend.fontsize'] = 10
    
    # Create comprehensive distribution plots
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()
    
    # 1. Redshift distribution
    if 'z' in catalog.columns:
        axes[0].hist(catalog['z'], bins=50, alpha=0.7, edgecolor='black')
        axes[0].set_xlabel('Redshift (z)')
        axes[0].set_ylabel('Number of galaxies')
        axes[0].set_title('Redshift Distribution')
        median_z = catalog['z'].median()
        axes[0].axvline(median_z, color='red', linestyle='--', label=f'Median: {median_z:.2f}')
        axes[0].legend()
    
    # 2. Stellar mass distribution
    if 'logm' in catalog.columns:
        axes[1].hist(catalog['logm'], bins=50, alpha=0.7, edgecolor='black')
        axes[1].set_xlabel('log(M*/M☉)')
        axes[1].set_ylabel('Number of galaxies')
        axes[1].set_title('Stellar Mass Distribution')
        median_mass = catalog['logm'].median()
        axes[1].axvline(median_mass, color='red', linestyle='--', label=f'Median: {median_mass:.2f}')
        axes[1].legend()
    
    # 3. Star formation rate distribution
    if 'sfr' in catalog.columns:
        sfr_log = np.log10(catalog['sfr'].replace(0, np.nan).dropna())
        axes[2].hist(sfr_log, bins=50, alpha=0.7, edgecolor='black')
        axes[2].set_xlabel('log(SFR) [M☉/yr]')
        axes[2].set_ylabel('Number of galaxies')
        axes[2].set_title('Star Formation Rate Distribution')
    
    # 4. Effective radius distribution
    if 'Re_kpc' in catalog.columns:
        valid_re = catalog['Re_kpc'].dropna()
        if len(valid_re) > 0:
            axes[3].hist(valid_re, bins=50, alpha=0.7, edgecolor='black')
            axes[3].set_xlabel('Effective Radius (kpc)')
            axes[3].set_ylabel('Number of galaxies')
            axes[3].set_title('Effective Radius Distribution')
        else:
            axes[3].text(0.5, 0.5, 'No valid Re_kpc data', ha='center', va='center', transform=axes[3].transAxes)
            axes[3].set_title('Effective Radius Distribution')
    
    # 5. Distance from Main Sequence
    if 'DeltaMS' in catalog.columns:
        axes[4].hist(catalog['DeltaMS'], bins=50, alpha=0.7, edgecolor='black')
        axes[4].set_xlabel('Δ(log SFR) from Main Sequence')
        axes[4].set_ylabel('Number of galaxies')
        axes[4].set_title('Distance from SF Main Sequence')
        axes[4].axvline(0, color='red', linestyle='--', label='Main Sequence')
        axes[4].legend()
    
    # 6. Sky distribution
    if 'ra' in catalog.columns and 'dec' in catalog.columns:
        axes[5].scatter(catalog['ra'], catalog['dec'], s=1, alpha=0.5)
        axes[5].set_xlabel('Right Ascension (deg)')
        axes[5].set_ylabel('Declination (deg)')
        axes[5].set_title('Sky Distribution')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Plot saved to: {save_path}")
    
    plt.show()
    
    return fig


def plot_fov_analysis(fov_df, nirspec_fov_arcsec=3.0, save_path=None):
    """
    Plot NIRSpec FoV constraint analysis results.
    
    Parameters:
    -----------
    fov_df : pandas.DataFrame
        FoV analysis results from check_nirspec_fov_constraints
    nirspec_fov_arcsec : float
        NIRSpec FoV size in arcseconds
    save_path : str, optional
        Path to save the plot
    """
    if len(fov_df) == 0:
        print("No data to plot")
        return None
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Effective radius in arcsec vs redshift
    colors = ['green' if fits else 'red' for fits in fov_df['fits_fov']]
    ax1.scatter(fov_df['z'], fov_df['Re_arcsec'], c=colors, alpha=0.6, s=20)
    ax1.axhline(y=nirspec_fov_arcsec/2, color='black', linestyle='--', label=f'NIRSpec FoV limit ({nirspec_fov_arcsec/2}")')
    ax1.set_xlabel('Redshift')
    ax1.set_ylabel('Effective Radius (arcsec)')
    ax1.set_title('Galaxy Size vs Redshift')
    ax1.legend()
    ax1.set_yscale('log')
    
    # Histogram of 2Re sizes
    ax2.hist(fov_df['diameter_needed'], bins=30, alpha=0.7, edgecolor='black')
    ax2.axvline(x=nirspec_fov_arcsec, color='red', linestyle='--', linewidth=2, label=f'NIRSpec FoV ({nirspec_fov_arcsec}")')
    ax2.set_xlabel('2 × Re (arcsec)')
    ax2.set_ylabel('Number of galaxies')
    ax2.set_title('Galaxy Size Distribution')
    ax2.legend()
    ax2.set_yscale('log')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"FoV analysis plot saved to: {save_path}")
    
    plt.show()
    
    return fig


def analyze_mass_function(catalog, mass_col='logm', bins=20):
    """Analyze stellar mass function of selected samples."""
    if mass_col not in catalog.columns:
        print(f"Column {mass_col} not found in catalog")
        return None
    
    valid_masses = catalog[mass_col].dropna()
    hist, bin_edges = np.histogram(valid_masses, bins=bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    plt.figure(figsize=(10, 6))
    plt.bar(bin_centers, hist, width=np.diff(bin_edges), alpha=0.7, edgecolor='black')
    plt.xlabel('log(M*/M☉)')
    plt.ylabel('Number of galaxies')
    plt.title('Stellar Mass Function')
    plt.grid(True, alpha=0.3)
    plt.show()
    
    return {'bin_centers': bin_centers, 'counts': hist, 'bin_edges': bin_edges}


def plot_main_sequence(catalog, mass_col='logm', sfr_col='logSFR'):
    """Plot star formation main sequence."""
    if mass_col not in catalog.columns or sfr_col not in catalog.columns:
        print(f"Required columns not found: {mass_col}, {sfr_col}")
        return None
    
    # Remove invalid data
    valid_data = catalog[[mass_col, sfr_col]].dropna()
    
    plt.figure(figsize=(10, 8))
    plt.scatter(valid_data[mass_col], valid_data[sfr_col], alpha=0.5, s=20)
    plt.xlabel('log(M*/M☉)')
    plt.ylabel('log(SFR) [M☉/yr]')
    plt.title('Star Formation Main Sequence')
    plt.grid(True, alpha=0.3)
    
    # Add main sequence reference line if DeltaMS is available
    if 'logSFR_MS' in catalog.columns:
        ms_line = catalog['logSFR_MS'].dropna()
        mass_for_ms = catalog.loc[ms_line.index, mass_col]
        plt.plot(mass_for_ms, ms_line, 'r-', linewidth=2, label='Main Sequence', alpha=0.8)
        plt.legend()
    
    plt.show()
    
    return valid_data


def compare_pair_vs_control(pair_sample, control_sample, parameters=['logm', 'z', 'logSFR', 'Re_kpc']):
    """Compare properties between pair and control samples."""
    if len(pair_sample) == 0 or len(control_sample) == 0:
        print("Empty samples provided")
        return None
    
    print("=== PAIR vs CONTROL COMPARISON ===")
    print(f"Pair sample size: {len(pair_sample)}")
    print(f"Control sample size: {len(control_sample)}")
    
    comparison_results = {}
    
    for param in parameters:
        if param in pair_sample.columns and param in control_sample.columns:
            pair_values = pair_sample[param].dropna()
            control_values = control_sample[param].dropna()
            
            if len(pair_values) > 0 and len(control_values) > 0:
                pair_median = pair_values.median()
                control_median = control_values.median()
                
                comparison_results[param] = {
                    'pair_median': pair_median,
                    'control_median': control_median,
                    'pair_std': pair_values.std(),
                    'control_std': control_values.std(),
                    'difference': pair_median - control_median
                }
                
                print(f"{param}: Pair={pair_median:.3f}, Control={control_median:.3f}, Diff={pair_median-control_median:.3f}")
    
    return comparison_results


if __name__ == "__main__":
    # Example usage
    print("Catalog Analysis Tools")
    print("======================\n")
    
    # Load example data
    data_path = "../../data/raw/pres_selection_data.csv"
    catalog = load_selection_data(data_path)
    
    if catalog is not None:
        # Analyze data quality
        quality_results = analyze_data_quality(catalog)
        
        # Analyze selection feasibility
        feasibility_results, valid_sample = analyze_selection_feasibility(catalog)
        
        # Check NIRSpec FoV constraints
        fov_results, success_rate = check_nirspec_fov_constraints(valid_sample)
        
        # Create plots
        plot_key_distributions(catalog, save_path="../../results/plots/key_distributions.png")
        
        if len(fov_results) > 0:
            plot_fov_analysis(fov_results, save_path="../../results/plots/fov_analysis.png")
    
    print("Run with data to perform full analysis")