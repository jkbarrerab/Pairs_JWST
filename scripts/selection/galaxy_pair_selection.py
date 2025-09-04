#!/usr/bin/env python3
"""
Galaxy pair selection for JWST NIRSpec observations.

This script implements the selection criteria for galaxy pairs in the COSMOSWeb sample:
- Maximum separation of 50 kpc projected in the sky
- Similar stellar masses
- Velocity difference < 300 km/s
- Star-forming galaxies (on/above main sequence)
- Size constraint: 2 effective radii within NIRSpec FoV
- Preferably observable Halpha line with NIRSpec
"""

import pandas as pd
import numpy as np
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
from astropy.coordinates import SkyCoord


def load_selection_data(filepath):
    """Load the selection criteria data from CSV file."""
    pass


def calculate_projected_separation(ra1, dec1, ra2, dec2, z):
    """Calculate projected physical separation between two galaxies."""
    pass


def check_mass_similarity(mass1, mass2, tolerance=0.3):
    """Check if two galaxies have similar stellar masses (within tolerance in dex)."""
    pass


def check_velocity_difference(v1, v2, max_diff=300):
    """Check if velocity difference is less than max_diff km/s."""
    pass


def check_star_forming_criterion(galaxy_data):
    """Check if galaxy is star-forming (on/above main sequence)."""
    pass


def check_size_constraint(r_eff, z):
    """Check if 2*R_eff fits within NIRSpec IFU FoV."""
    pass


def check_halpha_observable(z):
    """Check if Halpha line is observable with NIRSpec at given redshift."""
    pass


def select_galaxy_pairs(data, n_pairs=10):
    """Select galaxy pairs based on all criteria."""
    pass


def select_control_sample(data, n_control=10):
    """Select control sample of isolated galaxies."""
    pass


if __name__ == "__main__":
    # Load data and run selection
    data = load_selection_data("data/raw/selection_criteria_data.csv")
    pairs = select_galaxy_pairs(data)
    control = select_control_sample(data)