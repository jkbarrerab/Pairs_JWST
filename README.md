# JWST Galaxy Pairs Analysis

## Project Description
In this project I am going to use the file selection_criteria_data.csv to select 10 galaxies in pairs in the COSMOSWeb sample as well as 10 galaxies as control sample. The criteria is the following:

- The galaxy in the pair should have a maximum separation of 50 kpc projected in the sky with its companion
- The two galaxies in the pair should have similar stellar masses
- The difference in velocity should be less than 300 km/s within the pair
- The galaxies should be star-forming this is within or above the star formation main sequence
- The size of the galaxies in the samples should be that at least two effective radii of the galaxies is within the IFU NIRSpec FoV
- It is preferable that the galaxies have redshift such that it is possible to observe at least the Halpha line using the NIRSpec Instrument

## Repository Structure

```
├── .gitignore
├── README.md
├── config/
│   └── selection_config.yaml     # Configuration parameters for selection
├── data/
│   ├── raw/                      # COSMOSWeb catalog and selection data
│   └── processed/                # Cleaned/filtered data
├── scripts/
│   ├── selection/                # Galaxy pair selection algorithms
│   │   └── galaxy_pair_selection.py
│   ├── analysis/                 # Catalog analysis tools
│   │   └── catalog_analysis.py
│   └── visualization/            # Plotting utilities
│       └── plotting_utils.py
├── notebooks/                    # Jupyter notebooks for exploration
│   ├── 01_data_exploration.ipynb
│   ├── 02_galaxy_pair_selection.ipynb
│   └── 03_sample_analysis.ipynb
├── results/
│   ├── catalogs/                 # Selected samples (CSV)
│   ├── plots/                    # Figures and visualizations
│   └── measurements/             # Derived statistics
└── requirements.txt              # Python dependencies
```

## Getting Started

1. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

2. Place your data files in the appropriate directories:
   - `selection_criteria_data.csv` → `data/raw/`
   - `COSMOSWeb_mastercatalog_v1.fits` → `data/raw/`

3. Start with the exploration notebook:
   ```bash
   jupyter notebook notebooks/01_data_exploration.ipynb
   ```