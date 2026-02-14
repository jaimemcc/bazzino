# Archived Notebooks

This folder contains notebooks that have been archived as part of project organization efforts. Notebooks are archived when they represent completed development stages, exploratory work, or have been superseded by more recent implementations.

## Contents

### Assembly Notebooks (Data Processing Prototypes)
- **assemble_data.ipynb** - Initial proof-of-concept for reading TDT tank files and extracting photometry snips
- **assemble_dlc_data.ipynb** - Initial extraction of DeepLabCut behavioral data
- **Status:** Functionality has been refactored into `src/assemble_all_data.py` which serves as the canonical data assembly pipeline
- **Preserved for:** Reference to the development approach and original implementation details

### cluster_analysis_exploratory/ (Alternative Analysis Approaches)
Exploratory clustering and analysis notebooks representing alternative methods that were tested:
- **spectral_cluster_averaged_trials.ipynb** - Spectral clustering with trial-averaged snips
- **spectral_clustering_angular_velocity.ipynb** - Clustering analysis on angular velocity data (behavioral, not photometry)
- **cluster_analysis_averaged_trials.ipynb** - Analysis of 7-cluster model with averaged trials
- **clusters_diff_from_zero.ipynb** - Statistical testing of cluster responses
- **Status:** These represent exploratory work testing alternative clustering approaches
- **Current Production:** `spectral_clustering_all_trials.ipynb` (in main notebooks folder) implements the standard clustering approach with all trials

### Previous Implementation Attempts
- **fig1_and_2_sexes_divided.ipynb** - Figure generation split by sex (analysis approach changed)
- **figs_for_ang_vel_sexes_divided.ipynb** - Angular velocity figures split by sex
- **add_sex_to_xarrays.ipynb** - Sex data integration (analysis direction changed)
- **Status:** Sex-based analysis variations; current pipeline uses aggregated analysis

## Notes for Future Work

1. **Data Assembly Pipeline**: If you need to understand how data was originally assembled, or if you want to revert to an earlier approach, check assemble_data.ipynb and assemble_dlc_data.ipynb
2. **Clustering Methods**: If you want to revisit alternative clustering approaches, see cluster_analysis_exploratory/ 
3. **When to Use**: Archived notebooks are useful for:
   - Understanding the evolution of the analysis
   - Reverting to earlier approaches if needed
   - Documenting what alternatives were tried and why they weren't selected

## Archive Organization Strategy

- Exploratory work and alternative methods are preserved in dedicated subdirectories
- This maintains a clean active workspace while preserving research history
- Main notebooks/ folder now contains only the active, production-level analysis pipeline
