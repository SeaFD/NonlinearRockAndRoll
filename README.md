# NonlinearRockAndRoll

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2019b+-blue.svg)](https://www.mathworks.com/products/matlab.html)

## Overview

**NonlinearRockAndRoll** is a MATLAB-based modeling and analysis toolkit for parametric resonance in wave energy converters. This repository provides both hydrodynamic modeling capabilities and a real-time detection system for identifying parametric instabilities in floating structures.

### Key Features

- **Hydrodynamic Model**: Nonlinear hydrodynamic modeling framework with validated coefficients for wave energy converter dynamics
- **Detection System**: Real-time parametric instability detection using advanced signal processing techniques  
- **Validated Coefficients**: Pre-computed hydrodynamic parameters from experimental and CFD validation studies

## Contents

- `HydrodynamicModel.zip`: Complete MATLAB implementation of nonlinear hydrodynamic model with validated coefficients
- `DetectionSystem.zip`: Real-time parametric resonance detection algorithm with example datasets

## Associated Publications

This code supports peer-reviewed research on parametric resonance in wave energy converters. If you use this software in your research, please cite the relevant publications available in the repository documentation.

## Requirements

- MATLAB R2019b or later (tested up to R2024a)
- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox (for detection system)

## Quick Start

### 1. Clone and Extract

```bash
git clone https://github.com/SeaFD/NonlinearRockAndRoll.git
cd NonlinearRockAndRoll
```

### 2. Extract ZIP Files

```bash
unzip HydrodynamicModel.zip
unzip DetectionSystem.zip
```

### 3. Add to MATLAB Path

```matlab
addpath(genpath('HydrodynamicModel'))
addpath(genpath('DetectionSystem'))
```

### 4. Run Example

```matlab
% Load validated hydrodynamic coefficients
load('hydrodynamic_coeffs.mat')

% Set up simulation parameters
params.wave_height = 2.0;  % meters
params.wave_period = 8.0;  % seconds
params.duration = 300;     % seconds

% Run simulation
results = run_hydrodynamic_model(params);

% Visualize results
plot_response(results);
```

## Documentation

Detailed documentation is provided within each MATLAB function. Use MATLAB's `help` command:

```matlab
help run_hydrodynamic_model
help ParametricDetector
```

## Applications

This toolkit has been used for:

- Wave energy converter performance optimization
- Structural health monitoring of floating platforms
- Real-time control system development
- Parametric instability prediction in irregular seas

## Validation

The hydrodynamic model and detection system have been validated against:

- Experimental wave tank data from IST Lisbon facilities
- High-fidelity CFD simulations using OpenFOAM
- Full-scale measurements from operational wave energy devices

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Authors

**Joshua Davidson**
- ORCID: [0000-0001-5966-4272](https://orcid.org/0000-0001-5966-4272)
- GitHub: [@SeaFD](https://github.com/SeaFD)

## Related Projects

- [NWT_potentialFreeSurfaceFoam](https://github.com/SeaFD/NWT_potentialFreeSurfaceFoam) - OpenFOAM numerical wave tank
- [A-NWT](https://github.com/SeaFD/A-NWT) - Axisymmetric numerical wave tank

## Acknowledgments

This work was supported by the Marie Skłodowska-Curie Individual Fellowship (Grant No. 866249).
