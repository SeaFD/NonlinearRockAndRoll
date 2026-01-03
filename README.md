# NonlinearRockAndRoll

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2019b+-blue.svg)](https://www.mathworks.com/products/matlab.html)

## Overview

**NonlinearRockAndRoll** is a MATLAB-based modeling and analysis toolkit for parametric resonance in wave energy converters. This repository provides both hydrodynamic modeling capabilities and a real-time detection system for identifying parametric instabilities in floating structures.

### Key Features

- **Hydrodynamic Model**: Nonlinear hydrodynamic modeling framework with validated coefficients for wave energy converter dynamics
- **Detection System**: Real-time parametric instability detection using advanced signal processing techniques  
- **Validated Coefficients**: Pre-computed hydrodynamic parameters from experimental validation studies


## Associated Publication

This model and analysis toolkit is described in detail in the following publication:

**Davidson, J., & Kalmar-Nagy, T.** (2020). A real-time detection system for the onset of parametric resonance in wave energy converters. *Journal of Marine Science and Engineering*, 8(10), 819. [https://doi.org/10.3390/jmse8100819](https://doi.org/10.3390/jmse8100819)

**If you use this model and analysis toolkit in your research, please cite the above publication.**

## Requirements

- MATLAB R2019b or later (tested up to R2024a)
- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox (for detection system)

## Applications

This toolkit has been used for:

- Wave energy converter performance optimization
- Real-time control system development
- Parametric instability prediction in irregular seas

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Authors

**Josh Davidson**
- ORCID: [0000-0001-5966-4272](https://orcid.org/0000-0001-5966-4272)
- GitHub: [@SeaFD](https://github.com/SeaFD)

## Related Projects

- [NWT_potentialFreeSurfaceFoam](https://github.com/SeaFD/NWT_potentialFreeSurfaceFoam) - OpenFOAM numerical wave tank
- [A-NWT](https://github.com/SeaFD/A-NWT) - Axisymmetric numerical wave tank

## Acknowledgments

This work was supported by the Marie Skłodowska-Curie Individual Fellowship (Grant No. 866249).
