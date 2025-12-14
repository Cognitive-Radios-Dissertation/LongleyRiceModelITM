# Longley-Rice (ITM) Propagation Model in MATLAB

A production-grade implementation of the Longley-Rice Irregular Terrain Model (ITM) for Point-to-Point radio propagation prediction.

## Project Structure

This project is organized as follows:

- **`src/`**: Contains the core ITM algorithms and subroutines.
  - `longley_rice_p2p.m`: Main model interface.
  - `hzns.m`, `dlthx.m`, `zlsq1.m`, `qlrpfl.m`: Geometric preparatory functions.
  - `lrprop.m`: Core propagation physics (LOS, Diffraction, Scatter).
  - `avar.m`: Statistical variability calculations.

- **`scripts/`**: Simulation drivers and utility scripts.
  - `run_simulation.m`: Main script to run a path loss simulation over a terrain profile.
  - `create_mock_data.m`: Utility to generate synthetic terrain if real data is missing.

- **`data/`**: Input storage.
  - `X.txt` / `X.04`: Sample terrain profile (Distance vs. Elevation in meters, 10m spacing).

- **`results/`**: Output storage.
  - `pathloss_plot.png`: Generated graphs from the simulation.
  - `pathloss_results.csv`: Numerical results with distance-dependent path loss breakdown.

- **`docs/`**: Documentation and development logs.
  - `geminiExp.md`: Comprehensive implementation plan with self-critique and mathematical formulations.

## Usage

### Prerequisites
- MATLAB (R2018b or later recommended).

### Running a Simulation

1.  Open MATLAB.
2.  Navigate to the `scripts/` directory:
    ```matlab
    cd scripts
    ```
3.  Run the simulation:
    ```matlab
    run_simulation
    ```
4.  The script will:
    -   Load the terrain profile from `data/X.txt` or `data/X.04`.
    -   Calculate path loss for the specified parameters (default: 970 MHz, 52m/2.4m antennas, vertical polarization).
    -   Generate and save a plot to `results/pathloss_plot.png`.
    -   Save numerical results to `results/pathloss_results.csv`.
    -   Display summary statistics in the console.

## Configuration

Default parameters (can be modified in `scripts/run_simulation.m`):
- **Frequency**: 970 MHz
- **Tx Antenna Height**: 52 meters above ground
- **Rx Antenna Height**: 2.4 meters above ground
- **Polarization**: Vertical (pol = 1)
- **Surface Refractivity**: 301 N-units
- **Ground Permittivity**: 15
- **Ground Conductivity**: 0.005 S/m
- **Climate**: Continental Temperate (code 5)
- **Confidence Level**: 50% (median)

## Implementation Features

### Production-Ready Physics
-   **Effective Earth Radius**: K-factor calculation from surface refractivity
-   **Horizon Extraction**: Bidirectional scan with Earth curvature correction
-   **Terrain Irregularity**: Interdecile range of detrended elevations
-   **Diffraction**: Blended knife-edge and smooth-earth (Vogler) with terrain roughness weighting
-   **Troposcatter**: Physics-based forward scatter with climate corrections
-   **Ground Reflection**: Full Fresnel coefficients (permittivity, conductivity, polarization)
-   **Variability**: Climate-dependent time variability and frequency-dependent location variability

### Outputs
The implementation provides:
1. **Total path loss** in dB for complete path
2. **Path loss versus distance** at every terrain profile point
3. **Detailed breakdown** including:
   - Free space loss (L_bf)
   - Reference attenuation (A_ref)
   - Variability adjustments (Y_total)
   - Propagation mode at each point (1=LOS, 2=Diffraction, 3=Scatter)
   - Geometry parameters (horizons, effective heights, terrain irregularity)

### Configuration
Default parameters (configurable via options structure):
-   Operating frequency: 970 MHz
-   Polarization: 1 (Vertical)
-   Transmitter height: 52 meters AGL
-   Receiver height: 2.4 meters AGL
-   Surface refractivity: 301 N-units
-   Ground permittivity: 15
-   Ground conductivity: 0.005 S/m
-   Climate code: 5 (Continental Temperate)

## Validation

The implementation includes strict input validation for:
-   Frequency (20 - 20,000 MHz)
-   Ground Constants (Permittivity 1-100, Conductivity 0.00001-10 S/m)
-   Terrain Profile Length (Min 10 points)
-   Antenna Heights (0.5 - 3000 meters)

## Documentation

For detailed information about the implementation, including:
- Mathematical formulations for all ITM subroutines
- Self-critique analysis and improvement iterations
- Terrain processing pipeline
- Implementation decisions and rationale

See `docs/geminiExp.md` for the comprehensive implementation plan.

## License
[Your License Here]
