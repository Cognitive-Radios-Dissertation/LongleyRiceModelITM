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
  - `X.04`: Sample terrain profile (Distance vs. Elevation).

- **`results/`**: Output storage.
  - `pathloss_plot.png`: Generated graphs from the simulation.

- **`docs/`**: Documentation and development logs.

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
    -   Load or create the `data/X.04` terrain profile.
    -   Calculate path loss for the specified parameters (default: 970 MHz).
    -   Generate and save a plot to `results/pathloss_plot.png`.

## Validation

The implementation includes strict input validation for:
-   Frequency (20 - 20,000 MHz)
-   Ground Constants (Permittivity and Conductivity)
-   Terrain Profile Length (Min 10 points)

## License
[Your License Here]
