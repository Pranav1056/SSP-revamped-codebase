# SSP orbit determination

Python implementations of orbit determination, orbital elements, astrometry, photometry and orbit visualisation, developed from astronomy research exercises.

## Install and test

```bash
python -m pip install -e ".[dev]"
python -m pytest
```

For FITS photometry, install `.[photometry]`. For VPython visualisation, install `.[visualization]`.

## Run

```bash
python scripts/run_orbit_determination.py observations.txt
python scripts/run_plate_reduction.py examples/pranav_data.txt 1292.809093 1473.2126
python scripts/visualize_orbit.py elements.txt
```

The orbit-determination workflow expects three observations of nine whitespace-separated fields each (date, time, RA, Dec, three Sun-vector components). The included `JiInput.txt` belongs to earlier exercises and is not suitable for that workflow.

## Organisation

- `src/ssp/orbital_elements.py`: orbital elements.
- `src/ssp/kepler.py`: Kepler and Lagrange calculations.
- `src/ssp/coordinates.py`: Julian dates, coordinates and ephemerides.
- `src/ssp/astrometry.py`: least-squares plate reduction.
- `src/ssp/photometry.py`: aperture photometry.
- `scripts/`: runnable command-line scripts.
- `tests/`: basic numerical tests.
- `examples/`: historical sample data.
- `archive/`: original scripts and ZIP archives.

The original mathematical implementations have been preserved with minimal changes. Scientific accuracy and convergence for all input conditions have not been independently validated. See `archive/ZIP_INVENTORY.md` for archived ZIP contents.
