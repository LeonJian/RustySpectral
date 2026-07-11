# Installation

## Prerequisites

- Rust 1.70 or later
- `cargo` package manager

## Building from source

```bash
git clone https://github.com/anomalyco/RustySpectral.git
cd RustySpectral
cargo build --release
```

## Adding as a dependency

In your `Cargo.toml`:

```toml
[dependencies]
rustyspectral = "0.1.0"
ndarray = "0.16"
```

## Static Data

RustySpectral requires two sets of static ancillary data:

### Relative Spectral Response Data (RSR)

Standardized HDF5 files containing instrument relative spectral responses
for all supported satellite sensors. Downloaded automatically from Zenodo
when first needed.

```bash
cargo run --bin download_rsr
```

Options:
- `-o /path/to/dir` — destination directory
- `-d` — dry run (don't actually download)
- `-v` — verbose logging

### Atmospheric Correction Look-Up Tables (LUTs)

HDF5 files with pre-computed top-of-atmosphere Rayleigh reflectance for
various aerosol types, standard atmospheres, and viewing geometries.

```bash
# Download all LUTs
cargo run --bin download_atm_correction_luts

# Download specific aerosol types
cargo run --bin download_atm_correction_luts -- -a desert_aerosol,marine_clean_aerosol

# Dry run
cargo run --bin download_atm_correction_luts -- -d
```

Data is stored in `~/.local/share/pyspectral/` by default (platform standard
user data directory).

## Configuration

A YAML configuration file controls paths and download behavior.
Create `pyspectral.yaml` anywhere and point to it:

```bash
export PSP_CONFIG_FILE=/path/to/pyspectral.yaml
```

Example configuration:

```yaml
# Custom data directories
rsr_dir: /data/satellite/rsr
rayleigh_dir: /data/satellite/rayleigh_luts

# Disable automatic downloads (for offline/operational use)
download_from_internet: false

# Optional: per-platform Tb↔Radiance LUT filename overrides
# Meteosat-9-seviri:
#   tb2rad_lut_filename: /path/to/lut/tb2rad_lut_meteosat9_seviri_ir3.9.npz
```

If no `PSP_CONFIG_FILE` is set, defaults are:
- `rsr_dir`: `~/.local/share/pyspectral/`
- `rayleigh_dir`: `~/.local/share/pyspectral/`
- `download_from_internet`: `true`

### Operational (offline) deployment

1. Set `download_from_internet: true` initially, run both download binaries
2. Change to `download_from_internet: false`
3. Deploy with the populated data directories
