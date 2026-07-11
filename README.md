# RustySpectral

Rust port of [pyspectral](https://github.com/pytroll/pyspectral) — satellite sensor
spectral response analysis, solar irradiance computation, and atmospheric corrections
for visible and infrared satellite bands.

[![Crates.io](https://img.shields.io/crates/v/rustyspectral)](https://crates.io/crates/rustyspectral)
[![License](https://img.shields.io/crates/l/rustyspectral)](LICENSE)

## Overview

RustySpectral provides tools for:

- **Spectral response functions** (RSR) — load, query, and convert per-band/per-detector
  relative spectral responses for 56 satellite platforms in HDF5 format
- **Planck blackbody radiation** — radiance ↔ brightness temperature conversion
- **Solar irradiance** — in-band solar flux and spectral irradiance computed by
  convolving the ASTM E-490-00 solar spectrum with instrument RSR curves
- **Near-infrared reflectance** — derive the solar reflectance in the 3.7–3.9 µm
  window channel, removing thermal contribution with optional CO₂ absorption correction
- **Rayleigh scattering correction** — atmospheric correction of shortwave imager
  bands (400–800 nm) using pre-computed LUTs
- **IR limb-cooling correction** — parametric zenith-angle correction
  for "dirty window" infrared bands (DWD method)
- **SEVIRI radiance conversion** — nonlinear regression coefficient support

## Features

| Feature | Module | Status |
|---------|--------|--------|
| Blackbody Planck | `blackbody` | Full — radiance, rad2temp, both λ and ν space |
| Solar spectrum | `solar` | Full — load, interpolate, in-band flux/irradiance |
| Radiance ↔ Tb | `radiance_tb` | Full — RSR integration + SEVIRI regression |
| NIR reflectance | `reflectance` | Full — 3.9µm solar reflectance, CO₂ Rosenfeld |
| IR atm correction | `atm_correction_ir` | Full — DWD parametric zenith correction |
| Rayleigh correction | `rayleigh` | Full — 3D trilinear LUT interpolation, cloud relaxation |
| RSR reader | `rsr_reader` | Full — HDF5 I/O via hdf5-pure (no system HDF5 needed) |
| Raw RSR base | `raw_reader` | Full — customizable raw-reading framework |
| Band names | `bandnames` | Full — 12 sensors (MODIS, VIIRS, SEVIRI, ABI, AHI, etc.) |
| Config system | `config` | Full — YAML config, PSP_CONFIG_FILE env var |
| Data download | `download` | Full — Zenodo HTTP download for RSR and LUT data |
| CLI tools | `bin/` | Full — `download_rsr`, `download_atm_correction_luts` |

## Quick Start

### Prerequisites

Rust 1.70+ with `cargo`.

```bash
git clone https://github.com/.../RustySpectral.git
cd RustySpectral
cargo build --release
```

### Download Ancillary Data

Before using the RSR reader or Rayleigh correction, download the required static data:

```bash
# Download RSR data (HDF5 format response functions)
cargo run --bin download_rsr

# Download atmospheric correction LUTs
cargo run --bin download_atm_correction_luts

# Dry-run to see what would be downloaded
cargo run --bin download_rsr -- --dry_run
```

Static data are stored in `~/.local/share/pyspectral/` by default.
Set `PSP_CONFIG_FILE` to a custom YAML configuration to override paths.

### Configuration

Create a `pyspectral.yaml` file:

```yaml
rsr_dir: /data/rsr
rayleigh_dir: /data/rayleigh_luts
download_from_internet: true
```

Then export:

```bash
export PSP_CONFIG_FILE=/path/to/pyspectral.yaml
```

If no config is found, defaults are used (auto-download enabled, paths under
`~/.local/share/pyspectral/`).

## Usage Examples

### Load Spectral Responses

```rust
use rustyspectral::rsr_reader::RelativeSpectralResponse;

fn main() -> Result<(), String> {
    // Load by platform + instrument (auto-downloads if needed)
    let olci = RelativeSpectralResponse::new(
        Some("Sentinel-3A"),
        Some("olci"),
        None,
    )?;

    println!("Platform: {}", olci.platform_name);
    println!("Bands: {:?}", olci.band_names);

    // Central wavelength of first band
    if let Some(det1) = olci.rsr.get("Oa01").and_then(|d| d.get("det-1")) {
        println!("Central wavelength Oa01: {} µm", det1.central_wavelength);
    }

    // RSR integral
    let integral = olci.integral("Oa01");
    println!("RSR integral Oa01: {:?}", integral);

    Ok(())
}
```

### Solar Irradiance

```rust
use rustyspectral::solar::SolarIrradianceSpectrum;
use rustyspectral::rsr::Rsr;
use ndarray::arr1;

fn main() {
    let mut solar = SolarIrradianceSpectrum::new("data/e490_00a.dat", 0.005);

    // Solar constant
    let sc = solar.solar_constant();
    println!("Solar constant: {:.2} W/m²", sc);
    // Output: Solar constant: ~1365 W/m²

    // In-band solar flux for a 3.9µm-like RSR band
    let wvl = arr1(&[
        3.6124, 3.6164, 3.6265, 3.6364, 3.6465, 3.6565, 3.6664, 3.6765,
        3.6865, 3.6965, 3.7065, 3.7165, 3.7265, 3.7364, 3.7464, 3.7564,
        3.7664, 3.7763, 3.7863, 3.7964, 3.8064, 3.8164, 3.8264, 3.8365,
        3.8463, 3.8564, 3.8664, 3.8764, 3.8865, 3.8965, 3.9064, 3.9165,
        3.9265, 3.9364, 3.9465, 3.9535_f64,
    ]);
    let resp = arr1(&[
        0.01, 0.0118, 0.01987, 0.03226, 0.05028, 0.0849, 0.16645, 0.33792,
        0.59106, 0.81815, 0.96077, 0.92855, 0.86008, 0.8661, 0.87697, 0.85412,
        0.88922, 0.9541, 0.95687, 0.91037, 0.91058, 0.94256, 0.94719, 0.94808,
        1.0, 0.92676, 0.67429, 0.44715, 0.27762, 0.14852, 0.07141, 0.04151,
        0.02925, 0.02085, 0.01414, 0.01_f64,
    ]);
    let rsr = Rsr::new(wvl, resp);
    let flux = solar.inband_solarflux(&rsr, 1.0);
    println!("In-band solar flux: {:.6} W/m²", flux);
    // Output: 2.002928
}
```

### Blackbody Radiation

```rust
use rustyspectral::blackbody::*;

fn main() {
    let wavel = 11e-6;  // 11 µm
    let temp = 300.0;    // 300 K

    let radiance = planck(wavel, temp);
    println!("Radiance at 11µm, 300K: {:.3} W/m²/sr/m", radiance);
    // ~9573177 W/m²/sr/m

    // Round-trip: radiance → temperature
    let t = blackbody_rad2temp(wavel, radiance);
    println!("Round-trip temperature: {:.6} K", t);
    // 300.000000 K

    // Wavenumber space
    let wn = 90909.1;  // cm⁻¹
    let rad_wn = planck_wn(wn, temp);
    let t_wn = blackbody_wn_rad2temp(wn, rad_wn);
    println!("Wavenumber round-trip: {:.6} K", t_wn);
    // 300.000000 K
}
```

### SEVIRI Radiance Conversions

```rust
use rustyspectral::radiance_tb::{SeviriRadTbConverter, self};

fn main() {
    // Use the converter struct (auto-looks up regression parameters)
    let conv = SeviriRadTbConverter::new("Meteosat-9", "IR3.9")
        .expect("Band/platform not found");

    let tb = 300.0;
    let rad = conv.tb2radiance(tb);
    println!("IR3.9 radiance at {tb}K: {:.8} W/m²/sr/(cm⁻¹)⁻¹", rad);
    // ~9.797091e-6

    let tb_back = conv.radiance2tb(rad);
    println!("Round-trip: {:.6} K", tb_back);
    // 300.000000 K

    // Or use standalone functions
    let rad2 = radiance_tb::seviri_tb2radiance(300.0, 2568.832 * 100.0, 0.9954, 3.438);
    let tb2 = radiance_tb::seviri_radiance2tb(rad2, 2568.832 * 100.0, 0.9954, 3.438);
    println!("Standalone round-trip: {:.6} K", tb2);
}
```

### NIR Reflectance (3.9µm channel)

```rust
use rustyspectral::reflectance::ReflectanceCalculator;
use ndarray::arr1;

fn main() {
    // RSR curve for a 3.9µm-like band
    let wvl = ndarray::Array1::linspace(3.6e-6, 3.95e-6, 36);
    let resp = ndarray::Array1::ones(36);

    let calc = ReflectanceCalculator::new("EOS-Aqua", "modis")
        .with_rsr(wvl, resp)
        .with_solar_flux(2.002928);

    let sunz = arr1(&[80.0]);    // sun zenith angle (°)
    let tb_nir = arr1(&[290.0]); // 3.9µm brightness temp (K)
    let tb_thermal = arr1(&[282.0]); // 11µm brightness temp (K)

    let refl = calc.reflectance_from_tbs(&sunz, &tb_nir, &tb_thermal, None);
    println!("Reflectance: {:.6}", refl[0]);
    // ~0.251

    // With CO₂ absorption correction (Rosenfeld method)
    let tb_co2 = arr1(&[270.0]); // 13.4µm CO₂ band
    let refl_co2 = calc.reflectance_from_tbs(&sunz, &tb_nir, &tb_thermal, Some(&tb_co2));
    println!("Reflectance (CO₂ corrected): {:.6}", refl_co2[0]);

    // Emissive part
    let rad3x_t11 = calc.tb2radiance(&tb_thermal);
    let r3x = &refl;
    let rad3x = calc.tb2radiance(&tb_nir);
    let emissive = calc.emissive_part_3x(&rad3x_t11, r3x, &rad3x, false);
    println!("Emissive radiance: {:.6}", emissive[0]);
}
```

### Rayleigh Atmospheric Correction

```rust
use rustyspectral::rayleigh::Rayleigh;
use ndarray::arr1;

fn main() {
    let rcor = Rayleigh::new(
        "GOES-16",          // platform
        "abi",              // sensor
        Some("us_standard"), // atmosphere
        Some("marine_clean_aerosol"), // aerosol type
    );

    let sun_zenith = arr1(&[50.0, 60.0]);
    let sat_zenith = arr1(&[40.0, 50.0]);
    let azidiff = arr1(&[160.0, 160.0]);
    let redband = arr1(&[0.1, 0.15]); // optional cloud relaxation

    // Use wavelength directly (0.64µm → 640nm)
    let refl = rcor.get_reflectance(
        &sun_zenith,
        &sat_zenith,
        &azidiff,
        "0.64",          // wavelength in µm
        Some(&redband),
    );
    println!("Rayleigh reflectance: {:?}", refl);
}
```

### IR Atmospheric Correction (Limb Cooling)

```rust
use rustyspectral::atm_correction_ir::viewzen_corr;
use ndarray::Array2;

fn main() {
    let satz = Array2::from_elem((10, 10), 48.3);
    let tbs = Array2::from_elem((10, 10), 284.5);
    let corrected = viewzen_corr(&tbs, &satz);
    println!("Corrected Tb[0,0]: {:.6} K", corrected[[0, 0]]);
    // ~286 K (depends on exact input)
}
```

## Supported Platforms

RustySpectral supports the same platforms as pyspectral through identical
INSTRUMENTS mapping (56 satellite-platform combinations):

| Category | Sensors |
|----------|---------|
| Geostationary VIS/IR | SEVIRI (8–11), ABI (16–19), AHI (8–9), AGRI (4A–4B), AMI (2A), FCI (12/MTG-I1) |
| Polar VIS/IR | AVHRR/1–3 (TIROS-N through Metop-C), MODIS (Terra, Aqua), VIIRS (NPP, 20–21) |
| Ocean Color | OLCI (3A–3B), SLSTR (3A–3B), COCTS (1C) |
| Landsat | OLI/TIRS (8–9) |
| Sentinel-2 | MSI (2A–2C) |
| Metop-SG | MetImage |
| FY-3 | MERSI-1/2/-3/-RM, VIRR |
| DSCOVR | EPIC |
| Arctica-M | MSU-GSA |
| Electro-L | MSU-GS |
| GEO-KOMPSAT-2B | GOCI-2 |

**Band name dictionaries** are available for all 12 sensor types:
`generic`, `modis`, `seviri`, `viirs`, `avhrr3`, `abi`, `agri`, `ahi`, `ami`, `fci`, `slstr`, `vii`.

## Physical Constants

All constants match pyspectral exactly:

| Constant | Value | Unit |
|----------|-------|------|
| Planck constant `h` | 6.62606957e-34 | J·s |
| Boltzmann constant `k` | 1.3806488e-23 | J/K |
| Speed of light `c` | 2.99792458e8 | m/s |
| Terminator limit | 85.0 | degrees |
| RSR data version | v1.6.1 | — |

## Key Reference Values

Cross-validated against pyspectral Python:

| Measurement | Expected Value | Tolerance |
|-------------|---------------|-----------|
| Planck(11µm, 300K) | 9573176.935507433 | 1e-4 |
| Solar constant | ~1365 W/m² | range |
| inband_solarflux (MODIS B20) | 2.002927627 | 1e-3 |
| SEVIRI IR3.9 tb2radiance(300K) | 9.797091e-6 | 1e-8 |
| get_central_wave (MODIS B20) | 3.780281935 µm | 1e-6 |
| FWHM (ch3 band) | 0.06581801 | 1e-5 |
| viewzen_corr max diff | < 1e-3 | — |

## API Reference

### Core Modules

| Module | Structs | Key Functions |
|--------|---------|--------------|
| `blackbody` | — | `planck`, `planck_wn`, `blackbody_rad2temp`, `blackbody_wn_rad2temp`, `planck_array_wavelength`, `planck_array_wn`, `blackbody_rad2temp_array`, `blackbody_wn_rad2temp_array` |
| `solar` | `SolarIrradianceSpectrum` | `solar_constant`, `inband_solarflux`, `inband_solarirradiance`, `interpolate`, `set_wavespace_wavenumber` |
| `radiance_tb` | `RadTbConverter`, `SeviriRadTbConverter` | `radiance2tb`, `tb2radiance_simple`, `tb2radiance_array`, `tb2radiance_normalized`, `make_tb2rad_lut`, `seviri_radiance2tb`, `seviri_tb2radiance` |
| `reflectance` | `ReflectanceCalculator` | `reflectance_from_tbs`, `solar_radiance`, `tb2radiance`, `emissive_part`, `emissive_part_3x`, `derive_rad39_corr` |
| `atm_correction_ir` | `AtmosphericalCorrection` | `viewzen_corr` |
| `rayleigh` | `Rayleigh`, `RayleighConfigBase` | `get_reflectance`, `reduce_rayleigh_highzenith`, `normalize_sensor`, `clip_angles_inside_coordinate_range`, `trilinear_interpolate`, `rayleigh_interpolate_by_angles`, `check_and_download` |

### Data & Infrastructure

| Module | Structs | Key Functions |
|--------|---------|--------------|
| `rsr` | `Rsr`, `PerDetectorRsr`, `DetectorRsr` | `to_wavenumber`, `from_detectors` |
| `rsr_reader` | `RelativeSpectralResponse` | `integral`, `convert`, `resolve_band`, `get_bandname_from_wavelength`, `load_rsr_info_from_file`, `check_and_download` |
| `raw_reader` | `InstrumentRSR` | `get_options_from_config`, `get_bandfilenames` |
| `bandnames` | — | `get_bandnames` → 12 sensor dictionaries |
| `config` | `Config` | `get_config`, `recursive_dict_update` |
| `download` | — | `download_rsr`, `download_luts`, `get_rayleigh_lut_dir` |
| `utils` | (types: `RsrData`, `InstrumentValue`, `AtmCorrectionVersion`) | `trapezoid`, `get_central_wave`, `convert2wavenumber_rsr`, `sort_data`, `get_fullwidth_halfmax`, `get_bounds_integrated_energy`, `get_wave_range`, `get_bandname_from_wavelength`, `are_instruments_identical`, `check_and_adjust_instrument_name` |

### Global Constants (in `utils`)

| Constant | Description |
|----------|-------------|
| `INSTRUMENTS` | 56 platform→sensor mappings |
| `AEROSOL_TYPES` | 11 aerosol distribution types |
| `ATMOSPHERES` | 6 standard atmospheres |
| `ATM_CORRECTION_LUT_VERSION` | LUT versions for all aerosol types |
| `HTTPS_RAYLEIGH_LUTS` | Zenodo download URLs |
| `HTTP_PYSPECTRAL_RSR` | RSR data Zenodo URL |
| `RSR_DATA_VERSION` | `"v1.6.1"` |
| `WAVE_LENGTH` / `WAVE_NUMBER` | Wavespace identifiers |

## Project Structure

```
src/
├── lib.rs              # Crate root
├── bin/
│   ├── download_rsr.rs
│   └── download_atm_correction_luts.rs
├── blackbody.rs        # Planck radiation
├── solar.rs            # Solar irradiance spectrum
├── radiance_tb.rs      # Radiance ↔ brightness temperature
├── reflectance.rs      # NIR reflectance (3.9µm)
├── atm_correction_ir.rs # IR limb cooling
├── rayleigh.rs         # Rayleigh scattering correction
├── rsr.rs              # RSR data structures
├── rsr_reader.rs       # HDF5 RSR reader
├── raw_reader.rs       # Raw RSR reader base
├── config.rs           # YAML configuration
├── download.rs         # HTTP data download
├── utils.rs            # Utilities + global constants
└── bandnames.rs        # Sensor band name dictionaries
```

## Differences from pyspectral

- **Pure Rust** — no system HDF5 library required (uses `hdf5-pure`)
- **No dask/xarray** — ndarray is the only array backend
- **No masked arrays** — NaN is used for invalid values
- **Scalar-only Planck** — separate `planck_array_*` functions for array inputs
- **`RsrData` struct** — replaces Python's nested dicts for type safety
- **Builder pattern** — `ReflectanceCalculator` uses `with_*()` methods instead of a large constructor
- **No plotting** — plot with external tools
- **No `convert2hdf5`** — raw-to-internal format conversion not yet implemented
- **Atmosphere names** — normalized to underscores (`us_standard` instead of `us-standard`)

## Testing

```bash
cargo test
# 113 tests, zero warnings
```
