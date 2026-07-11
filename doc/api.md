# API Reference

## `blackbody` — Planck blackbody radiation

| Function | Signature | Description |
|----------|-----------|-------------|
| `planck` | `(wave: f64, temperature: f64) → f64` | Spectral radiance in wavelength space (W·m⁻²·sr⁻¹·m⁻¹) |
| `planck_wn` | `(wavenumber: f64, temperature: f64) → f64` | Spectral radiance in wavenumber space |
| `blackbody` | `(wave: f64, temperature: f64) → f64` | Alias for `planck` |
| `blackbody_wn` | `(wavenumber: f64, temperature: f64) → f64` | Alias for `planck_wn` |
| `planck_array_wavelength` | `(wave: f64, temperature: &Array2<f64>) → Array2<f64>` | Vectorized Planck |
| `planck_array_wn` | `(wavenumber: f64, temperature: &Array2<f64>) → Array2<f64>` | Vectorized Planck (wn) |
| `blackbody_rad2temp` | `(wavelength: f64, radiance: f64) → f64` | Inverse Planck: radiance → BT |
| `blackbody_wn_rad2temp` | `(wavenumber: f64, radiance: f64) → f64` | Inverse Planck (wn): radiance → BT |
| `blackbody_rad2temp_array` | `(wavelength: f64, radiance: &Array2<f64>) → Array2<f64>` | Vectorized inverse |
| `blackbody_wn_rad2temp_array` | `(wavenumber: f64, radiance: &Array2<f64>) → Array2<f64>` | Vectorized inverse (wn) |

**Constants:** `H_PLANCK`, `K_BOLTZMANN`, `C_SPEED`

**Edge cases:** zero temperature → NaN, zero/negative radiance → NaN

## `solar` — Solar irradiance spectrum

### `SolarIrradianceSpectrum`

| Method | Signature | Description |
|--------|-----------|-------------|
| `new` | `(filename: &str, dlambda: f64) → Self` | Load spectrum from file |
| `solar_constant` | `(&self) → f64` | Total integrated irradiance (W/m²) |
| `inband_solarflux` | `(&mut self, rsr: &Rsr, scale: f64) → f64` | In-band flux convolved with RSR |
| `inband_solarirradiance` | `(&mut self, rsr: &Rsr, scale: f64) → f64` | Normalized spectral irradiance |
| `interpolate` | `(&mut self, dlambda: f64, range: Option<(f64, f64)>)` | Interpolate onto regular grid |
| `set_wavespace_wavenumber` | `(&mut self)` | Convert to wavenumber space |

**Fields:** `wavelength: Array1<f64>`, `irradiance: Array1<f64>`,
`ipol_wavelength: Option<Array1<f64>>`, `ipol_irradiance: Option<Array1<f64>>`,
`wavenumber: Option<Array1<f64>>`

## `radiance_tb` — Radiance ↔ brightness temperature

### Standalone Functions

| Function | Signature | Description |
|----------|-----------|-------------|
| `radiance2tb` | `(radiance: f64, wavelength: f64) → f64` | Radiance → BT via inverse Planck |
| `tb2radiance_simple` | `(tb: f64, wavelength: &Array1<f64>, response: &Array1<f64>) → f64` | BT → integrated radiance |
| `tb2radiance_array` | `(tb: &Array1<f64>, wavelength: &Array1<f64>, response: &Array1<f64>) → Array1<f64>` | Vectorized BT → radiance |
| `tb2radiance_normalized` | `(tb: f64, wavelength: &Array1<f64>, response: &Array1<f64>) → f64` | BT → normalized radiance |
| `make_tb2rad_lut` | `(wavelength: &Array1<f64>, response: &Array1<f64>, tb_resolution: f64) → (Array1<f64>, Array1<f64>)` | Generate TB→Radiance LUT |
| `seviri_radiance2tb` | `(radiance: f64, vc: f64, alpha: f64, beta: f64) → f64` | SEVIRI radiance → TB |
| `seviri_tb2radiance` | `(tb: f64, vc: f64, alpha: f64, beta: f64) → f64` | SEVIRI TB → radiance |
| `get_seviri_params` | `() → HashMap<&str, HashMap<&str, (f64, f64, f64)>>` | All SEVIRI regression parameters |

**Constants:** `TB_MIN = 150.0`, `TB_MAX = 360.0`, `EPSILON = 0.01`,
`SEVIRI` (lazy static map of all regression parameters)

### `RadTbConverter`

| Method | Signature | Description |
|--------|-----------|-------------|
| `new` | `(platform_name, instrument, band, wavelength, response) → Self` | Create converter with RSR |
| `with_detector` | `(self, detector: &str) → Self` | Set detector name |
| `tb2radiance` | `(&self, tb: &Array1<f64>, normalized: bool) → Array1<f64>` | TB → radiance via RSR |
| `radiance2tb` | `(&self, rad: &Array1<f64>) → Array1<f64>` | Radiance → BT |
| `make_tb2rad_lut` | `(&self, tb_resolution: f64, normalized: bool) → (Array1<f64>, Array1<f64>)` | Generate LUT |

**Fields:** `platform_name`, `instrument`, `band`, `wavelength`, `response`,
`central_wavelength`, `rsr_integral`, `detector`

### `SeviriRadTbConverter`

| Method | Signature | Description |
|--------|-----------|-------------|
| `new` | `(platform_name: &str, band: &str) → Option<Self>` | Auto-lookup regression params |
| `radiance2tb` | `(&self, rad: f64) → f64` | SEVIRI radiance → TB |
| `tb2radiance` | `(&self, tb: f64) → f64` | SEVIRI TB → radiance |

**Fields:** `platform_name`, `band`, `vc`, `alpha`, `beta`

## `reflectance` — NIR reflectance

### `ReflectanceCalculator`

| Method | Signature | Description |
|--------|-----------|-------------|
| `new` | `(platform_name: &str, instrument: &str) → Self` | Create with defaults |
| `with_rsr` | `(self, wavelength: Array1<f64>, response: Array1<f64>) → Self` | Set RSR curve |
| `with_solar_flux` | `(self, solar_flux: f64) → Self` | Set in-band solar flux |
| `with_central_wavelength` | `(self, cw: f64) → Self` | Override central wavelength |
| `with_sunz_threshold` | `(self, threshold: f64) → Self` | Sun zenith threshold |
| `with_masking_limit` | `(self, limit: Option<f64>) → Self` | Masking limit (None = no mask) |
| `reflectance_from_tbs` | `(&self, sun_zenith, tb_near_ir, tb_thermal, tb_ir_co2) → Array1<f64>` | Main reflectance computation |
| `solar_radiance` | `(&self, sun_zenith: &Array1<f64>) → Array1<f64>` | Compute solar radiance |
| `tb2radiance` | `(&self, tb: &Array1<f64>) → Array1<f64>` | BT → radiance (RSR or central λ) |
| `derive_rad39_corr` | `(&self, bt11: &Array1<f64>, bt13: &Array1<f64>) → Array1<f64>` | CO₂ correction factor |
| `emissive_part` | `(&self, tb_nir: &Array1<f64>, tb_thermal: Option<&Array1<f64>>) → Array1<f64>` | Emissive radiance |
| `emissive_part_3x` | `(&self, rad3x_t11, r3x, rad3x, tb: bool) → Array1<f64>` | Full emissive part (stateful) |

**Fields:** `platform_name`, `instrument`, `solar_flux`, `central_wavelength`,
`sunz_threshold`, `masking_limit`

**Constants:** `TERMINATOR_LIMIT = 85.0`, `TB_MIN = 150.0`, `TB_MAX = 360.0`

## `atm_correction_ir` — IR limb cooling

### `AtmosphericalCorrection`

| Method | Signature | Description |
|--------|-----------|-------------|
| `new` | `(platform_name: &str, sensor: &str) → Self` | Create correction instance |
| `get_correction` | `(&self, sat_zenith: &Array2<f64>, bandname: &str, data: &Array2<f64>) → Array2<f64>` | Apply correction |

### `viewzen_corr`

| Function | Signature | Description |
|----------|-----------|-------------|
| `viewzen_corr` | `(data: &Array2<f64>, view_zen: &Array2<f64>) → Array2<f64>` | DWD parametric correction |

**Inner functions** (identical to pyspectral): `tau0(T)`, `tau(T)`, `delta(z)`, `ratio(v, v0, v_ref)`

## `rayleigh` — Rayleigh scattering correction

### `Rayleigh`

| Method | Signature | Description |
|--------|-----------|-------------|
| `new` | `(platform_name, sensor, atmosphere, aerosol_type) → Self` | Initialize with LUT |
| `get_reflectance` | `(&self, sun_zenith, sat_zenith, azidiff, band_or_wavelength, redband) → Array1<f64>` | Main correction (0–100%) |

### Standalone Functions

| Function | Signature | Description |
|----------|-----------|-------------|
| `clip_angles_inside_coordinate_range` | `(zenith: &Array1<f64>, secant_max: f64) → Array1<f64>` | Clip to LUT range |
| `clip_angles_inside_coordinate_range_scalar` | `(zenith: f64, secant_max: f64) → f64` | Clip scalar angle |
| `reduce_rayleigh_highzenith` | `(zenith, rayref, thresh, maxzen, strength) → Array1<f64>` | High-zenith reduction |
| `get_wavelength_index_and_factor` | `(wvl_coord: &Array1<f64>, wvl: f64) → (usize, f64)` | LUT wavelength indexing |
| `get_wavelength_adjusted_lut` | `(refl: &Array4<f64>, wvl_coord: &Array1<f64>, wvl: f64) → Array3<f64>` | Wavelength interpolation in 4D LUT |
| `trilinear_interpolate` | `(grid, sunz_sec, azidiff, satz_sec, coords...) → f64` | 3D multilinear interpolation |
| `rayleigh_interpolate_by_angles` | `(sunz, satz, azid, refl, coords...) → Array1<f64>` | Full interpolation loop |
| `normalize_sensor` | `(platform_name: &str, sensor: &str) → String` | Validate sensor name |
| `get_reflectance_lut_from_file` | `(path: &Path) → Result<(Array1<f64>, Array1<f64>, Array1<f64>), String>` | Read LUT coordinates |
| `read_reflectance_lut_4d` | `(path: &Path) → Result<Array4<f64>, String>` | Read 4D reflectance data |
| `read_wavelength_lut_coord` | `(path: &Path) → Result<Array1<f64>, String>` | Read wavelength coordinates |
| `check_and_download` | `(dry_run: bool, aerosol_types: Option<&[String]>)` | Check/update LUTs |

### `RayleighConfigBase`

| Method | Signature | Description |
|--------|-----------|-------------|
| `new` | `(aerosol_type: &str, atm_type: &str) → Self` | Version-check initialization |

**Fields:** `aerosol_type`, `atm_type`, `do_download`, `lutfiles_version_uptodate`

## `rsr` — RSR data structures

| Struct | Fields | Methods |
|--------|--------|---------|
| `Rsr` | `wavelength: Array1<f64>`, `response: Array1<f64>`, `central_wavelength: f64` | `new`, `from_detectors`, `to_wavenumber` |
| `PerDetectorRsr` | `name: String`, `wavelength`, `response`, `central_wavelength`, `wavenumber: Option<Array1<f64>>` | `new` |
| `DetectorRsr` | `detectors: Vec<PerDetectorRsr>` | — |

## `rsr_reader` — RSR HDF5 reader

### `RelativeSpectralResponse`

| Method | Signature | Description |
|--------|-----------|-------------|
| `new` | `(platform_name, instrument, filename) → Result<Self, String>` | Load RSR from HDF5 |
| `integral` | `(&self, band_name: &str) → HashMap<String, f64>` | RSR integrals per detector |
| `convert` | `(&mut self)` | λ → ν conversion |
| `get_bandname_from_wavelength` | `(&self, wavel, epsilon, multiple_bands) → Option<Vec<String>>` | Find band by wavelength |
| `resolve_band` | `(&self, key: &str) → Option<String>` | Sensor→generic band name fallback |

**Fields:** `platform_name`, `instrument`, `description`, `band_names`,
`rsr: RSRDict`, `unit`, `si_scale`, `filename`

### Functions

| Function | Signature | Description |
|----------|-----------|-------------|
| `load_rsr_info_from_file` | `(filename: &Path) → Result<RsrFileInfo, String>` | Read HDF5 into struct |
| `check_and_download` | `(dest_dir: Option<&Path>, dry_run: bool)` | Version check + download |

## `raw_reader` — Raw RSR base

### `InstrumentRSR`

| Method | Signature | Description |
|--------|-----------|-------------|
| `new` | `(bandname: &str, platform_name: &str, bandnames: &[String]) → Self` | Create reader |
| `get_options_from_config` | `(&mut self)` | Load path config |
| `get_bandfilenames` | `(&mut self)` | Resolve band files |

**Fields:** `platform_name`, `instrument`, `bandname`, `bandnames`, `filenames`,
`output_dir`, `path`, `filename`

## `bandnames` — Band name dictionaries

| Function | Signature | Description |
|----------|-----------|-------------|
| `get_bandnames` | `() → HashMap<&str, BandNames>` | All 12 sensor dictionaries |

**Type alias:** `BandNames = HashMap<&str, &str>`

## `config` — Configuration system

| Function | Signature | Description |
|----------|-----------|-------------|
| `get_config` | `(config_file: Option<&Path>) → Config` | Read YAML config |
| `recursive_dict_update` | `(base: Value, update: &Value) → Value` | Deep merge YAML values |

### `Config`

**Fields:** `rsr_dir: PathBuf`, `rayleigh_dir: PathBuf`,
`download_from_internet: bool`, `raw: HashMap<String, Value>`

## `download` — Data download

| Function | Signature | Description |
|----------|-----------|-------------|
| `download_rsr` | `(dest_dir: Option<&Path>, dry_run: bool) → io::Result<()>` | Download RSR tarball |
| `download_luts` | `(aerosol_types: Option<&[String]>, dry_run: bool) → io::Result<()>` | Download LUT tarballs |
| `get_rayleigh_lut_dir` | `(config: &Config, aerosol_type: &str) → PathBuf` | LUT directory path |

## `utils` — Utilities and constants

### Core Functions

| Function | Signature | Description |
|----------|-----------|-------------|
| `trapezoid` | `(y: &Array1<f64>, x: &Array1<f64>) → f64` | Trapezoidal integration |
| `get_central_wave` | `(wav, resp, weight: &Array1<f64>) → f64` | Weighted central wavelength |
| `convert2wavenumber_rsr` | `(wavelength, response) → (Array1<f64>, Array1<f64>)` | λ → ν conversion |
| `sort_data` | `(x, y) → (Array1<f64>, Array1<f64>)` | Sort + deduplicate |
| `get_fullwidth_halfmax` | `(rsp, wvl) → f64` | FWHM bandwidth |
| `get_bounds_integrated_energy` | `(rsp, wvl, ener_perc_lim) → (f64, f64)` | Energy bounds |
| `get_wave_range` | `(wvl, resp, threshold) → (min, cwl, max)` | Wavelength range |
| `get_bandname_from_wavelength` | `(sensor, wavel, rsr, epsilon, multi) → Option<Vec<String>>` | Band name lookup |
| `are_instruments_identical` | `(name1, name2) → bool` | Instrument name comparison |
| `check_and_adjust_instrument_name` | `(platform, instrument) → String` | Normalize/validate |
| `get_instruments` | `() → HashMap<&str, InstrumentValue>` | All INSTRUMENTS entries |
| `get_atm_correction_lut_version` | `() → HashMap<&str, AtmCorrectionVersion>` | LUT version info |
| `get_https_rayleigh_luts` | `() → HashMap<&str, &str>` | LUT download URLs |
| `get_rayleigh_lut_dir` | `(base_dir: &PathBuf, aerosol_type: &str) → PathBuf` | LUT directory path |

### Types

| Type | Description |
|------|-------------|
| `RsrData` | Struct: `wavelength`, `response`, `central_wavelength` |
| `InstrumentValue` | Enum: `Single(String)` or `List(Vec<String>)` |
| `AtmCorrectionVersion` | Struct: `version: &'static str`, `filename: &'static str` |

### Constants

| Constant | Type | Value |
|----------|------|-------|
| `WAVE_LENGTH` | `&str` | `"wavelength"` |
| `WAVE_NUMBER` | `&str` | `"wavenumber"` |
| `INSTRUMENTS` | Lazy<HashMap> | 56 platform→sensor mappings |
| `AEROSOL_TYPES` | `&[&str]` | 11 aerosol types |
| `ATMOSPHERES` | `&[(&str, usize)]` | 6 atmospheres |
| `ATM_CORRECTION_LUT_VERSION` | Lazy<HashMap> | 11 version entries |
| `HTTPS_RAYLEIGH_LUTS` | Lazy<HashMap> | 11 download URLs |
| `HTTP_PYSPECTRAL_RSR` | `&str` | Zenodo RSR URL |
| `RSR_DATA_VERSION` | `&str` | `"v1.6.1"` |
| `RSR_DATA_VERSION_FILENAME` | `&str` | `"PYSPECTRAL_RSR_VERSION"` |
