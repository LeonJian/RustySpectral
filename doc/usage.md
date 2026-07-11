# Usage

## Spectral Response Functions

### Loading RSR data

```rust
use rustyspectral::rsr_reader::RelativeSpectralResponse;

// By platform + instrument
let olci = RelativeSpectralResponse::new(
    Some("Sentinel-3A"), Some("olci"), None
)?;

// By explicit HDF5 filename
let modis = RelativeSpectralResponse::new(
    None, None, Some(std::path::Path::new("/data/rsr_modis_EOS-Aqua.h5"))
)?;
```

### Accessing RSR data

```rust
// List all band names
println!("{:?}", olci.band_names);

// Get per-detector response data
if let Some(detectors) = olci.rsr.get("Oa01") {
    if let Some(det1) = detectors.get("det-1") {
        println!("Wavelengths: {:?}", det1.wavelength);   // µm
        println!("Responses: {:?}", det1.response);
        println!("Central wavelength: {} µm", det1.central_wavelength);
    }
}

// Band name resolution (sensor → generic fallback)
let resolved = olci.resolve_band("VIS006");
```

### RSR operations

```rust
// Band integral
let integral_map = olci.integral("Oa01");
for (det, val) in &integral_map {
    println!("Integral {}: {}", det, val);
}

// Convert to wavenumber space
olci.convert();
println!("Unit: {}", olci.unit);           // "cm-1"
println!("Scale: {}", olci.si_scale);      // 100.0

// Find band by wavelength
let bands = olci.get_bandname_from_wavelength(0.67, 0.05, false);
println!("Band at 0.67µm: {:?}", bands);
```

### Low-level RSR structs (no HDF5)

```rust
use rustyspectral::rsr::{Rsr, PerDetectorRsr, DetectorRsr};
use ndarray::arr1;

// Single detector
let wvl = arr1(&[3.6, 3.7, 3.8, 3.9, 4.0_f64]);
let resp = arr1(&[0.0, 0.5, 1.0, 0.5, 0.0_f64]);
let rsr = Rsr::new(wvl, resp);
println!("Central wavelength: {} µm", rsr.central_wavelength);

// Convert to wavenumber
let (wn, wresp) = rsr.to_wavenumber();

// Multi-detector
let det1 = PerDetectorRsr::new("det-1".into(), wvl1, resp1);
let det2 = PerDetectorRsr::new("det-2".into(), wvl2, resp2);
let multi = DetectorRsr { detectors: vec![det1, det2] };
```

## Solar Irradiance

```rust
use rustyspectral::solar::SolarIrradianceSpectrum;
use rustyspectral::rsr::Rsr;

// Load the built-in ASTM E-490-00 spectrum
let mut solar = SolarIrradianceSpectrum::new("data/e490_00a.dat", 0.005);

// Total solar constant
let sc = solar.solar_constant();
println!("Solar constant: {:.2} W/m²", sc);  // ~1365

// Interpolate onto regular grid
solar.interpolate(0.001, Some((0.200, 0.240)));
if let Some(wvl) = &solar.ipol_wavelength {
    println!("Interpolated grid: {} points, {:.3}–{:.3} µm",
        wvl.len(), wvl[0], wvl[wvl.len()-1]);
}

// In-band solar flux (integrated over RSR)
let rsr = Rsr::new(wavelengths, responses);
let flux = solar.inband_solarflux(&rsr, 1.0);
println!("In-band flux: {:.6} W/m²", flux);  // 2.002928 for MODIS B20

// In-band spectral irradiance (normalized by RSR integral)
let irradiance = solar.inband_solarirradiance(&rsr, 1.0);

// Wavenumber-space
solar.set_wavespace_wavenumber();
let sc_wn = solar.solar_constant();
```

## Blackbody Radiation

```rust
use rustyspectral::blackbody::*;

// Planck function (wavelength space)
let wavel = 11e-6;   // 11 µm in meters
let temp = 300.0;    // Kelvin
let rad = planck(wavel, temp);

// Planck function (wavenumber space)
let wn = 90909.1;    // cm⁻¹ at 11µm
let rad_wn = planck_wn(wn, temp);

// Inverse: radiance → brightness temperature
let t = blackbody_rad2temp(wavel, rad);       // wavelength
let t_wn = blackbody_wn_rad2temp(wn, rad_wn); // wavenumber
// Both return ~300.0 K

// Convenience aliases
let rad1 = blackbody(wavel, temp);
let rad2 = blackbody_wn(wn, temp);

// Array temperatures
use ndarray::arr2;
let temps = arr2(&[[300.0, 301.0], [299.0, 298.0]]);
let rad_arr = planck_array_wavelength(wavel, &temps);  // 2x2 output
let t_arr = blackbody_rad2temp_array(wavel, &rad_arr);

// Edge cases (safe)
planck(wavel, 0.0);              // → NaN
blackbody_rad2temp(wavel, 0.0);  // → NaN
```

## Radiance ↔ Brightness Temperature

### RSR-based conversion

```rust
use rustyspectral::radiance_tb::*;
use ndarray::arr1;

// Setup an RSR curve
let wavelength = arr1(&[3.6e-6, 3.7e-6, 3.8e-6, 3.9e-6, 4.0e-6]);
let response = arr1(&[0.1, 0.5, 1.0, 0.5, 0.1]);

// TB → Radiance (via RSR integration + Planck)
let rad = tb2radiance_normalized(300.0, &wavelength, &response);

// Radiance → TB (via inverse Planck at central wavelength)
let central_wl = 3.8e-6;
let tb = radiance2tb(rad, central_wl);

// Array inputs
let tbs = arr1(&[200.0, 270.0, 300.0, 350.0]);
let rads = tb2radiance_array(&tbs, &wavelength, &response);

// Generate a TB→Radiance lookup table
let (lut_tb, lut_rad) = make_tb2rad_lut(&wavelength, &response, 0.1);
// lut spans 150–360K at 0.1K resolution
```

### SEVIRI regression conversion

```rust
use rustyspectral::radiance_tb::*;

// Standalone functions
let vc = 2568.832 * 100.0;  // central wavenumber in m⁻¹
let alpha = 0.9954;
let beta = 3.438;

let rad = seviri_tb2radiance(300.0, vc, alpha, beta);
let tb = seviri_radiance2tb(rad, vc, alpha, beta);
// Round-trip: tb ≈ 300.0 ± 1e-4

// Or use the convenience struct
let conv = SeviriRadTbConverter::new("Meteosat-9", "IR3.9").unwrap();
let rad = conv.tb2radiance(300.0);
let tb = conv.radiance2tb(rad);
```

### RSR-based converter struct

```rust
use rustyspectral::radiance_tb::RadTbConverter;
use ndarray::arr1;

let conv = RadTbConverter::new(
    "EOS-Aqua", "modis", "20",
    wavelength.clone(), response.clone(),
);

let rad = conv.tb2radiance(&arr1(&[300.0]), true);  // normalized
let tb = conv.radiance2tb(&rad);

// Generate LUT
let (lut_tb, lut_rad) = conv.make_tb2rad_lut(0.1, true);
```

## Near-Infrared Reflectance (3.9µm)

```rust
use rustyspectral::reflectance::ReflectanceCalculator;
use ndarray::{arr1, Array1};

// Create with RSR
let wvl = Array1::linspace(3.6e-6, 3.95e-6, 36);
let resp = Array1::ones(36);

let calc = ReflectanceCalculator::new("EOS-Aqua", "modis")
    .with_rsr(wvl, resp)
    .with_solar_flux(2.002928);

// Single pixel
let sunz = arr1(&[80.0]);
let tb_nir = arr1(&[290.0]);    // 3.9µm brightness temperature
let tb_thermal = arr1(&[282.0]); // 11µm brightness temperature

let refl = calc.reflectance_from_tbs(&sunz, &tb_nir, &tb_thermal, None);
println!("Reflectance: {:.6}", refl[0]);  // ~0.251

// With CO₂ absorption correction (SEVIRI Rosenfeld method)
let tb_co2 = arr1(&[270.0]);  // 13.4µm CO₂ channel
let refl_co2 = calc.reflectance_from_tbs(
    &sunz, &tb_nir, &tb_thermal, Some(&tb_co2)
);

// Multiple pixels
let sunz = arr1(&[30.0, 50.0, 70.0]);
let tb_nir = arr1(&[295.0, 290.0, 285.0]);
let tb_thermal = arr1(&[283.0, 280.0, 277.0]);
let refl = calc.reflectance_from_tbs(&sunz, &tb_nir, &tb_thermal, None);
// refl.len() == 3, values in [0, 1]

// Solar radiance
let sr = calc.solar_radiance(&sunz);
// sr = solar_flux * cos(sunz_rad) / π

// CO₂ correction factors
let corr = calc.derive_rad39_corr(&tb_thermal, &arr1(&[270.0]));

// Emissive part (requires reflectance to be computed first)
let rad3x_t11 = calc.tb2radiance(&tb_thermal);
let rad3x = calc.tb2radiance(&tb_nir);
let emissive = calc.emissive_part_3x(&rad3x_t11, &refl, &rad3x, false);
// emissive = rad3x_t11 * (1 - r3x), NaN→rad3x fallback

// Solar radiance scaling
let sr = calc.solar_radiance(&sunz);

// With custom sunz_threshold and masking_limit
let calc = ReflectanceCalculator::new("NOAA-19", "avhrr3")
    .with_central_wavelength(3.78e-6)
    .with_solar_flux(2.0)
    .with_sunz_threshold(88.0)
    .with_masking_limit(None);  // no masking
```

## Rayleigh Atmospheric Correction

```rust
use rustyspectral::rayleigh::Rayleigh;
use ndarray::arr1;

// Initialize
let rcor = Rayleigh::new(
    "GOES-16",
    "abi",
    Some("us_standard"),              // atmosphere type
    Some("marine_clean_aerosol"),     // aerosol distribution
);

// By band name (requires RSR data)
let refl = rcor.get_reflectance(
    &arr1(&[50.0]),  // sun zenith (°)
    &arr1(&[40.0]),  // sat zenith (°)
    &arr1(&[160.0]), // azimuth difference (°)
    "ch2",           // band name
    None,            // no cloud relaxation
);

// By wavelength (micrometers, bypasses RSR lookup)
let refl = rcor.get_reflectance(
    &arr1(&[50.0, 60.0]),
    &arr1(&[40.0, 50.0]),
    &arr1(&[160.0, 160.0]),
    "0.64",          // 640 nm
    None,
);

// With cloud relaxation (red-band reflectance)
let redband = arr1(&[0.10, 0.15]);
let refl = rcor.get_reflectance(
    &arr1(&[50.0, 60.0]),
    &arr1(&[40.0, 50.0]),
    &arr1(&[160.0, 160.0]),
    "0.64",
    Some(&redband),
);

// Standalone functions
use rustyspectral::rayleigh::*;
let clipped = clip_angles_inside_coordinate_range(&arr1(&[79.0, 69.0, 32.0, f64::NAN]), 2.75);

let reduced = reduce_rayleigh_highzenith(
    &arr1(&[70.0, 65.0, 60.0]),
    &arr1(&[50.0, 50.0, 50.0]),
    70.0, 90.0, 1.0,
);
```

### Available atmosphere types

| Name (Python-style) | Name (Rust internal) |
|---------------------|---------------------|
| `subarctic summer` | `subarctic_summer` |
| `subarctic winter` | `subarctic_winter` |
| `midlatitude summer` | `midlatitude_summer` |
| `midlatitude winter` | `midlatitude_winter` |
| `tropical` | `tropical` |
| `us-standard` | `us_standard` |

Both naming conventions (hyphen/space and underscore) are accepted by
the Rust API.

### Available aerosol types

`antarctic_aerosol`, `continental_average_aerosol`, `continental_clean_aerosol`,
`continental_polluted_aerosol`, `desert_aerosol`, `marine_clean_aerosol`,
`marine_polluted_aerosol`, `marine_tropical_aerosol`, `rayleigh_only`,
`rural_aerosol`, `urban_aerosol`

## IR Atmospheric Correction

```rust
use rustyspectral::atm_correction_ir::*;
use ndarray::Array2;

// Using the struct
let atm = AtmosphericalCorrection::new("Suomi-NPP", "viirs");
let corrected = atm.get_correction(&sat_zenith, "M4", &brightness_temp);

// Or standalone function
let corrected = viewzen_corr(&brightness_temp, &sat_zenith);
```

The correction is based on the DWD parametric method:
- z=0°: `ΔT = tau0(T)`
- 0<z<90°: `ΔT = tau(T) · delta(z)`
- z≥90°: no correction

## Band Name Dictionaries

```rust
use rustyspectral::bandnames::get_bandnames;

let names = get_bandnames();

// Get SEVIRI band names
let seviri = names.get("seviri").unwrap();
assert_eq!(seviri.get("VIS006"), Some(&"VIS0.6"));
assert_eq!(seviri.get("IR_108"), Some(&"IR10.8"));

// Generic numeric channels
let generic = names.get("generic").unwrap();
assert_eq!(generic.get("20"), Some(&"ch20"));

// All sensors: generic, modis, seviri, viirs, avhrr3, abi,
//              agri, ahi, ami, fci, slstr, vii
```

## Utility Functions

```rust
use rustyspectral::utils::*;
use ndarray::arr1;

// Trapezoidal integration
let x = arr1(&[0.0, 1.0, 2.0, 3.0]);
let y = arr1(&[0.0, 1.0, 4.0, 9.0]);
let integral = trapezoid(&y, &x);  // 9.5

// Wavelength→wavenumber conversion
let (wn, wresp) = convert2wavenumber_rsr(&wavelength, &response);

// Weighted central wavelength
let wvl = arr1(&[0.5, 0.6, 0.7_f64]);
let resp = arr1(&[0.0, 1.0, 0.0_f64]);
let weight = arr1(&[1.0, 1.0, 1.0_f64]);
let cw = get_central_wave(&wvl, &resp, &weight);  // ~0.6

// Rayleigh-weighting example: weight = 1/λ⁴
let rayleigh_weight = wvl.mapv(|w| 1.0 / w.powi(4));
let cw_rayleigh = get_central_wave(&wvl, &resp, &rayleigh_weight);

// Sort and deduplicate
let (xs, ys) = sort_data(&x_vals, &y_vals);

// FWHM bandwidth
let fwhm = get_fullwidth_halfmax(&response, &wavelength);

// Integrated energy bounds
let (low, high) = get_bounds_integrated_energy(&response, &wavelength, 1.0);

// Wavelength range above threshold
let (min_wl, cwl, max_wl) = get_wave_range(&wavelength, &response, 0.15);

// Find band name from wavelength
let band = get_bandname_from_wavelength("modis", 0.67, &rsr, 0.1, false);

// Instrument name tools
assert!(are_instruments_identical("avhrr/1", "avhrr-1"));
let norm = check_and_adjust_instrument_name("GOES-16", "abi");  // "abi"

// Access global constants
println!("{:?}", AEROSOL_TYPES);
println!("{:?}", ATMOSPHERES);
```
