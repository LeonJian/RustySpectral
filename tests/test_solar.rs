use approx::{assert_abs_diff_eq, assert_relative_eq};
use ndarray::array;

use rustyspectral::rsr::Rsr;
use rustyspectral::solar::SolarIrradianceSpectrum;

const SOLAR_SPECTRUM_FILE: &str = "data/e490_00a.dat";

fn test_rsr() -> Rsr {
    let wvl = array![
        3.6123999,
        3.6163599,
        3.6264927,
        3.6363862,
        3.646468,
        3.6564937,
        3.6664478,
        3.6765388,
        3.6865413,
        3.6964585,
        3.7065142,
        3.716509,
        3.7264658,
        3.7364102,
        3.7463682,
        3.7563652,
        3.7664226,
        3.7763396,
        3.7863384,
        3.7964207,
        3.8063589,
        3.8163606,
        3.8264089,
        3.8364836,
        3.8463381,
        3.8563975,
        3.8664163,
        3.8763755,
        3.8864797,
        3.8964978,
        3.9064275,
        3.9164873,
        3.9264729,
        3.9364026,
        3.9465107,
        3.9535347_f64,
    ];
    let resp = array![
        0.01, 0.0118, 0.01987, 0.03226, 0.05028, 0.0849, 0.16645, 0.33792, 0.59106, 0.81815,
        0.96077, 0.92855, 0.86008, 0.8661, 0.87697, 0.85412, 0.88922, 0.9541, 0.95687, 0.91037,
        0.91058, 0.94256, 0.94719, 0.94808, 1.0, 0.92676, 0.67429, 0.44715, 0.27762, 0.14852,
        0.07141, 0.04151, 0.02925, 0.02085, 0.01414, 0.01_f64,
    ];
    Rsr::new(wvl, resp)
}

#[test]
fn test_load_solar_spectrum() {
    let solar = SolarIrradianceSpectrum::new(SOLAR_SPECTRUM_FILE, 0.005);
    assert_eq!(solar.wavelength.len(), 1697);
    assert_eq!(solar.irradiance.len(), 1697);
    assert!(solar.wavelength[0] < solar.wavelength[1]);
}

#[test]
fn test_solar_constant() {
    let solar = SolarIrradianceSpectrum::new(SOLAR_SPECTRUM_FILE, 0.005);
    let sc = solar.solar_constant();
    // Solar constant ~1361 W/m^2 (range 1360-1362)
    println!("SOLAR CONSTANT: {}", sc);
    assert!(sc > 1364.0 && sc < 1368.0, "sc = {}", sc);
}

#[test]
fn test_interpolate() {
    let mut solar = SolarIrradianceSpectrum::new(SOLAR_SPECTRUM_FILE, 0.005);
    solar.interpolate(0.001, Some((0.200, 0.240)));
    let ipol = solar.ipol_wavelength.as_ref().unwrap();
    assert!(!ipol.is_empty());
    // Start and end matching
    assert_relative_eq!(ipol[0], 0.200, epsilon = 1e-6);
    assert_relative_eq!(ipol[ipol.len() - 1], 0.240, epsilon = 1e-6);
    // Check the interpolation grid
    let expected_n = f64::round((0.240 - 0.200) / 0.001) as usize + 1;
    assert_eq!(ipol.len(), expected_n);
}

#[test]
fn test_inband_solarflux() {
    let mut solar = SolarIrradianceSpectrum::new(SOLAR_SPECTRUM_FILE, 0.005);
    let rsr = test_rsr();
    let flux = solar.inband_solarflux(&rsr, 1.0);
    // Python reference: 2.002927627
    assert_abs_diff_eq!(flux, 2.002_927_627, epsilon = 1e-3);
}

#[test]
fn test_inband_solarflux_wavenumber() {
    let mut solar = SolarIrradianceSpectrum::new(SOLAR_SPECTRUM_FILE, 0.005);
    solar.set_wavespace_wavenumber();
    let rsr = test_rsr();
    let flux = solar.inband_solarflux(&rsr, 1.0);
    assert!(flux > 0.0);
}

#[test]
fn test_inband_solarirradiance() {
    let mut solar = SolarIrradianceSpectrum::new(SOLAR_SPECTRUM_FILE, 0.005);
    let rsr = test_rsr();
    let irradiance = solar.inband_solarirradiance(&rsr, 1.0);
    assert!(irradiance > 0.0);
    // Irradiance should differ from flux (flux is integrated, irradiance is normalized)
    let flux = solar.inband_solarflux(&rsr, 1.0);
    assert!(irradiance != flux);
}

#[test]
fn test_solar_constant_wavenumber() {
    let mut solar = SolarIrradianceSpectrum::new(SOLAR_SPECTRUM_FILE, 0.005);
    solar.set_wavespace_wavenumber();
    let sc = solar.solar_constant();
    assert!(sc > 0.0);
}
