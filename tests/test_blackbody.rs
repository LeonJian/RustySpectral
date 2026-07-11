use approx::assert_abs_diff_eq;
use ndarray::{arr2, array};

use rustyspectral::blackbody::*;

// Reference values from pyspectral Python tests
const RAD_11MICRON_300K: f64 = 9573176.935507433;
const RAD_11MICRON_301K: f64 = 9714686.576498277;
const WN_RAD_11MICRON_300K: f64 = 0.00115835441353;
const WN_RAD_11MICRON_301K: f64 = 0.00117547716523;

const H: f64 = 6.626_069_57e-34;
const K: f64 = 1.380_648_8e-23;
const C: f64 = 2.997_924_58e8;

#[test]
fn test_physical_constants() {
    assert_abs_diff_eq!(H_PLANCK, H, epsilon = 1e-40);
    assert_abs_diff_eq!(K_BOLTZMANN, K, epsilon = 1e-30);
    assert_abs_diff_eq!(C_SPEED, C, epsilon = 1e-1);
}

#[test]
fn test_planck_wavelength_single_300k() {
    let wavel = 11e-6;
    let rad = planck(wavel, 300.0);
    assert_abs_diff_eq!(rad, RAD_11MICRON_300K, epsilon = 1e-4);
}

#[test]
fn test_planck_wavelength_single_301k() {
    let wavel = 11e-6;
    let rad = planck(wavel, 301.0);
    assert_abs_diff_eq!(rad, RAD_11MICRON_301K, epsilon = 1e-4);
}

#[test]
fn test_planck_wn_single_300k() {
    let wavenumber = 90909.1;
    let rad = planck_wn(wavenumber, 300.0);
    assert_abs_diff_eq!(rad, WN_RAD_11MICRON_300K, epsilon = 3e-10);
}

#[test]
fn test_planck_wn_single_301k() {
    let wavenumber = 90909.1;
    let rad = planck_wn(wavenumber, 301.0);
    assert_abs_diff_eq!(rad, WN_RAD_11MICRON_301K, epsilon = 3e-10);
}

#[test]
fn test_rad2temp_single_300k() {
    let wavel = 11e-6;
    let rad = planck(wavel, 300.0);
    let t = blackbody_rad2temp(wavel, rad);
    assert_abs_diff_eq!(t, 300.0, epsilon = 1e-8);
}

#[test]
fn test_rad2temp_single_301k() {
    let wavel = 11e-6;
    let rad = planck(wavel, 301.0);
    let t = blackbody_rad2temp(wavel, rad);
    assert_abs_diff_eq!(t, 301.0, epsilon = 1e-8);
}

#[test]
fn test_rad2temp_wn_single_300k() {
    let wavenumber = 90909.1;
    let rad = planck_wn(wavenumber, 300.0);
    let t = blackbody_wn_rad2temp(wavenumber, rad);
    assert_abs_diff_eq!(t, 300.0, epsilon = 1e-8);
}

#[test]
fn test_rad2temp_wn_single_301k() {
    let wavenumber = 90909.1;
    let rad = planck_wn(wavenumber, 301.0);
    let t = blackbody_wn_rad2temp(wavenumber, rad);
    assert_abs_diff_eq!(t, 301.0, epsilon = 1e-8);
}

#[test]
fn test_planck_wavelength_array() {
    let wavel = 10e-6;
    let temps = arr2(&[[300.0, 301.0], [299.0, 298.0], [279.0, 286.0]]);
    let result = planck_array_wavelength(wavel, &temps);
    assert_eq!(result.shape(), &[3, 2]);
    // Sanity: 300K > 298K radiance
    assert!(result[[0, 0]] > result[[1, 1]]);
    // Sanity: 301K > 300K radiance
    assert!(result[[0, 1]] > result[[0, 0]]);
}

#[test]
fn test_rad2temp_array() {
    let wavenumber = 90909.1;
    let radiances = arr2(&[[0.001, 0.0009], [0.0012, 0.0018]]);
    let temps = blackbody_wn_rad2temp_array(wavenumber, &radiances);
    let expected = arr2(&[
        [290.327_691_6, 283.761_154_41],
        [302.418_133_0, 333.141_416_4],
    ]);
    assert_abs_diff_eq!(temps[[0, 0]], expected[[0, 0]], epsilon = 1e-5);
    assert_abs_diff_eq!(temps[[0, 1]], expected[[0, 1]], epsilon = 1e-5);
    assert_abs_diff_eq!(temps[[1, 0]], expected[[1, 0]], epsilon = 1e-5);
    assert_abs_diff_eq!(temps[[1, 1]], expected[[1, 1]], epsilon = 1e-5);
}

#[test]
fn test_zero_radiance_does_not_panic() {
    let wavel = 11e-6;
    let t = blackbody_rad2temp(wavel, 0.0);
    assert!(t.is_nan());
}

#[test]
fn test_zero_temperature_does_not_panic() {
    let wavel = 11e-6;
    let rad = planck(wavel, 0.0);
    assert!(rad.is_nan());
}

#[test]
fn test_planck_wn_zero_temperature_does_not_panic() {
    let wavenumber = 90909.1;
    let rad = planck_wn(wavenumber, 0.0);
    assert!(rad.is_nan());
}

#[test]
fn test_wn_rad2temp_zero_radiance_does_not_panic() {
    let wavenumber = 90909.1;
    let t = blackbody_wn_rad2temp(wavenumber, 0.0);
    assert!(t.is_nan());
}

#[test]
fn test_planck_wavelength_vector() {
    let wavels = array![10e-6, 11e-6, 12e-6];
    let rads = planck(wavels[0], 300.0);
    assert!(rads.is_finite());
}
