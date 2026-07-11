use approx::{assert_abs_diff_eq, assert_relative_eq};
use ndarray::array;

use rustyspectral::radiance_tb::*;

const TEST_WVL: &[f64] = &[
    3.6123999, 3.6163599, 3.6264927, 3.6363862, 3.646468, 3.6564937, 3.6664478, 3.6765388,
    3.6865413, 3.6964585, 3.7065142, 3.716509, 3.7264658, 3.7364102, 3.7463682, 3.7563652,
    3.7664226, 3.7763396, 3.7863384, 3.7964207, 3.8063589, 3.8163606, 3.8264089, 3.8364836,
    3.8463381, 3.8563975, 3.8664163, 3.8763755, 3.8864797, 3.8964978, 3.9064275, 3.9164873,
    3.9264729, 3.9364026, 3.9465107, 3.9535347,
];

const TEST_RESP: &[f64] = &[
    0.01, 0.0118, 0.01987, 0.03226, 0.05028, 0.0849, 0.16645, 0.33792, 0.59106, 0.81815, 0.96077,
    0.92855, 0.86008, 0.8661, 0.87697, 0.85412, 0.88922, 0.9541, 0.95687, 0.91037, 0.91058,
    0.94256, 0.94719, 0.94808, 1.0, 0.92676, 0.67429, 0.44715, 0.27762, 0.14852, 0.07141, 0.04151,
    0.02925, 0.02085, 0.01414, 0.01,
];

fn test_wavelength() -> ndarray::Array1<f64> {
    ndarray::Array1::from_vec(TEST_WVL.iter().map(|v| v * 1e-6).collect())
}

fn test_response() -> ndarray::Array1<f64> {
    ndarray::Array1::from_vec(TEST_RESP.to_vec())
}

#[test]
fn test_radiance2tb_simple() {
    let wvl = test_wavelength();
    let resp = test_response();
    let central_wvl = 3.780_281_935e-6;
    let rad = tb2radiance_normalized(300.0, &wvl, &resp);
    let tb = radiance2tb(rad, central_wvl);
    assert_abs_diff_eq!(tb, 300.0, epsilon = 0.5);
}

#[test]
fn test_tb2radiance_single() {
    let wvl = test_wavelength();
    let resp = test_response();
    let rad = tb2radiance_simple(300.0, &wvl, &resp);
    assert!(rad > 0.0);
    // Higher Tb => higher radiance
    let rad2 = tb2radiance_simple(350.0, &wvl, &resp);
    assert!(rad2 > rad);
}

#[test]
fn test_tb2radiance_array() {
    let wvl = test_wavelength();
    let resp = test_response();
    let tbs = array![200.0, 270.0, 300.0, 302.0, 350.0_f64];
    let rads = tb2radiance_array(&tbs, &wvl, &resp);
    assert_eq!(rads.len(), 5);
    assert!(rads[0] > 0.0);
    assert!(rads[0] < rads[1]);
    assert!(rads[1] < rads[2]);
    assert!(rads[3] < rads[4]);
}

#[test]
fn test_tb2radiance_normalized() {
    let wvl = test_wavelength();
    let resp = test_response();
    let rad = tb2radiance_normalized(300.0, &wvl, &resp);
    assert!(rad > 0.0);
}

#[test]
fn test_make_tb2rad_lut() {
    let wvl = test_wavelength();
    let resp = test_response();
    let (lut_tb, lut_rad) = make_tb2rad_lut(&wvl, &resp, 0.1);
    assert!(lut_tb.len() > 0);
    assert_eq!(lut_tb.len(), lut_rad.len());
    // Monotonically increasing
    for i in 1..lut_rad.len() {
        assert!(lut_rad[i] > lut_rad[i - 1]);
    }
    // Range check
    assert!(lut_tb[0] >= 150.0);
    assert!(lut_tb[lut_tb.len() - 1] <= 360.0);
}

#[test]
fn test_seviri_tb2radiance() {
    // SEVIRI IR3.9 Meteosat-9 parameters
    let vc = 2568.832 * 100.0; // wavenumber in SI
    let alpha = 0.9954;
    let beta = 3.438;
    let rad = seviri_tb2radiance(300.0, vc, alpha, beta);
    assert!(rad > 0.0);
    // SEVIRI regression should give ~9.8e-6
    assert_relative_eq!(rad, 9.797_091e-6, epsilon = 1e-8);
}

#[test]
fn test_seviri_radiance2tb() {
    let vc = 2568.832 * 100.0;
    let alpha = 0.9954;
    let beta = 3.438;
    let rad = seviri_tb2radiance(300.0, vc, alpha, beta);
    let tb = seviri_radiance2tb(rad, vc, alpha, beta);
    assert_relative_eq!(tb, 300.0, epsilon = 1e-4);
}

#[test]
fn test_seviri_multiple_temperatures() {
    let vc = 2568.832 * 100.0;
    let alpha = 0.9954;
    let beta = 3.438;
    for tb in [200.0, 250.0, 300.0, 350.0].iter() {
        let rad = seviri_tb2radiance(*tb, vc, alpha, beta);
        let tb_back = seviri_radiance2tb(rad, vc, alpha, beta);
        assert_relative_eq!(tb_back, *tb, epsilon = 1e-4);
    }
}
