use approx::assert_abs_diff_eq;
use ndarray::array;
use rustyspectral::utils::*;

fn test_rsr_data() -> (ndarray::Array1<f64>, ndarray::Array1<f64>) {
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
    (wvl, resp)
}

fn ch3_rsr() -> (ndarray::Array1<f64>, ndarray::Array1<f64>) {
    let wvl = array![
        0.55518544,
        0.56779468,
        0.58099002,
        0.59481323,
        0.60931027,
        0.62453163,
        0.64053291,
        0.65737575,
        0.67512828,
        0.69386619,
        0.71367401_f64,
    ];
    let resp = array![
        2.61000005e-05,
        1.07899999e-04,
        3.26119992e-03,
        2.90650606e-01,
        9.02460396e-01,
        9.60878074e-01,
        9.97266889e-01,
        9.94823873e-01,
        7.18220174e-01,
        8.31819978e-03,
        9.34999989e-05_f64,
    ];
    (wvl, resp)
}

#[test]
fn test_trapezoid() {
    let x = array![0.0, 1.0, 2.0, 3.0];
    let y = array![0.0, 1.0, 4.0, 9.0];
    let result = trapezoid(&y, &x);
    assert_abs_diff_eq!(result, 9.5, epsilon = 1e-10);
}

#[test]
fn test_central_wavelength() {
    let (wvl, resp) = test_rsr_data();
    let cw = get_central_wave(&wvl, &resp, 1.0);
    // From Python: get_central_wave gives ~3.78028 (microns)
    assert_abs_diff_eq!(cw, 3.780_281_935, epsilon = 1e-6);
}

#[test]
fn test_convert2wavenumber() {
    let (wvl, resp) = test_rsr_data();
    let (wnum, wresp) = convert2wavenumber_rsr(&wvl, &resp);

    // Check first/last values match Python reference
    let expected_first = 1.0 / (1e-4 * wvl[35]);
    let expected_last = 1.0 / (1e-4 * wvl[0]);
    assert_abs_diff_eq!(wnum[0], expected_first, epsilon = 1e-3);
    assert_abs_diff_eq!(wnum[35], expected_last, epsilon = 1e-3);

    // Response should be flipped
    assert_abs_diff_eq!(wresp[0], resp[35], epsilon = 1e-10);
    assert_abs_diff_eq!(wresp[35], resp[0], epsilon = 1e-10);
}

#[test]
fn test_sort_data() {
    let x = array![1.0, 5.6, 30.0, 2.1, 108.2, 57.8, 1e9, 2.1_f64];
    let y = array![45.0, 92.0, 20.0, 10.0, 15.0, 67.0, 108.0, 15.0_f64];

    let (x_sorted, y_sorted) = sort_data(&x, &y);

    let expected_x = array![1.0, 2.1, 5.6, 30.0, 57.8, 108.2, 1e9_f64];
    let expected_y = array![45.0, 10.0, 92.0, 20.0, 67.0, 15.0, 108.0_f64];

    for i in 0..expected_x.len() {
        assert_abs_diff_eq!(x_sorted[i], expected_x[i], epsilon = 1e-10);
        assert_abs_diff_eq!(y_sorted[i], expected_y[i], epsilon = 1e-10);
    }
}

#[test]
fn test_sort_data_no_duplicates() {
    let x = array![3.0, 1.0, 2.0_f64];
    let y = array![30.0, 10.0, 20.0_f64];
    let (x_sorted, y_sorted) = sort_data(&x, &y);
    let expected_x = array![1.0, 2.0, 3.0_f64];
    let expected_y = array![10.0, 20.0, 30.0_f64];
    for i in 0..3 {
        assert_abs_diff_eq!(x_sorted[i], expected_x[i]);
        assert_abs_diff_eq!(y_sorted[i], expected_y[i]);
    }
}

#[test]
fn test_fwhm() {
    let (wvl, resp) = ch3_rsr();
    let fwhm = get_fullwidth_halfmax(&resp, &wvl);
    assert_abs_diff_eq!(fwhm, 0.065_818_01, epsilon = 1e-5);
}

#[test]
fn test_fwhm_with_known_rsr() {
    let (wvl, resp) = test_rsr_data();
    let fwhm = get_fullwidth_halfmax(&resp, &wvl);
    assert!(fwhm > 0.0);
}

#[test]
fn test_integrated_energy_one_percent() {
    let (wvl, resp) = ch3_rsr();
    let (low, high) = get_bounds_integrated_energy(&resp, &wvl, 1.0);
    assert_abs_diff_eq!(low, 0.594_813_23, epsilon = 1e-5);
    assert_abs_diff_eq!(high, 0.657_375_75, epsilon = 1e-5);
}

#[test]
fn test_integrated_energy_ten_percent() {
    let (wvl, resp) = ch3_rsr();
    let (low, high) = get_bounds_integrated_energy(&resp, &wvl, 10.0);
    assert_abs_diff_eq!(low, 0.609_310_27, epsilon = 1e-5);
    assert_abs_diff_eq!(high, 0.657_375_75, epsilon = 1e-5);
}

#[test]
fn test_get_wave_range() {
    let (wvl, resp) = ch3_rsr();
    let (min_wvl, cwl, max_wvl) = get_wave_range(&wvl, &resp, 0.15);
    assert_abs_diff_eq!(min_wvl, 0.59481323, epsilon = 1e-5);
    assert_abs_diff_eq!(max_wvl, 0.67512828, epsilon = 1e-5);
    assert!(cwl > min_wvl && cwl < max_wvl);
}

#[test]
fn test_get_wave_range_higher_threshold() {
    let (wvl, resp) = ch3_rsr();
    let (min_wvl, _cwl, max_wvl) = get_wave_range(&wvl, &resp, 0.5);
    assert_abs_diff_eq!(min_wvl, 0.60931027, epsilon = 1e-5);
    assert_abs_diff_eq!(max_wvl, 0.67512828, epsilon = 1e-5);
}

#[test]
fn test_are_instruments_identical_same() {
    assert!(are_instruments_identical("abi", "abi"));
}

#[test]
fn test_are_instruments_identical_different() {
    assert!(!are_instruments_identical("abi", "viirs"));
}

#[test]
fn test_are_instruments_identical_avhrr_slash() {
    assert!(are_instruments_identical("avhrr/1", "avhrr-1"));
}

#[test]
fn test_are_instruments_identical_avhrr_reverse() {
    assert!(are_instruments_identical("avhrr-3", "avhrr/3"));
}

#[test]
fn test_are_instruments_identical_not_match() {
    assert!(!are_instruments_identical("avhrr/1", "avhrr/3"));
}
