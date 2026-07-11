use approx::assert_relative_eq;
use ndarray::array;

use rustyspectral::rsr::*;

#[test]
fn test_rsr_single_detector() {
    let wvl = array![3.6, 3.7, 3.8, 3.9, 4.0_f64];
    let resp = array![0.0, 0.5, 1.0, 0.5, 0.0_f64];
    let rsr = Rsr::new(wvl.clone(), resp.clone());
    assert_eq!(rsr.wavelength.len(), 5);
    assert_eq!(rsr.response.len(), 5);
    assert_relative_eq!(rsr.wavelength[2], 3.8);
    assert_relative_eq!(rsr.response[2], 1.0);
    // Central wavelength should be ~3.8 (peak of symmetric triangular-like response)
    assert_relative_eq!(rsr.central_wavelength, 3.8, epsilon = 1e-10);
}

#[test]
fn test_rsr_to_wavenumber() {
    let wvl = array![3.6, 3.7, 3.8, 3.9, 4.0_f64];
    let resp = array![0.0, 0.5, 1.0, 0.5, 0.0_f64];
    let rsr = Rsr::new(wvl, resp);
    let (wnum, wresp) = rsr.to_wavenumber();
    assert_eq!(wnum.len(), 5);
    assert_eq!(wresp.len(), 5);
    // Wavenumbers should be descending (1/wavelength relationship, flipped)
    assert!(wnum[0] < wnum[4]);
    // Responses flipped
    assert_relative_eq!(wresp[0], 0.0, epsilon = 1e-10);
    assert_relative_eq!(wresp[2], 1.0, epsilon = 1e-10);
}

#[test]
fn test_per_detector_rsr() {
    let wvl = array![3.6, 3.7, 3.8, 3.9, 4.0_f64];
    let resp = array![0.0, 0.5, 1.0, 0.5, 0.0_f64];
    let det = PerDetectorRsr::new("det-1".to_string(), wvl, resp);
    assert_eq!(det.name, "det-1");
    assert_relative_eq!(det.central_wavelength, 3.8, epsilon = 1e-10);
}

#[test]
fn test_detector_rsr_multi() {
    let wvl1 = array![3.6, 3.7, 3.8_f64];
    let resp1 = array![0.0, 1.0, 0.0_f64];
    let det1 = PerDetectorRsr::new("det-1".to_string(), wvl1, resp1);

    let wvl2 = array![4.6, 4.7, 4.8_f64];
    let resp2 = array![0.0, 1.0, 0.0_f64];
    let det2 = PerDetectorRsr::new("det-2".to_string(), wvl2, resp2);

    let drsr = DetectorRsr {
        detectors: vec![det1, det2],
    };
    assert_eq!(drsr.detectors.len(), 2);
    assert_eq!(drsr.detectors[0].name, "det-1");
    assert_eq!(drsr.detectors[1].name, "det-2");
    assert_relative_eq!(drsr.detectors[0].central_wavelength, 3.7, epsilon = 1e-10);
    assert_relative_eq!(drsr.detectors[1].central_wavelength, 4.7, epsilon = 1e-10);
}
