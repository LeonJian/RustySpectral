use approx::assert_relative_eq;
use ndarray::{arr1, Array1};
use rustyspectral::reflectance::*;

#[test]
fn test_get_as_array() {
    let result = get_as_array(42.0, 5);
    assert_eq!(result.len(), 5);
    for v in result.iter() {
        assert_relative_eq!(*v, 42.0);
    }
}

#[test]
fn test_terminator_limit() {
    assert_relative_eq!(TERMINATOR_LIMIT, 85.0);
}

#[test]
fn test_reflectance_from_tbs() {
    let wvl = Array1::linspace(3.6e-6, 3.95e-6, 36);
    let resp = Array1::ones(36);
    let calc = ReflectanceCalculator::new("NOAA-19", "avhrr3")
        .with_rsr(wvl, resp)
        .with_solar_flux(2.0);
    let sunz = arr1(&[30.0]);
    let tb_nir = arr1(&[290.0]);
    let tb_thermal = arr1(&[280.0]);
    let refl = calc.reflectance_from_tbs(&sunz, &tb_nir, &tb_thermal, None);
    assert!(refl[0] >= 0.0);
    assert!(refl[0] <= 1.0);
}

#[test]
fn test_reflectance_with_co2() {
    let wvl = Array1::linspace(3.6e-6, 3.95e-6, 36);
    let resp = Array1::ones(36);
    let calc = ReflectanceCalculator::new("NOAA-19", "avhrr3")
        .with_rsr(wvl, resp)
        .with_solar_flux(2.0);
    let sunz = arr1(&[30.0]);
    let tb_nir = arr1(&[290.0]);
    let tb_thermal = arr1(&[280.0]);
    let tb_co2 = arr1(&[270.0]);
    let refl = calc.reflectance_from_tbs(&sunz, &tb_nir, &tb_thermal, Some(&tb_co2));
    assert!(refl[0] >= 0.0);
    assert!(refl[0] <= 1.0);
}

#[test]
fn test_reflectance_terminator() {
    let wvl = Array1::linspace(3.6e-6, 3.95e-6, 36);
    let resp = Array1::ones(36);
    let calc = ReflectanceCalculator::new("NOAA-19", "avhrr3")
        .with_rsr(wvl, resp)
        .with_solar_flux(2.0);
    let sunz = arr1(&[TERMINATOR_LIMIT]);
    let tb_nir = arr1(&[290.0]);
    let tb_thermal = arr1(&[280.0]);
    let refl = calc.reflectance_from_tbs(&sunz, &tb_nir, &tb_thermal, None);
    assert!(refl[0].is_nan());
}

#[test]
fn test_reflectance_beyond_terminator() {
    let wvl = Array1::linspace(3.6e-6, 3.95e-6, 36);
    let resp = Array1::ones(36);
    let calc = ReflectanceCalculator::new("NOAA-19", "avhrr3")
        .with_rsr(wvl, resp)
        .with_solar_flux(2.0);
    let sunz = arr1(&[95.0]);
    let tb_nir = arr1(&[290.0]);
    let tb_thermal = arr1(&[280.0]);
    let refl = calc.reflectance_from_tbs(&sunz, &tb_nir, &tb_thermal, None);
    assert!(refl[0].is_nan());
}

#[test]
fn test_reflectance_default_limits() {
    let wvl = Array1::linspace(3.6e-6, 3.95e-6, 36);
    let resp = Array1::ones(36);
    let calc = ReflectanceCalculator::new("NOAA-19", "avhrr3")
        .with_rsr(wvl, resp)
        .with_solar_flux(2.0);
    let sunz = arr1(&[30.0]);
    let tb_nir = arr1(&[290.0]);
    let tb_thermal = arr1(&[260.0]);
    let refl = calc.reflectance_from_tbs(&sunz, &tb_nir, &tb_thermal, None);
    assert!(refl[0] >= 0.0);
    assert!(refl[0] <= 1.0);
}

#[test]
fn test_emissive_part() {
    let wvl = Array1::linspace(3.6e-6, 3.95e-6, 36);
    let resp = Array1::ones(36);
    let calc = ReflectanceCalculator::new("NOAA-19", "avhrr3").with_rsr(wvl, resp);
    let tb_nir = arr1(&[290.0]);
    let emiss_rad = calc.emissive_part(&tb_nir);
    assert!(emiss_rad[0] > 0.0);
}

#[test]
fn test_solar_radiance() {
    let calc = ReflectanceCalculator::new("test", "test").with_solar_flux(1400.0);
    let sunz = arr1(&[0.0, 60.0]);
    let sr = calc.solar_radiance(&sunz);
    assert_relative_eq!(sr[0], 1400.0 / std::f64::consts::PI, epsilon = 1e-10);
    assert_relative_eq!(sr[1], 1400.0 * 0.5 / std::f64::consts::PI, epsilon = 1e-10);
}

#[test]
fn test_multiple_pixels() {
    let wvl = Array1::linspace(3.6e-6, 3.95e-6, 36);
    let resp = Array1::ones(36);
    let calc = ReflectanceCalculator::new("test", "test")
        .with_rsr(wvl, resp)
        .with_solar_flux(2.0);
    let sunz = arr1(&[10.0, 30.0, 60.0]);
    let tb_nir = arr1(&[300.0, 295.0, 290.0]);
    let tb_thermal = arr1(&[285.0, 283.0, 280.0]);
    let refl = calc.reflectance_from_tbs(&sunz, &tb_nir, &tb_thermal, None);
    assert_eq!(refl.len(), 3);
    for &r in refl.iter() {
        assert!(r >= 0.0 && r <= 1.0);
    }
}
