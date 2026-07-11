use approx::assert_abs_diff_eq;
use ndarray::array;
use rustyspectral::rayleigh::*;

#[test]
fn test_clip_angles_with_nans() {
    let zenith = array![79.0, 69.0, 32.0, f64::NAN];
    let result = clip_angles_inside_coordinate_range(&zenith, 2.75);
    let expected = array![68.676_313_74, 68.676_313_74, 32.0, 0.0_f64];
    for i in 0..4 {
        assert_abs_diff_eq!(result[i], expected[i], epsilon = 1e-4);
    }
}

#[test]
fn test_reduce_rayleigh_no_reduction() {
    let sun_zenith = array![70.0, 65.0, 60.0_f64];
    let in_rayleigh = array![50.0, 50.0, 50.0_f64];
    let result = reduce_rayleigh_highzenith(&sun_zenith, &in_rayleigh, 70.0, 90.0, 1.0);
    for i in 0..3 {
        assert_abs_diff_eq!(result[i], in_rayleigh[i], epsilon = 1e-6);
    }
}

#[test]
fn test_reduce_rayleigh_moderate() {
    let sun_zenith = array![70.0, 65.0, 60.0_f64];
    let in_rayleigh = array![50.0, 50.0, 50.0_f64];
    let result = reduce_rayleigh_highzenith(&sun_zenith, &in_rayleigh, 30.0, 90.0, 1.0);
    let expected = array![16.666_666_67, 20.833_333_33, 25.0_f64];
    for i in 0..3 {
        assert_abs_diff_eq!(result[i], expected[i], epsilon = 1e-3);
    }
}

#[test]
fn test_reduce_rayleigh_extreme() {
    let sun_zenith = array![70.0, 65.0, 60.0_f64];
    let in_rayleigh = array![50.0, 50.0, 50.0_f64];
    let result = reduce_rayleigh_highzenith(&sun_zenith, &in_rayleigh, 30.0, 90.0, 1.5);
    let expected = array![0.0, 6.25, 12.5_f64];
    for i in 0..3 {
        assert_abs_diff_eq!(result[i], expected[i], epsilon = 1e-3);
    }
}

#[test]
fn test_wavelength_index_and_factor() {
    let wvl_coord = array![631.0_f64, 636.0_f64];
    let (idx, factor) = get_wavelength_index_and_factor(&wvl_coord, 634.0);
    assert_eq!(idx, 1);
    assert_abs_diff_eq!(factor, (636.0 - 634.0) / (636.0 - 631.0), epsilon = 1e-10);
}

#[test]
fn test_normalize_sensor_known() {
    assert_eq!(normalize_sensor("GOES-16", "abi"), "abi");
    assert_eq!(normalize_sensor("NOAA-19", "avhrr/3"), "avhrr3");
    assert_eq!(normalize_sensor("FY-4A", "agri"), "agri");
    assert_eq!(normalize_sensor("Himawari-8", "ahi"), "ahi");
    assert_eq!(normalize_sensor("NOAA-20", "viirs"), "viirs");
    assert_eq!(normalize_sensor("Meteosat-9", "seviri"), "seviri");
}

#[test]
fn test_normalize_sensor_unknown_platform_uses_sensor() {
    assert_eq!(normalize_sensor("Unknown", "myinstr"), "myinstr");
}

#[test]
fn test_aerosol_types_list() {
    let types = rustyspectral::utils::AEROSOL_TYPES;
    assert_eq!(types.len(), 11);
    assert!(types.contains(&"marine_clean_aerosol"));
    assert!(types.contains(&"desert_aerosol"));
    assert!(types.contains(&"rayleigh_only"));
}

#[test]
fn test_atmospheres_list() {
    let atms = rustyspectral::utils::ATMOSPHERES;
    assert_eq!(atms.len(), 6);
    let names: Vec<&str> = atms.iter().map(|(n, _)| *n).collect();
    assert!(names.contains(&"midlatitude_summer"));
    assert!(names.contains(&"us_standard"));
}
