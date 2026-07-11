use approx::assert_relative_eq;
use rustyspectral::reflectance::*;

#[test]
fn test_get_as_array_scalar() {
    let result = get_as_array(42.0, None);
    assert_eq!(result.len(), 1);
    assert_relative_eq!(result[0], 42.0);
}

#[test]
fn test_get_as_array_with_shape() {
    let result = get_as_array(3.0, Some(&[5]));
    assert_eq!(result.len(), 5);
    for v in result.iter() {
        assert_relative_eq!(*v, 3.0);
    }
}

#[test]
fn test_terminator_limit() {
    assert_relative_eq!(TERMINATOR_LIMIT, 89.0);
}

#[test]
fn test_reflectance_from_tbs() {
    let calc = ReflectanceCalculator::new("NOAA-19", "avhrr3");
    // Near-nadir, moderate solar zenith, warm scene
    let refl = calc.reflectance_from_tbs(30.0, 290.0, 280.0, None);
    assert!(refl >= 0.0 && refl <= 1.5);
}

#[test]
fn test_reflectance_with_co2() {
    let calc = ReflectanceCalculator::new("NOAA-19", "avhrr3");
    let refl = calc.reflectance_from_tbs(30.0, 290.0, 280.0, Some(270.0));
    assert!(refl >= 0.0 && refl <= 1.5);
}

#[test]
fn test_reflectance_terminal() {
    let calc = ReflectanceCalculator::new("NOAA-19", "avhrr3");
    // At terminator, reflectance should be 0 (no solar illumination)
    let refl = calc.reflectance_from_tbs(TERMINATOR_LIMIT, 290.0, 280.0, None);
    assert!(refl == 0.0);
}

#[test]
fn test_reflectance_beyond_terminator() {
    let calc = ReflectanceCalculator::new("NOAA-19", "avhrr3");
    let refl = calc.reflectance_from_tbs(95.0, 290.0, 280.0, None);
    assert!(refl == 0.0);
}

#[test]
fn test_reflectance_default_limits() {
    let calc = ReflectanceCalculator::new("NOAA-19", "avhrr3");
    let refl = calc.reflectance_from_tbs(30.0, 290.0, 260.0, None);
    // Should be clamped between 0 and 1
    assert!(refl >= 0.0 && refl <= 1.0);
}

#[test]
fn test_emissive_part() {
    let calc = ReflectanceCalculator::new("NOAA-19", "avhrr3");
    let emiss_rad = calc.emissive_part(290.0);
    assert!(emiss_rad > 0.0);
}
