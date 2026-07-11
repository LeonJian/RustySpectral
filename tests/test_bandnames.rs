use rustyspectral::bandnames::get_bandnames;

#[test]
fn test_generic_bandnames() {
    let names = get_bandnames();
    let generic = names.get("generic").unwrap();
    assert_eq!(generic.get("VIS006"), Some(&"VIS0.6"));
    assert_eq!(generic.get("IR_039"), Some(&"IR3.9"));
    assert_eq!(generic.get("I01"), Some(&"I1"));
    assert_eq!(generic.get("M04"), Some(&"M4"));
    assert_eq!(generic.get("C01"), Some(&"ch1"));
    assert_eq!(generic.get("C16"), Some(&"ch16"));
}

#[test]
fn test_generic_numeric_channels() {
    let names = get_bandnames();
    let generic = names.get("generic").unwrap();
    assert_eq!(generic.get("1"), Some(&"ch1"));
    assert_eq!(generic.get("20"), Some(&"ch20"));
    assert_eq!(generic.get("36"), Some(&"ch36"));
}

#[test]
fn test_seviri_bandnames() {
    let names = get_bandnames();
    let seviri = names.get("seviri").unwrap();
    assert_eq!(seviri.get("VIS006"), Some(&"VIS0.6"));
    assert_eq!(seviri.get("IR_108"), Some(&"IR10.8"));
    assert_eq!(seviri.get("HRV"), Some(&"HRV"));
}

#[test]
fn test_viirs_bandnames() {
    let names = get_bandnames();
    let viirs = names.get("viirs").unwrap();
    assert_eq!(viirs.get("I01"), Some(&"I1"));
    assert_eq!(viirs.get("M04"), Some(&"M4"));
}

#[test]
fn test_abi_bandnames() {
    let names = get_bandnames();
    let abi = names.get("abi").unwrap();
    assert_eq!(abi.get("C01"), Some(&"ch1"));
    assert_eq!(abi.get("C16"), Some(&"ch16"));
}

#[test]
fn test_ahi_bandnames() {
    let names = get_bandnames();
    let ahi = names.get("ahi").unwrap();
    assert_eq!(ahi.get("B01"), Some(&"ch1"));
    assert_eq!(ahi.get("B16"), Some(&"ch16"));
}

#[test]
fn test_avhrr3_bandnames() {
    let names = get_bandnames();
    let avhrr3 = names.get("avhrr3").unwrap();
    assert_eq!(avhrr3.get("1"), Some(&"ch1"));
    assert_eq!(avhrr3.get("4"), Some(&"ch4"));
    assert_eq!(avhrr3.get("3b"), Some(&"ch3b"));
}

#[test]
fn test_fci_bandnames() {
    let names = get_bandnames();
    let fci = names.get("fci").unwrap();
    assert_eq!(fci.get("vis_04"), Some(&"VIS0.4"));
    assert_eq!(fci.get("ir_38"), Some(&"IR3.8_HR"));
}

#[test]
fn test_modis_keeps_numeric() {
    let names = get_bandnames();
    let modis = names.get("modis").unwrap();
    assert_eq!(modis.get("1"), Some(&"1"));
    assert_eq!(modis.get("20"), Some(&"20"));
}

#[test]
fn test_unknown_sensor_defaults_to_generic() {
    let names = get_bandnames();
    assert!(names.contains_key("generic"));
    assert!(!names.contains_key("nonexistent"));
}
