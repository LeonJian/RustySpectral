use ndarray::Array1;
use std::collections::HashMap;
use std::path::{Path, PathBuf};

pub fn trapezoid(y: &Array1<f64>, x: &Array1<f64>) -> f64 {
    let n = y.len();
    if n < 2 {
        return 0.0;
    }
    let mut sum = 0.0;
    for i in 1..n {
        sum += (y[i] + y[i - 1]) * (x[i] - x[i - 1]);
    }
    0.5 * sum
}

pub fn get_central_wave(wav: &Array1<f64>, resp: &Array1<f64>, weight: &Array1<f64>) -> f64 {
    let numerator = trapezoid(&(resp * wav * weight), wav);
    let denominator = trapezoid(&(resp * weight), wav);
    if denominator == 0.0 {
        return f64::NAN;
    }
    numerator / denominator
}

pub fn convert2wavenumber_rsr(
    wavelength: &Array1<f64>,
    response: &Array1<f64>,
) -> (Array1<f64>, Array1<f64>) {
    let n = wavelength.len();
    let mut wavenumber = Array1::zeros(n);
    let mut resp_out = Array1::zeros(n);

    for i in 0..n {
        let j = n - 1 - i;
        wavenumber[i] = 1.0 / (1e-4 * wavelength[j]);
        resp_out[i] = response[j];
    }

    (wavenumber, resp_out)
}

pub fn sort_data(x_vals: &Array1<f64>, y_vals: &Array1<f64>) -> (Array1<f64>, Array1<f64>) {
    let n = x_vals.len();

    let mut indices: Vec<usize> = (0..n).collect();
    indices.sort_by(|&a, &b| {
        x_vals[a]
            .partial_cmp(&x_vals[b])
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut sorted_x = Array1::zeros(n);
    let mut sorted_y = Array1::zeros(n);
    for (i, &idx) in indices.iter().enumerate() {
        sorted_x[i] = x_vals[idx];
        sorted_y[i] = y_vals[idx];
    }

    let mut keep = vec![true; n];
    for i in 1..n {
        if sorted_x[i] <= sorted_x[i - 1] {
            keep[i] = false;
        }
    }

    let deduped: Vec<_> = (0..n).filter(|&i| keep[i]).collect();
    let m = deduped.len();
    let mut result_x = Array1::zeros(m);
    let mut result_y = Array1::zeros(m);
    for (i, &idx) in deduped.iter().enumerate() {
        result_x[i] = sorted_x[idx];
        result_y[i] = sorted_y[idx];
    }

    (result_x, result_y)
}

pub fn get_fullwidth_halfmax(rsp: &Array1<f64>, wvl: &Array1<f64>) -> f64 {
    let half_max = 0.5;
    let indices: Vec<usize> = rsp
        .iter()
        .enumerate()
        .filter(|(_, &v)| v >= half_max)
        .map(|(i, _)| i)
        .collect();

    if indices.len() < 2 {
        return f64::NAN;
    }

    wvl[indices[indices.len() - 1]] - wvl[indices[0]]
}

pub fn get_bounds_integrated_energy(
    rsp: &Array1<f64>,
    wvl: &Array1<f64>,
    ener_perc_lim: f64,
) -> (f64, f64) {
    let n = rsp.len();
    let mut crs = Array1::zeros(n);
    crs[0] = rsp[0];
    for i in 1..n {
        crs[i] = crs[i - 1] + rsp[i];
    }
    let max_val = crs[n - 1];
    for i in 0..n {
        crs[i] = crs[i] / max_val * 100.0;
    }

    let low: Vec<usize> = crs
        .iter()
        .enumerate()
        .filter(|(_, &v)| v >= ener_perc_lim)
        .map(|(i, _)| i)
        .collect();

    let high: Vec<usize> = crs
        .iter()
        .enumerate()
        .filter(|(_, &v)| v <= (100.0 - ener_perc_lim))
        .map(|(i, _)| i)
        .collect();

    let min_wvl = wvl[low[0]];
    let max_wvl = wvl[high[high.len() - 1]];

    (min_wvl, max_wvl)
}

pub fn get_wave_range(wvl: &Array1<f64>, resp: &Array1<f64>, threshold: f64) -> (f64, f64, f64) {
    let cwl = get_central_wave(wvl, resp, &Array1::from_elem(wvl.len(), 1.0));

    let pts: Vec<usize> = resp
        .iter()
        .enumerate()
        .filter(|(_, &v)| v > threshold)
        .map(|(i, _)| i)
        .collect();

    let min_wvl = wvl[pts[0]];
    let max_wvl = wvl[pts[pts.len() - 1]];

    (min_wvl, cwl, max_wvl)
}

pub fn are_instruments_identical(name1: &str, name2: &str) -> bool {
    if name1 == name2 {
        return true;
    }
    let translate = |s: &str| -> String {
        match s {
            "avhrr-1" => "avhrr/1".to_string(),
            "avhrr-2" => "avhrr/2".to_string(),
            "avhrr-3" => "avhrr/3".to_string(),
            _ => s.to_string(),
        }
    };
    translate(name1) == translate(name2)
}

pub fn check_and_adjust_instrument_name(platform_name: &str, instrument: &str) -> String {
    if let Some(expected) = INSTRUMENTS.get(platform_name) {
        match expected {
            InstrumentValue::Single(s) => {
                let normalized = s.replace("/", "");
                if normalized != instrument {
                    log::warn!(
                        "Inconsistent sensor/satellite input - sensor set to {}",
                        normalized
                    );
                }
                normalized
            }
            InstrumentValue::List(list) => {
                let norm_instr = instrument.replace("/", "");
                for s in list {
                    if s.replace("/", "") == norm_instr {
                        return norm_instr;
                    }
                }
                list.join("/")
            }
        }
    } else {
        instrument.to_lowercase().replace("/", "").replace("-", "")
    }
}

pub fn get_bandname_from_wavelength(
    sensor: &str,
    wavelength: f64,
    rsr: &HashMap<String, HashMap<String, RsrData>>,
    epsilon: f64,
    multiple_bands: bool,
) -> Option<Vec<String>> {
    let band_names = crate::bandnames::get_bandnames();
    let _sensor_names = band_names.get(sensor).cloned();

    let mut matches = Vec::new();
    for (band_name, detectors) in rsr {
        if let Some(det1) = detectors.get("det-1") {
            let cw = det1.central_wavelength;
            if (cw - wavelength).abs() < epsilon {
                matches.push(band_name.clone());
            }
        }
    }

    if matches.is_empty() {
        None
    } else if multiple_bands {
        Some(matches)
    } else {
        Some(vec![matches[0].clone()])
    }
}

// --- Constants from pyspectral.utils ---

pub const WAVE_LENGTH: &str = "wavelength";
pub const WAVE_NUMBER: &str = "wavenumber";

pub const HTTP_PYSPECTRAL_RSR: &str =
    "https://zenodo.org/records/19373017/files/pyspectral_rsr_data.tgz";
pub const RSR_DATA_VERSION_FILENAME: &str = "PYSPECTRAL_RSR_VERSION";
pub const RSR_DATA_VERSION: &str = "v1.6.1";

#[derive(Debug, Clone)]
pub enum InstrumentValue {
    Single(String),
    List(Vec<String>),
}

pub fn get_instruments() -> HashMap<&'static str, InstrumentValue> {
    let mut m: HashMap<&'static str, InstrumentValue> = HashMap::new();

    macro_rules! s {
        ($k:expr, $v:expr) => {
            m.insert($k, InstrumentValue::Single($v.into()));
        };
    }
    macro_rules! l { ($k:expr, $($v:expr),+) => { m.insert($k, InstrumentValue::List(vec![$($v.into()),+])); }; }

    s!("Envisat", "aatsr");
    s!("GOES-16", "abi");
    s!("GOES-17", "abi");
    s!("GOES-18", "abi");
    s!("GOES-19", "abi");
    s!("FY-4A", "agri");
    l!("FY-4B", "agri", "ghi");
    s!("Himawari-8", "ahi");
    s!("Himawari-9", "ahi");
    s!("GEO-KOMPSAT-2A", "ami");
    s!("GEO-KOMPSAT-2B", "goci-2");
    s!("NOAA-10", "avhrr/1");
    s!("NOAA-6", "avhrr/1");
    s!("NOAA-8", "avhrr/1");
    s!("TIROS-N", "avhrr/1");
    s!("NOAA-11", "avhrr/2");
    s!("NOAA-12", "avhrr/2");
    s!("NOAA-14", "avhrr/2");
    s!("NOAA-7", "avhrr/2");
    s!("NOAA-9", "avhrr/2");
    s!("Metop-A", "avhrr/3");
    s!("Metop-B", "avhrr/3");
    s!("Metop-C", "avhrr/3");
    s!("NOAA-15", "avhrr/3");
    s!("NOAA-16", "avhrr/3");
    s!("NOAA-17", "avhrr/3");
    s!("NOAA-18", "avhrr/3");
    s!("NOAA-19", "avhrr/3");
    s!("HY-1C", "cocts");
    s!("Meteosat-12", "fci");
    s!("MTG-I1", "fci");
    s!("Metop-SG-A1", "metimage");
    s!("EOS-Aqua", "modis");
    s!("EOS-Terra", "modis");
    s!("Aqua", "modis");
    s!("Terra", "modis");
    s!("Sentinel-2A", "msi");
    s!("Sentinel-2B", "msi");
    s!("Sentinel-2C", "msi");
    s!("Arctica-M-N1", "msu-gsa");
    s!("Electro-L-N2", "msu-gs");
    l!("Sentinel-3A", "olci", "slstr");
    l!("Sentinel-3B", "olci", "slstr");
    s!("Landsat-8", "oli_tirs");
    s!("Landsat-9", "oli_tirs");
    s!("Meteosat-10", "seviri");
    s!("Meteosat-11", "seviri");
    s!("Meteosat-8", "seviri");
    s!("Meteosat-9", "seviri");
    s!("NOAA-20", "viirs");
    s!("NOAA-21", "viirs");
    s!("Suomi-NPP", "viirs");
    l!("FY-3A", "virr", "mersi-1");
    l!("FY-3B", "virr", "mersi-1");
    l!("FY-3C", "virr", "mersi-1");
    s!("FY-3D", "mersi-2");
    s!("FY-3F", "mersi-3");
    s!("FY-3G", "mersi-rm");
    s!("DSCOVR", "epic");

    m
}

pub static INSTRUMENTS: once_cell::sync::Lazy<HashMap<&'static str, InstrumentValue>> =
    once_cell::sync::Lazy::new(get_instruments);

pub const AEROSOL_TYPES: &[&str] = &[
    "antarctic_aerosol",
    "continental_average_aerosol",
    "continental_clean_aerosol",
    "continental_polluted_aerosol",
    "desert_aerosol",
    "marine_clean_aerosol",
    "marine_polluted_aerosol",
    "marine_tropical_aerosol",
    "rayleigh_only",
    "rural_aerosol",
    "urban_aerosol",
];

pub const ATMOSPHERES: &[(&str, usize)] = &[
    ("subarctic_summer", 4),
    ("subarctic_winter", 5),
    ("midlatitude_summer", 6),
    ("midlatitude_winter", 7),
    ("tropical", 8),
    ("us_standard", 9),
];

#[derive(Debug, Clone)]
pub struct AtmCorrectionVersion {
    pub version: &'static str,
    pub filename: &'static str,
}

pub fn get_atm_correction_lut_version() -> HashMap<&'static str, AtmCorrectionVersion> {
    let mut m = HashMap::new();
    macro_rules! v {
        ($k:expr, $ver:expr, $fn:expr) => {
            m.insert(
                $k,
                AtmCorrectionVersion {
                    version: $ver,
                    filename: $fn,
                },
            );
        };
    }
    v!("antarctic_aerosol", "v1.0.1", "PYSPECTRAL_ATM_CORR_LUT_AA");
    v!(
        "continental_average_aerosol",
        "v1.0.1",
        "PYSPECTRAL_ATM_CORR_LUT_CAA"
    );
    v!(
        "continental_clean_aerosol",
        "v1.0.1",
        "PYSPECTRAL_ATM_CORR_LUT_CCA"
    );
    v!(
        "continental_polluted_aerosol",
        "v1.0.1",
        "PYSPECTRAL_ATM_CORR_LUT_CPA"
    );
    v!("desert_aerosol", "v1.0.1", "PYSPECTRAL_ATM_CORR_LUT_DA");
    v!(
        "marine_clean_aerosol",
        "v1.0.1",
        "PYSPECTRAL_ATM_CORR_LUT_MCA"
    );
    v!(
        "marine_polluted_aerosol",
        "v1.0.1",
        "PYSPECTRAL_ATM_CORR_LUT_MPA"
    );
    v!(
        "marine_tropical_aerosol",
        "v1.0.1",
        "PYSPECTRAL_ATM_CORR_LUT_MTA"
    );
    v!("rural_aerosol", "v1.0.1", "PYSPECTRAL_ATM_CORR_LUT_RA");
    v!("urban_aerosol", "v1.0.1", "PYSPECTRAL_ATM_CORR_LUT_UA");
    v!("rayleigh_only", "v1.0.1", "PYSPECTRAL_ATM_CORR_LUT_RO");
    m
}

pub static ATM_CORRECTION_LUT_VERSION: once_cell::sync::Lazy<
    HashMap<&'static str, AtmCorrectionVersion>,
> = once_cell::sync::Lazy::new(get_atm_correction_lut_version);

pub fn get_https_rayleigh_luts() -> HashMap<&'static str, &'static str> {
    let mut m = HashMap::new();
    let base = "https://zenodo.org/records/";
    m.insert(
        "antarctic_aerosol",
        "19372152/files/pyspectral_atm_correction_lut_aa.tgz",
    );
    m.insert(
        "continental_average_aerosol",
        "19372152/files/pyspectral_atm_correction_lut_caa.tgz",
    );
    m.insert(
        "continental_clean_aerosol",
        "19372152/files/pyspectral_atm_correction_lut_cca.tgz",
    );
    m.insert(
        "continental_polluted_aerosol",
        "19372152/files/pyspectral_atm_correction_lut_cpa.tgz",
    );
    m.insert(
        "desert_aerosol",
        "19372152/files/pyspectral_atm_correction_lut_da.tgz",
    );
    m.insert(
        "marine_clean_aerosol",
        "19372152/files/pyspectral_atm_correction_lut_mca.tgz",
    );
    m.insert(
        "marine_polluted_aerosol",
        "19372152/files/pyspectral_atm_correction_lut_mpa.tgz",
    );
    m.insert(
        "marine_tropical_aerosol",
        "19372152/files/pyspectral_atm_correction_lut_mta.tgz",
    );
    m.insert(
        "rural_aerosol",
        "19372152/files/pyspectral_atm_correction_lut_ra.tgz",
    );
    m.insert(
        "urban_aerosol",
        "19372152/files/pyspectral_atm_correction_lut_ua.tgz",
    );
    m.insert(
        "rayleigh_only",
        "19372152/files/pyspectral_atm_correction_lut_ro.tgz",
    );

    let _ = base;
    m
}

pub static HTTPS_RAYLEIGH_LUTS: once_cell::sync::Lazy<HashMap<&'static str, &'static str>> =
    once_cell::sync::Lazy::new(get_https_rayleigh_luts);

#[derive(Debug, Clone)]
pub struct RsrData {
    pub wavelength: Array1<f64>,
    pub response: Array1<f64>,
    pub central_wavelength: f64,
}

pub fn get_rayleigh_lut_dir(base_dir: &Path, aerosol_type: &str) -> PathBuf {
    base_dir.join(aerosol_type)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_instruments_contains_expected() {
        let instr = get_instruments();
        assert!(instr.contains_key("GOES-16"));
        assert!(instr.contains_key("Meteosat-10"));
        assert!(instr.contains_key("Sentinel-3A"));
    }

    #[test]
    fn test_instruments_list() {
        let instr = get_instruments();
        match instr.get("Sentinel-3A").unwrap() {
            InstrumentValue::List(l) => {
                assert!(l.contains(&"olci".to_string()));
                assert!(l.contains(&"slstr".to_string()));
            }
            _ => panic!("Expected list"),
        }
    }

    #[test]
    fn test_check_and_adjust_known() {
        // platform=GOES-16 instrument=abi should return "abi"
        let result = check_and_adjust_instrument_name("GOES-16", "abi");
        assert_eq!(result, "abi");
    }

    #[test]
    fn test_check_and_adjust_unknown() {
        let result = check_and_adjust_instrument_name("Unknown", "MyInst/2");
        assert_eq!(result, "myinst2");
    }

    #[test]
    fn test_aerosol_types_count() {
        assert_eq!(AEROSOL_TYPES.len(), 11);
    }

    #[test]
    fn test_atmospheres_count() {
        assert_eq!(ATMOSPHERES.len(), 6);
    }

    #[test]
    fn test_lut_version_entries() {
        let v = get_atm_correction_lut_version();
        assert_eq!(v.len(), 11);
    }

    #[test]
    fn test_https_rayleigh_luts_entries() {
        let m = get_https_rayleigh_luts();
        assert_eq!(m.len(), 11);
    }

    fn test_rsr_data() -> (ndarray::Array1<f64>, ndarray::Array1<f64>) {
        let wvl = ndarray::arr1(&[
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
        ]);
        let resp = ndarray::arr1(&[
            0.01, 0.0118, 0.01987, 0.03226, 0.05028, 0.0849, 0.16645, 0.33792, 0.59106, 0.81815,
            0.96077, 0.92855, 0.86008, 0.8661, 0.87697, 0.85412, 0.88922, 0.9541, 0.95687, 0.91037,
            0.91058, 0.94256, 0.94719, 0.94808, 1.0, 0.92676, 0.67429, 0.44715, 0.27762, 0.14852,
            0.07141, 0.04151, 0.02925, 0.02085, 0.01414, 0.01_f64,
        ]);
        (wvl, resp)
    }

    fn ch3_rsr() -> (ndarray::Array1<f64>, ndarray::Array1<f64>) {
        let wvl = ndarray::arr1(&[
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
        ]);
        let resp = ndarray::arr1(&[
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
        ]);
        (wvl, resp)
    }

    #[test]
    fn test_trapezoid() {
        let x = ndarray::arr1(&[0.0, 1.0, 2.0, 3.0]);
        let y = ndarray::arr1(&[0.0, 1.0, 4.0, 9.0]);
        let result = trapezoid(&y, &x);
        approx::assert_abs_diff_eq!(result, 9.5, epsilon = 1e-10);
    }

    #[test]
    fn test_central_wavelength() {
        let (wvl, resp) = test_rsr_data();
        let cw = get_central_wave(&wvl, &resp, &Array1::from_elem(wvl.len(), 1.0));
        approx::assert_abs_diff_eq!(cw, 3.780_281_935, epsilon = 1e-6);
    }

    #[test]
    fn test_convert2wavenumber() {
        let (wvl, resp) = test_rsr_data();
        let (wnum, wresp) = convert2wavenumber_rsr(&wvl, &resp);
        let expected_first = 1.0 / (1e-4 * wvl[35]);
        let expected_last = 1.0 / (1e-4 * wvl[0]);
        approx::assert_abs_diff_eq!(wnum[0], expected_first, epsilon = 1e-3);
        approx::assert_abs_diff_eq!(wnum[35], expected_last, epsilon = 1e-3);
        approx::assert_abs_diff_eq!(wresp[0], resp[35], epsilon = 1e-10);
        approx::assert_abs_diff_eq!(wresp[35], resp[0], epsilon = 1e-10);
    }

    #[test]
    fn test_sort_data() {
        let x = ndarray::arr1(&[1.0, 5.6, 30.0, 2.1, 108.2, 57.8, 1e9, 2.1_f64]);
        let y = ndarray::arr1(&[45.0, 92.0, 20.0, 10.0, 15.0, 67.0, 108.0, 15.0_f64]);
        let (x_sorted, y_sorted) = sort_data(&x, &y);
        let expected_x = ndarray::arr1(&[1.0, 2.1, 5.6, 30.0, 57.8, 108.2, 1e9_f64]);
        let expected_y = ndarray::arr1(&[45.0, 10.0, 92.0, 20.0, 67.0, 15.0, 108.0_f64]);
        for i in 0..expected_x.len() {
            approx::assert_abs_diff_eq!(x_sorted[i], expected_x[i], epsilon = 1e-10);
            approx::assert_abs_diff_eq!(y_sorted[i], expected_y[i], epsilon = 1e-10);
        }
    }

    #[test]
    fn test_sort_data_no_duplicates() {
        let x = ndarray::arr1(&[3.0, 1.0, 2.0_f64]);
        let y = ndarray::arr1(&[30.0, 10.0, 20.0_f64]);
        let (x_sorted, y_sorted) = sort_data(&x, &y);
        let expected_x = ndarray::arr1(&[1.0, 2.0, 3.0_f64]);
        let expected_y = ndarray::arr1(&[10.0, 20.0, 30.0_f64]);
        for i in 0..3 {
            approx::assert_abs_diff_eq!(x_sorted[i], expected_x[i]);
            approx::assert_abs_diff_eq!(y_sorted[i], expected_y[i]);
        }
    }

    #[test]
    fn test_fwhm() {
        let (wvl, resp) = ch3_rsr();
        let fwhm = get_fullwidth_halfmax(&resp, &wvl);
        approx::assert_abs_diff_eq!(fwhm, 0.065_818_01, epsilon = 1e-5);
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
        approx::assert_abs_diff_eq!(low, 0.594_813_23, epsilon = 1e-5);
        approx::assert_abs_diff_eq!(high, 0.657_375_75, epsilon = 1e-5);
    }

    #[test]
    fn test_integrated_energy_ten_percent() {
        let (wvl, resp) = ch3_rsr();
        let (low, high) = get_bounds_integrated_energy(&resp, &wvl, 10.0);
        approx::assert_abs_diff_eq!(low, 0.609_310_27, epsilon = 1e-5);
        approx::assert_abs_diff_eq!(high, 0.657_375_75, epsilon = 1e-5);
    }

    #[test]
    fn test_get_wave_range() {
        let (wvl, resp) = ch3_rsr();
        let (min_wvl, cwl, max_wvl) = get_wave_range(&wvl, &resp, 0.15);
        approx::assert_abs_diff_eq!(min_wvl, 0.594_813_23, epsilon = 1e-5);
        approx::assert_abs_diff_eq!(max_wvl, 0.675_128_28, epsilon = 1e-5);
        assert!(cwl > min_wvl && cwl < max_wvl);
    }

    #[test]
    fn test_get_wave_range_higher_threshold() {
        let (wvl, resp) = ch3_rsr();
        let (min_wvl, _cwl, max_wvl) = get_wave_range(&wvl, &resp, 0.5);
        approx::assert_abs_diff_eq!(min_wvl, 0.609_310_27, epsilon = 1e-5);
        approx::assert_abs_diff_eq!(max_wvl, 0.675_128_28, epsilon = 1e-5);
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
}
