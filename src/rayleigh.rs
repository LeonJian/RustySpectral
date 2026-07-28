use ndarray::{Array1, Array3, Array4, Axis};
use std::fs;
use std::path::{Path, PathBuf};

use log::{debug, info, warn};

use crate::config::get_config;
use crate::download::{download_luts, get_rayleigh_lut_dir};
use crate::rsr_reader::RelativeSpectralResponse;
use crate::utils::{
    self, get_central_wave, AEROSOL_TYPES, ATMOSPHERES, ATM_CORRECTION_LUT_VERSION,
};

struct CachedLutData {
    reflectance_4d: Array4<f64>,
    wvl_coord: Array1<f64>,
    azid_coord: Array1<f64>,
    satz_sec_coord: Array1<f64>,
    sunz_sec_coord: Array1<f64>,
}

pub fn clip_angles_inside_coordinate_range(
    zenith_angle: &Array1<f64>,
    zenith_secant_max: f64,
) -> Array1<f64> {
    let clip_angle = (1.0 / zenith_secant_max).acos().to_degrees();
    zenith_angle.mapv(|z| {
        if z.is_nan() {
            0.0
        } else {
            z.clamp(0.0, clip_angle)
        }
    })
}

#[inline]
pub fn clip_angles_inside_coordinate_range_scalar(
    zenith_angle: f64,
    zenith_secant_max: f64,
) -> f64 {
    if zenith_angle.is_nan() {
        0.0f64
    } else {
        let clip_angle = (1.0 / zenith_secant_max).acos().to_degrees();
        zenith_angle.clamp(0.0, clip_angle)
    }
}

pub fn reduce_rayleigh_highzenith(
    zenith: &Array1<f64>,
    rayref: &Array1<f64>,
    thresh_zen: f64,
    maxzen: f64,
    strength: f64,
) -> Array1<f64> {
    let factor: Array1<f64> = zenith.mapv(|z| {
        if z < thresh_zen {
            0.0
        } else {
            (z - thresh_zen) / (maxzen - thresh_zen)
        }
    });
    let factor = 1.0 - strength * &factor;
    let factor = factor.mapv(|f| f.clamp(0.0, 1.0));
    rayref * &factor
}

#[inline]
pub fn get_wavelength_index_and_factor(wvl_coord: &Array1<f64>, wvl: f64) -> (usize, f64) {
    let idx = match wvl_coord.iter().position(|&v| v > wvl) {
        Some(0) => 1,
        Some(i) => i,
        None => wvl_coord.len() - 1,
    };
    let wavelength_index = idx;
    let wvl1 = wvl_coord[wavelength_index - 1];
    let wvl2 = wvl_coord[wavelength_index];
    let wavelength_factor = (wvl2 - wvl) / (wvl2 - wvl1);
    (wavelength_index, wavelength_factor)
}

pub fn get_wavelength_adjusted_lut(
    rayleigh_refl: &Array4<f64>,
    wvl_coord: &Array1<f64>,
    wvl: f64,
) -> Array3<f64> {
    let (wi, wf) = get_wavelength_index_and_factor(wvl_coord, wvl);
    let slice1 = rayleigh_refl.index_axis(Axis(0), wi - 1);
    let slice2 = rayleigh_refl.index_axis(Axis(0), wi);
    wf * &slice1 + (1.0 - wf) * &slice2
}

pub fn trilinear_interpolate(
    grid: &Array3<f64>,
    sunz_sec: f64,
    azidiff_in: f64,
    satz_sec: f64,
    sunz_coord: &Array1<f64>,
    azid_coord: &Array1<f64>,
    satz_coord: &Array1<f64>,
) -> f64 {
    let azidiff = 180.0 - azidiff_in;

    let si = find_interval_index(sunz_coord, sunz_sec);
    let ai = find_interval_index(azid_coord, azidiff);
    let ti = find_interval_index(satz_coord, satz_sec);

    let s0 = sunz_coord[si];
    let s1 = sunz_coord[si + 1];
    let a0 = azid_coord[ai];
    let a1 = azid_coord[ai + 1];
    let t0 = satz_coord[ti];
    let t1 = satz_coord[ti + 1];

    let sd = (sunz_sec - s0) / (s1 - s0);
    let ad = (azidiff - a0) / (a1 - a0);
    let td = (satz_sec - t0) / (t1 - t0);

    let c000 = grid[(si, ai, ti)];
    let c001 = grid[(si, ai, ti + 1)];
    let c010 = grid[(si, ai + 1, ti)];
    let c011 = grid[(si, ai + 1, ti + 1)];
    let c100 = grid[(si + 1, ai, ti)];
    let c101 = grid[(si + 1, ai, ti + 1)];
    let c110 = grid[(si + 1, ai + 1, ti)];
    let c111 = grid[(si + 1, ai + 1, ti + 1)];

    let c00 = c000 * (1.0 - td) + c001 * td;
    let c01 = c010 * (1.0 - td) + c011 * td;
    let c10 = c100 * (1.0 - td) + c101 * td;
    let c11 = c110 * (1.0 - td) + c111 * td;

    let c0 = c00 * (1.0 - ad) + c01 * ad;
    let c1 = c10 * (1.0 - ad) + c11 * ad;

    (c0 * (1.0 - sd) + c1 * sd) * 100.0
}

#[inline]
fn find_interval_index(coords: &Array1<f64>, value: f64) -> usize {
    let idx = coords
        .iter()
        .position(|&v| v > value)
        .unwrap_or(coords.len() - 1);
    (idx.saturating_sub(1)).min(coords.len() - 2)
}

#[allow(clippy::too_many_arguments)]
pub fn rayleigh_interpolate_by_angles(
    sun_zenith: &Array1<f64>,
    sat_zenith: &Array1<f64>,
    azidiff: &Array1<f64>,
    rayleigh_refl: &Array4<f64>,
    wvl_coord: &Array1<f64>,
    wvl: f64,
    sunz_sec_coord: &Array1<f64>,
    satz_sec_coord: &Array1<f64>,
    azid_coord: &Array1<f64>,
) -> Array1<f64> {
    let grid3 = get_wavelength_adjusted_lut(rayleigh_refl, wvl_coord, wvl);

    let n = sun_zenith.len();
    let mut result = Array1::zeros(n);
    for i in 0..n {
        let sz = clip_angles_inside_coordinate_range_scalar(
            sun_zenith[i],
            sunz_sec_coord[sunz_sec_coord.len() - 1],
        );
        let satz = clip_angles_inside_coordinate_range_scalar(
            sat_zenith[i],
            satz_sec_coord[satz_sec_coord.len() - 1],
        );
        let sunzsec = 1.0 / sz.to_radians().cos();
        let satzsec = 1.0 / satz.to_radians().cos();

        result[i] = trilinear_interpolate(
            &grid3,
            sunzsec,
            azidiff[i],
            satzsec,
            sunz_sec_coord,
            azid_coord,
            satz_sec_coord,
        );
    }
    result
}

pub fn normalize_sensor(platform_name: &str, sensor: &str) -> String {
    let instruments = utils::get_instruments();
    let instr = match instruments.get(platform_name) {
        Some(val) => match val {
            utils::InstrumentValue::Single(s) => {
                if s != sensor {
                    warn!("Inconsistent sensor/satellite input - sensor set to {}", s);
                }
                s.clone()
            }
            utils::InstrumentValue::List(list) => {
                if !list.contains(&sensor.to_string()) {
                    panic!(
                        "This satellite has multiple sensors, you must explicitly state which to use."
                    );
                }
                sensor.to_string()
            }
        },
        None => sensor.to_string(),
    };
    instr.replace("/", "")
}

fn normalize_atmosphere(name: &str) -> String {
    name.replace(" ", "_").replace("-", "_")
}

pub struct RayleighConfigBase {
    pub aerosol_type: String,
    pub atm_type: String,
    pub do_download: bool,
    pub lutfiles_version_uptodate: bool,
}

impl RayleighConfigBase {
    pub fn new(aerosol_type: &str, atm_type: &str) -> Self {
        let config = get_config(None);

        let atm_normalized = normalize_atmosphere(atm_type);
        let atm_valid = ATMOSPHERES
            .iter()
            .any(|(name, _)| normalize_atmosphere(name) == atm_normalized);
        if !atm_valid {
            panic!(
                "Atmosphere type not supported! Need to be one of {:?}",
                ATMOSPHERES.iter().map(|(n, _)| *n).collect::<Vec<_>>()
            );
        }

        let aerosol_valid = AEROSOL_TYPES.contains(&aerosol_type);
        if !aerosol_valid {
            panic!(
                "Aerosol type not supported! Need to be one of {:?}",
                AEROSOL_TYPES
            );
        }

        let lutfiles_version_uptodate = Self::check_version(aerosol_type);

        RayleighConfigBase {
            aerosol_type: aerosol_type.to_string(),
            atm_type: atm_normalized,
            do_download: config.download_from_internet,
            lutfiles_version_uptodate,
        }
    }

    fn check_version(aerosol_type: &str) -> bool {
        let config = get_config(None);
        let lut_dir = get_rayleigh_lut_dir(&config, aerosol_type);

        let ver_info = match ATM_CORRECTION_LUT_VERSION.get(aerosol_type) {
            Some(v) => v,
            None => return false,
        };

        let version_file = lut_dir.join(ver_info.filename);
        if !version_file.exists() {
            return false;
        }

        let current = fs::read_to_string(&version_file)
            .unwrap_or_default()
            .trim()
            .to_string();
        current == ver_info.version
    }
}

pub struct Rayleigh {
    pub base: RayleighConfigBase,
    pub sensor: String,
    pub platform_name: String,
    pub reflectance_lut_filename: PathBuf,
    lut_data: CachedLutData,
}

impl Rayleigh {
    pub fn new(
        platform_name: &str,
        sensor: &str,
        atmosphere: Option<&str>,
        aerosol_type: Option<&str>,
    ) -> Self {
        let atm_type = atmosphere.unwrap_or("us_standard");
        let aer_type = aerosol_type.unwrap_or("marine_clean_aerosol");

        let base = RayleighConfigBase::new(aer_type, atm_type);
        let sensor_norm = normalize_sensor(platform_name, sensor);

        let config = get_config(None);
        let rayleigh_dir = get_rayleigh_lut_dir(&config, aer_type);
        let ext = atm_type.replace("_", " ");
        let lutname = format!("rayleigh_lut_{}.h5", ext.replace(" ", "_"));
        let reflectance_lut_filename = rayleigh_dir.join(&lutname);

        debug!("LUT filename: {}", reflectance_lut_filename.display());

        if !base.lutfiles_version_uptodate && base.do_download {
            info!("Rayleigh LUT files not up to date, will download from internet...");
            let types = vec![aer_type.to_string()];
            if let Err(e) = download_luts(Some(&types), false) {
                warn!("Failed to download LUTs: {}", e);
            }
        }

        let lut_data = Self::load_lut_data(&reflectance_lut_filename);

        Rayleigh {
            base,
            sensor: sensor_norm,
            platform_name: platform_name.to_string(),
            reflectance_lut_filename,
            lut_data,
        }
    }

    fn load_lut_data(lut_filename: &Path) -> CachedLutData {
        let h5file = hdf5_pure::File::open(lut_filename).expect("Failed to open Rayleigh LUT file");
        let root = h5file.root();

        let read_coord = |name: &str| -> Array1<f64> {
            let ds = root
                .dataset(name)
                .unwrap_or_else(|_| panic!("Dataset '{}' not found in LUT file", name));
            Array1::from_vec(
                ds.read_f64()
                    .unwrap_or_else(|_| panic!("Failed to read '{}'", name)),
            )
        };

        let azid_coord = read_coord("azimuth_difference");
        let satz_sec_coord = read_coord("satellite_zenith_secant");
        let sunz_sec_coord = read_coord("sun_zenith_secant");
        let wvl_coord = read_coord("wavelengths");

        let refl_ds = root
            .dataset("reflectance")
            .expect("Dataset 'reflectance' not found in LUT file");
        let shape = refl_ds.shape().expect("Failed to get reflectance shape");
        let (nw, ns, na, nt) = (
            shape[0] as usize,
            shape[1] as usize,
            shape[2] as usize,
            shape[3] as usize,
        );
        let refl_values: Vec<f64> = refl_ds
            .read_f64()
            .expect("Failed to read reflectance dataset");
        let reflectance_4d = Array4::from_shape_vec((nw, ns, na, nt), refl_values)
            .expect("Reflectance shape mismatch");

        CachedLutData {
            reflectance_4d,
            wvl_coord,
            azid_coord,
            satz_sec_coord,
            sunz_sec_coord,
        }
    }

    pub fn get_reflectance(
        &self,
        sun_zenith: &Array1<f64>,
        sat_zenith: &Array1<f64>,
        azidiff: &Array1<f64>,
        band_name_or_wavelength: &str,
        redband: Option<&Array1<f64>>,
    ) -> Array1<f64> {
        let wvl_nm: f64;
        let band_name: String;

        if let Ok(wvl_um) = band_name_or_wavelength.parse::<f64>() {
            warn!(
                "A wavelength is provided instead of band name - disregard the RSRs. \
                 Effective wavelength: {} (micro meter)",
                wvl_um
            );
            wvl_nm = wvl_um * 1000.0;
            let _band_name = format!("{:.6}um", wvl_um);
        } else {
            band_name = band_name_or_wavelength.to_string();

            let cwvl = match self.get_rsr_wavelength_from_band_name(&band_name) {
                Some(w) => {
                    debug!("Band name: {}  Effective wavelength: {}um", band_name, w);
                    w * 1000.0
                }
                None => {
                    warn!(
                        "Effective wavelength for band {} outside nominal 400-800 nm range!",
                        band_name
                    );
                    info!("Setting the rayleigh/aerosol reflectance contribution to zero!");
                    let n = sun_zenith.len();
                    return Array1::zeros(n);
                }
            };
            wvl_nm = cwvl;
        }

        let n = sun_zenith.len();
        let mut result = Array1::zeros(n);

        match self.interp_rayleigh_refl_by_angles(sun_zenith, sat_zenith, azidiff, wvl_nm) {
            Ok(res) => {
                let mut final_res = res;
                if let Some(rb) = redband {
                    final_res = self.relax_rayleigh_refl_correction_where_cloudy(rb, &final_res);
                }
                final_res = final_res.mapv(|v| v.clamp(0.0, 100.0));
                result = final_res;
            }
            Err(_) => {
                warn!("Failed to interpolate Rayleigh reflectance, returning zeros.");
            }
        }

        result
    }

    fn get_rsr_wavelength_from_band_name(&self, band_name: &str) -> Option<f64> {
        let rsr = match RelativeSpectralResponse::new(
            Some(&self.platform_name),
            Some(&self.sensor),
            None,
        ) {
            Ok(r) => r,
            Err(e) => {
                warn!(
                    "No spectral responses for platform {} and sensor {}: {}",
                    self.platform_name, self.sensor, e
                );
                return None;
            }
        };

        if let Some(detectors) = rsr.rsr.get(band_name) {
            if let Some(det1) = detectors.get("det-1") {
                let wvl = &det1.wavelength;
                let resp = &det1.response;
                let weight: Array1<f64> = wvl.mapv(|w| 1.0 / w.powi(4));
                let cwvl = get_central_wave(wvl, resp, &weight);
                return Some(cwvl);
            }
        }

        let band_names = crate::bandnames::get_bandnames();
        if let Some(sensor_names) = band_names.get(self.sensor.as_str()) {
            let mapped = sensor_names.get(band_name).copied().unwrap_or(band_name);
            if let Some(detectors) = rsr.rsr.get(mapped) {
                if let Some(det1) = detectors.get("det-1") {
                    let wvl = &det1.wavelength;
                    let resp = &det1.response;
                    let weight: Array1<f64> = wvl.mapv(|w| 1.0 / w.powi(4));
                    let cwvl = get_central_wave(wvl, resp, &weight);
                    return Some(cwvl);
                }
            }
        }

        None
    }
    fn interp_rayleigh_refl_by_angles(
        &self,
        sun_zenith: &Array1<f64>,
        sat_zenith: &Array1<f64>,
        azidiff: &Array1<f64>,
        wvl_nm: f64,
    ) -> Result<Array1<f64>, String> {
        let lut = &self.lut_data;
        let grid3 = get_wavelength_adjusted_lut(&lut.reflectance_4d, &lut.wvl_coord, wvl_nm);

        let n = sun_zenith.len();
        let mut result = Array1::zeros(n);

        let sunz_sec_max = lut.sunz_sec_coord[lut.sunz_sec_coord.len() - 1];
        let satz_sec_max = lut.satz_sec_coord[lut.satz_sec_coord.len() - 1];

        for i in 0..n {
            let sz = clip_angles_inside_coordinate_range_scalar(sun_zenith[i], sunz_sec_max);
            let satz = clip_angles_inside_coordinate_range_scalar(sat_zenith[i], satz_sec_max);
            let sunzsec = 1.0 / sz.to_radians().cos().max(0.0001);
            let satzsec = 1.0 / satz.to_radians().cos().max(0.0001);

            result[i] = trilinear_interpolate(
                &grid3,
                sunzsec,
                azidiff[i],
                satzsec,
                &lut.sunz_sec_coord,
                &lut.azid_coord,
                &lut.satz_sec_coord,
            );
        }

        Ok(result)
    }

    fn relax_rayleigh_refl_correction_where_cloudy(
        &self,
        redband: &Array1<f64>,
        rayleigh_refl: &Array1<f64>,
    ) -> Array1<f64> {
        let n = redband.len();
        let mut result = Array1::zeros(n);
        for i in 0..n {
            let rb = redband[i];
            let rr = rayleigh_refl[i];
            if rb < 20.0 {
                result[i] = rr;
            } else {
                result[i] = (1.0 - (rb - 20.0) / 80.0) * rr;
            }
        }
        result
    }
}

#[allow(clippy::type_complexity)]
pub fn get_reflectance_lut_from_file(
    lut_filename: &Path,
) -> Result<(Array1<f64>, Array1<f64>, Array1<f64>), String> {
    if !lut_filename.exists() {
        return Err(format!(
            "Rayleigh LUT file does not exist! Filename = {}",
            lut_filename.display()
        ));
    }

    let h5file = hdf5_pure::File::open(lut_filename)
        .map_err(|e| format!("Failed to open LUT file: {}", e))?;

    let root = h5file.root();

    let azidiff = read_lut_coord(&root, "azimuth_difference")?;
    let satz_sec = read_lut_coord(&root, "satellite_zenith_secant")?;
    let sunz_sec = read_lut_coord(&root, "sun_zenith_secant")?;

    Ok((azidiff, satz_sec, sunz_sec))
}

fn read_lut_coord(root: &hdf5_pure::Group, name: &str) -> Result<Array1<f64>, String> {
    let ds = root
        .dataset(name)
        .map_err(|e| format!("Dataset '{}' not found: {}", name, e))?;
    let values: Vec<f64> = ds
        .read_f64()
        .map_err(|e| format!("Failed to read '{}': {}", name, e))?;
    Ok(Array1::from_vec(values))
}

pub fn read_reflectance_lut_4d(lut_filename: &Path) -> Result<Array4<f64>, String> {
    let h5file = hdf5_pure::File::open(lut_filename)
        .map_err(|e| format!("Failed to open LUT file: {}", e))?;
    let root = h5file.root();
    let ds = root
        .dataset("reflectance")
        .map_err(|e| format!("Dataset 'reflectance' not found: {}", e))?;
    let shape = ds
        .shape()
        .map_err(|e| format!("Failed to get shape: {}", e))?;

    if shape.len() != 4 {
        return Err(format!("Expected 4D reflectance, got {}D", shape.len()));
    }
    let (nw, ns, na, nt) = (
        shape[0] as usize,
        shape[1] as usize,
        shape[2] as usize,
        shape[3] as usize,
    );
    let values: Vec<f64> = ds
        .read_f64()
        .map_err(|e| format!("Failed to read reflectance: {}", e))?;

    let arr = Array4::from_shape_vec((nw, ns, na, nt), values)
        .map_err(|e| format!("Shape mismatch: {}", e))?;
    Ok(arr)
}

pub fn read_wavelength_lut_coord(lut_filename: &Path) -> Result<Array1<f64>, String> {
    let h5file = hdf5_pure::File::open(lut_filename)
        .map_err(|e| format!("Failed to open LUT file: {}", e))?;
    let root = h5file.root();
    read_lut_coord(&root, "wavelengths")
}

pub fn check_and_download(dry_run: bool, aerosol_types: Option<&[String]>) {
    let types: Vec<String> = aerosol_types
        .map(|v| v.to_vec())
        .unwrap_or_else(|| AEROSOL_TYPES.iter().map(|s| s.to_string()).collect());

    let mut needed: Vec<String> = Vec::new();
    for aerosol_type in &types {
        let base = RayleighConfigBase::new(aerosol_type, "us_standard");
        if base.lutfiles_version_uptodate {
            info!(
                "Atm correction LUTs for {} already the latest!",
                aerosol_type
            );
        } else {
            needed.push(aerosol_type.clone());
        }
    }

    if !needed.is_empty() {
        info!("Downloading LUTs for: {:?}", needed);
        if let Err(e) = download_luts(Some(&needed), dry_run) {
            warn!("Failed to download LUTs: {}", e);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use ndarray::arr1;

    #[test]
    fn test_clip_angles_with_nans() {
        let zenith = arr1(&[f64::NAN, 30.0, 45.0]);
        let result = clip_angles_inside_coordinate_range(&zenith, 2.5);
        assert_eq!(result[0], 0.0);
        assert!(result[1] > 0.0);
    }

    #[test]
    fn test_clip_angles_scalar() {
        let result = clip_angles_inside_coordinate_range_scalar(100.0, 2.0);
        assert!(result < 90.0);
        assert!(result > 0.0);
    }

    #[test]
    fn test_reduce_rayleigh_highzenith() {
        let zenith = arr1(&[60.0, 70.0, 80.0]);
        let rayref = arr1(&[10.0, 10.0, 10.0]);
        let result = reduce_rayleigh_highzenith(&zenith, &rayref, 70.0, 90.0, 1.0);
        assert!(result[0] >= 9.0);
        assert!(result[2] <= 10.0);
    }

    #[test]
    fn test_get_wavelength_index_and_factor() {
        let coords = arr1(&[400.0, 500.0, 600.0]);
        let (idx, factor) = get_wavelength_index_and_factor(&coords, 450.0);
        assert_eq!(idx, 1);
        assert!(factor > 0.0 && factor < 1.0);
    }

    #[test]
    fn test_normalize_sensor_known() {
        let name = normalize_sensor("GOES-16", "abi");
        assert_eq!(name, "abi");
    }

    #[test]
    fn test_normalize_sensor_with_slash() {
        let name = normalize_sensor("Metop-A", "avhrr/3");
        assert_eq!(name, "avhrr3");
    }

    #[test]
    fn test_aerosol_types() {
        assert_eq!(AEROSOL_TYPES.len(), 11);
    }

    #[test]
    fn test_atmospheres() {
        assert_eq!(ATMOSPHERES.len(), 6);
    }

    #[test]
    fn test_rayleigh_config_base_valid() {
        let base = RayleighConfigBase::new("marine_clean_aerosol", "us_standard");
        assert_eq!(base.aerosol_type, "marine_clean_aerosol");
    }

    #[test]
    #[should_panic]
    fn test_rayleigh_config_base_invalid_aerosol() {
        RayleighConfigBase::new("nonexistent_aerosol", "us_standard");
    }

    #[test]
    #[should_panic]
    fn test_rayleigh_config_base_invalid_atmosphere() {
        RayleighConfigBase::new("marine_clean_aerosol", "nonexistent_atm");
    }

    #[test]
    fn test_get_reflectance_lut_file_not_found() {
        let result = get_reflectance_lut_from_file(Path::new("/nonexistent/path.h5"));
        assert!(result.is_err());
    }

    #[test]
    fn test_reduce_rayleigh_no_reduction() {
        let sun_zenith = arr1(&[70.0, 65.0, 60.0]);
        let in_rayleigh = arr1(&[50.0, 50.0, 50.0]);
        let result = reduce_rayleigh_highzenith(&sun_zenith, &in_rayleigh, 70.0, 90.0, 1.0);
        for i in 0..3 {
            assert_abs_diff_eq!(result[i], in_rayleigh[i], epsilon = 1e-6);
        }
    }

    #[test]
    fn test_reduce_rayleigh_moderate() {
        let sun_zenith = arr1(&[70.0, 65.0, 60.0]);
        let in_rayleigh = arr1(&[50.0, 50.0, 50.0]);
        let result = reduce_rayleigh_highzenith(&sun_zenith, &in_rayleigh, 30.0, 90.0, 1.0);
        let expected = arr1(&[16.666_666_67, 20.833_333_33, 25.0]);
        for i in 0..3 {
            assert_abs_diff_eq!(result[i], expected[i], epsilon = 1e-3);
        }
    }

    #[test]
    fn test_reduce_rayleigh_extreme() {
        let sun_zenith = arr1(&[70.0, 65.0, 60.0]);
        let in_rayleigh = arr1(&[50.0, 50.0, 50.0]);
        let result = reduce_rayleigh_highzenith(&sun_zenith, &in_rayleigh, 30.0, 90.0, 1.5);
        let expected = arr1(&[0.0, 6.25, 12.5]);
        for i in 0..3 {
            assert_abs_diff_eq!(result[i], expected[i], epsilon = 1e-3);
        }
    }

    #[test]
    fn test_wavelength_index_factor() {
        let wvl_coord = arr1(&[631.0, 636.0]);
        let (idx, factor) = get_wavelength_index_and_factor(&wvl_coord, 634.0);
        assert_eq!(idx, 1);
        assert_abs_diff_eq!(factor, (636.0 - 634.0) / (636.0 - 631.0), epsilon = 1e-10);
    }

    #[test]
    fn test_normalize_sensor_full() {
        assert_eq!(normalize_sensor("GOES-16", "abi"), "abi");
        assert_eq!(normalize_sensor("NOAA-19", "avhrr/3"), "avhrr3");
        assert_eq!(normalize_sensor("FY-4A", "agri"), "agri");
        assert_eq!(normalize_sensor("Himawari-8", "ahi"), "ahi");
        assert_eq!(normalize_sensor("NOAA-20", "viirs"), "viirs");
        assert_eq!(normalize_sensor("Meteosat-9", "seviri"), "seviri");
    }

    #[test]
    fn test_normalize_sensor_unknown_platform() {
        assert_eq!(normalize_sensor("Unknown", "myinstr"), "myinstr");
    }

    #[test]
    fn test_aerosol_types_check() {
        let types = crate::utils::AEROSOL_TYPES;
        assert_eq!(types.len(), 11);
        assert!(types.contains(&"marine_clean_aerosol"));
        assert!(types.contains(&"desert_aerosol"));
        assert!(types.contains(&"rayleigh_only"));
    }

    #[test]
    fn test_atmospheres_check() {
        let atms = crate::utils::ATMOSPHERES;
        assert_eq!(atms.len(), 6);
        let names: Vec<&str> = atms.iter().map(|(n, _)| *n).collect();
        assert!(names.contains(&"midlatitude_summer"));
        assert!(names.contains(&"us_standard"));
    }

    #[test]
    fn test_get_wavelength_index_and_factor_edge_cases() {
        let coords = arr1(&[400.0, 500.0, 600.0]);
        let (idx, factor) = get_wavelength_index_and_factor(&coords, 350.0);
        assert_eq!(idx, 1);
        assert_abs_diff_eq!(factor, 1.5, epsilon = 1e-10);

        let (idx, factor) = get_wavelength_index_and_factor(&coords, 650.0);
        assert_eq!(idx, 2);
        assert_abs_diff_eq!(factor, -0.5, epsilon = 1e-10);
    }

    #[test]
    fn test_trilinear_interpolate_corner() {
        let grid = Array3::from_shape_fn((2, 2, 2), |(i, j, k)| (i + j + k) as f64);
        let sunz_coord = arr1(&[1.0, 2.0]);
        let azid_coord = arr1(&[0.0, 180.0]);
        let satz_coord = arr1(&[1.0, 2.0]);

        // sunz=1.0 (sd=0), azidiff_in=0 => azidiff=180 (ad=1), satz=1.0 (td=0)
        // maps to grid[0,1,0] = 1 * 100 = 100
        let result =
            trilinear_interpolate(&grid, 1.0, 0.0, 1.0, &sunz_coord, &azid_coord, &satz_coord);
        assert_abs_diff_eq!(result, 100.0, epsilon = 1e-10);

        // sunz=2.0 (sd=1), azidiff_in=0 => azidiff=180 (ad=1), satz=2.0 (td=1)
        // maps to grid[1,1,1] = 3 * 100 = 300
        let result =
            trilinear_interpolate(&grid, 2.0, 0.0, 2.0, &sunz_coord, &azid_coord, &satz_coord);
        assert_abs_diff_eq!(result, 300.0, epsilon = 1e-10);
    }

    #[test]
    fn test_trilinear_interpolate_midpoint() {
        let grid = Array3::from_shape_fn((2, 2, 2), |(i, j, k)| (i + j + k) as f64);
        let sunz_coord = arr1(&[1.0, 2.0]);
        let azid_coord = arr1(&[0.0, 180.0]);
        let satz_coord = arr1(&[1.0, 2.0]);

        let result =
            trilinear_interpolate(&grid, 1.5, 90.0, 1.5, &sunz_coord, &azid_coord, &satz_coord);
        assert_abs_diff_eq!(result, 150.0, epsilon = 1e-10);
    }

    #[test]
    fn test_find_interval_index() {
        let coords = arr1(&[1.0, 2.0, 3.0, 4.0, 5.0]);
        assert_eq!(find_interval_index(&coords, 0.5), 0);
        assert_eq!(find_interval_index(&coords, 1.5), 0);
        assert_eq!(find_interval_index(&coords, 2.5), 1);
        assert_eq!(find_interval_index(&coords, 5.5), 3);
    }

    #[test]
    fn test_clip_angles_scalar_nan() {
        let result = clip_angles_inside_coordinate_range_scalar(f64::NAN, 2.0);
        assert_eq!(result, 0.0);
    }

    #[test]
    fn test_reduce_rayleigh_highzenith_full_reduction() {
        let zenith = arr1(&[85.0, 86.0, 87.0]);
        let rayref = arr1(&[50.0, 50.0, 50.0]);
        let result = reduce_rayleigh_highzenith(&zenith, &rayref, 85.0, 90.0, 1.0);
        let expected = arr1(&[50.0, 40.0, 30.0]);
        for i in 0..3 {
            assert_abs_diff_eq!(result[i], expected[i], epsilon = 1e-3);
        }
    }

    #[test]
    fn test_reduce_rayleigh_highzenith_below_threshold() {
        let zenith = arr1(&[10.0, 20.0, 30.0]);
        let rayref = arr1(&[10.0, 20.0, 30.0]);
        let result = reduce_rayleigh_highzenith(&zenith, &rayref, 40.0, 90.0, 1.0);
        for i in 0..3 {
            assert_abs_diff_eq!(result[i], rayref[i], epsilon = 1e-6);
        }
    }
}
