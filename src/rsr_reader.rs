use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

use log::{debug, info, warn};
use ndarray::Array1;

use crate::config::get_config;
use crate::utils::{
    self, check_and_adjust_instrument_name, get_central_wave, trapezoid, RsrData, RSR_DATA_VERSION,
    RSR_DATA_VERSION_FILENAME,
};

pub type RSRDict = HashMap<String, HashMap<String, RsrData>>;

pub struct RelativeSpectralResponse {
    pub platform_name: String,
    pub instrument: String,
    pub description: String,
    pub band_names: Vec<String>,
    pub rsr: RSRDict,
    pub unit: String,
    pub si_scale: f64,
    pub filename: PathBuf,
    wavespace: String,
}

impl RelativeSpectralResponse {
    pub fn new(
        platform_name: Option<&str>,
        instrument: Option<&str>,
        filename: Option<&Path>,
    ) -> Result<Self, String> {
        let (filepath, platform_name, instrument) =
            Self::sanitize_inputs(filename, platform_name, instrument)?;

        ensure_rsr_data(&filepath)?;

        let rsr_info = load_rsr_info_from_file(&filepath)?;

        let platform_name = match platform_name {
            Some(p) if !p.is_empty() => p.to_string(),
            _ => rsr_info.platform_name.clone(),
        };
        let instrument = match instrument {
            Some(i) if !i.is_empty() => i.to_string(),
            _ => rsr_info.instrument.clone(),
        };

        let mut rv = RelativeSpectralResponse {
            platform_name,
            instrument,
            description: rsr_info.description,
            band_names: rsr_info.band_names.clone(),
            rsr: rsr_info.rsr.clone(),
            unit: "1e-6 m".to_string(),
            si_scale: 1e-6,
            filename: filepath,
            wavespace: utils::WAVE_LENGTH.to_string(),
        };

        rv.band_names = rv.rsr.keys().cloned().collect();
        rv.band_names.sort();
        Ok(rv)
    }

    fn sanitize_inputs(
        filename: Option<&Path>,
        platform_name: Option<&str>,
        instrument: Option<&str>,
    ) -> Result<(PathBuf, Option<String>, Option<String>), String> {
        match filename {
            Some(f) => {
                if platform_name.is_some() || instrument.is_some() {
                    return Err(
                        "Either provide filename, or platform_name+instrument, not both"
                            .to_string(),
                    );
                }
                Ok((f.to_path_buf(), None, None))
            }
            None => {
                let pn =
                    platform_name.ok_or("platform_name is required when filename is not given")?;
                let instr =
                    instrument.ok_or("instrument is required when filename is not given")?;
                let adj_instr = check_and_adjust_instrument_name(pn, instr);

                let config = get_config(None);
                let fn_str = format!("rsr_{}_{}.h5", adj_instr, pn);
                let f = config.rsr_dir.join(&fn_str);
                debug!("Constructed RSR filename: {}", f.display());
                Ok((f, Some(pn.to_string()), Some(adj_instr)))
            }
        }
    }

    pub fn integral(&self, band_name: &str) -> HashMap<String, f64> {
        let mut result = HashMap::new();
        if let Some(detectors) = self.rsr.get(band_name) {
            for (det_name, det_data) in detectors {
                let integral = trapezoid(&det_data.response, &det_data.wavelength);
                result.insert(det_name.clone(), integral);
            }
        }
        result
    }

    pub fn convert(&mut self) {
        if self.wavespace == utils::WAVE_LENGTH {
            let mut new_rsr: RSRDict = HashMap::new();
            for (band_name, detectors) in &self.rsr {
                let mut new_detectors: HashMap<String, RsrData> = HashMap::new();
                for (det_name, det_data) in detectors {
                    let (wavenumber, response) =
                        utils::convert2wavenumber_rsr(&det_data.wavelength, &det_data.response);
                    let central_wn = get_central_wave(
                        &wavenumber,
                        &response,
                        &Array1::from_elem(wavenumber.len(), 1.0),
                    );
                    new_detectors.insert(
                        det_name.clone(),
                        RsrData {
                            wavelength: wavenumber,
                            response,
                            central_wavelength: central_wn,
                        },
                    );
                }
                new_rsr.insert(band_name.clone(), new_detectors);
            }
            self.rsr = new_rsr;
            self.wavespace = utils::WAVE_NUMBER.to_string();
            self.unit = "cm-1".to_string();
            self.si_scale = 100.0;
        } else {
            panic!("Conversion from wavenumber to wavelength not yet supported");
        }
    }

    pub fn get_bandname_from_wavelength(
        &self,
        wavelength: f64,
        epsilon: f64,
        multiple_bands: bool,
    ) -> Option<Vec<String>> {
        utils::get_bandname_from_wavelength(
            &self.instrument,
            wavelength,
            &self.rsr,
            epsilon,
            multiple_bands,
        )
    }

    pub fn resolve_band(&self, key: &str) -> Option<String> {
        if self.rsr.contains_key(key) {
            return Some(key.to_string());
        }
        let band_names = crate::bandnames::get_bandnames();
        if let Some(sensor_names) = band_names.get(self.instrument.as_str()) {
            if let Some(mapped) = sensor_names.get(key) {
                let mapped_str = mapped.to_string();
                if self.rsr.contains_key(&mapped_str) {
                    return Some(mapped_str);
                }
            }
        }
        if let Some(generic_names) = band_names.get("generic") {
            if let Some(mapped) = generic_names.get(key) {
                let mapped_str = mapped.to_string();
                if self.rsr.contains_key(&mapped_str) {
                    return Some(mapped_str);
                }
            }
        }
        None
    }
}

#[derive(Debug, Clone)]
pub struct RsrFileInfo {
    pub platform_name: String,
    pub instrument: String,
    pub description: String,
    pub band_names: Vec<String>,
    pub rsr: RSRDict,
}

pub fn load_rsr_info_from_file(filename: &Path) -> Result<RsrFileInfo, String> {
    if !filename.exists() {
        return Err(format!("RSR file not found: {}", filename.display()));
    }

    let h5file = hdf5_pure::File::open(filename)
        .map_err(|e| format!("Failed to open HDF5 file {}: {}", filename.display(), e))?;

    let root = h5file.root();
    let attrs = root.attrs().unwrap_or_default();

    let platform_name = get_platform_name(&attrs);
    let instrument = get_instrument(&attrs, &platform_name);
    let description = attr_to_string(&attrs, "description").unwrap_or_default();
    let band_names: Vec<String> = root.groups().unwrap_or_default();
    let rsr = get_relative_spectral_responses(&root, &band_names)?;

    Ok(RsrFileInfo {
        platform_name,
        instrument,
        description,
        band_names,
        rsr,
    })
}

fn get_platform_name(attrs: &HashMap<String, hdf5_pure::AttrValue>) -> String {
    let oscar_names: HashMap<&str, &str> = HashMap::from([
        ("eos-2", "EOS-Aqua"),
        ("eos-1", "EOS-Terra"),
        ("npp", "Suomi-NPP"),
        ("jpss-1", "NOAA-20"),
        ("jpss-2", "NOAA-21"),
        ("metop-a", "Metop-A"),
        ("metop-b", "Metop-B"),
        ("metop-c", "Metop-C"),
        ("meteosat-8", "Meteosat-8"),
        ("meteosat-9", "Meteosat-9"),
        ("meteosat-10", "Meteosat-10"),
        ("meteosat-11", "Meteosat-11"),
        ("himawari-8", "Himawari-8"),
        ("himawari-9", "Himawari-9"),
        ("goes-16", "GOES-16"),
        ("goes-17", "GOES-17"),
        ("goes-18", "GOES-18"),
        ("goes-19", "GOES-19"),
    ]);

    if let Some(name) = attr_to_string(attrs, "platform_name") {
        if let Some(mapped) = oscar_names.get(name.to_lowercase().as_str()) {
            return mapped.to_string();
        }
        return name;
    }

    let platform = attr_to_string(attrs, "platform").unwrap_or_default();
    let sat_number = attr_to_int(attrs, "sat_number");

    if !platform.is_empty() {
        if let Some(n) = sat_number {
            let combined = format!("{}-{}", platform, n).to_lowercase();
            if let Some(mapped) = oscar_names.get(combined.as_str()) {
                return mapped.to_string();
            }
            return format!("{}{}", platform, n);
        }
        return platform;
    }

    "unknown".to_string()
}

fn get_instrument(attrs: &HashMap<String, hdf5_pure::AttrValue>, platform_name: &str) -> String {
    if let Some(name) = attr_to_string(attrs, "sensor") {
        return name.to_lowercase();
    }

    let instruments = utils::get_instruments();
    if let Some(val) = instruments.get(platform_name) {
        match val {
            utils::InstrumentValue::Single(s) => return s.clone(),
            utils::InstrumentValue::List(l) => return l.first().cloned().unwrap_or_default(),
        }
    }

    "unknown".to_string()
}

fn get_relative_spectral_responses(
    root: &hdf5_pure::Group,
    band_names: &[String],
) -> Result<RSRDict, String> {
    let mut rsr: RSRDict = HashMap::new();

    for band_name in band_names {
        let band_group = match root.group(band_name) {
            Ok(g) => g,
            Err(_) => continue,
        };

        let mut detectors: HashMap<String, RsrData> = HashMap::new();
        let nd = get_number_of_detectors(&band_group);

        if nd > 1 {
            for i in 1..=nd {
                let det_name = format!("det-{}", i);
                if let Ok(det_group) = band_group.group(&det_name) {
                    if let Some(data) = read_detector_data(&det_group) {
                        detectors.insert(det_name, data);
                    }
                }
            }
        } else {
            if let Some(data) = read_detector_data(&band_group) {
                detectors.insert("det-1".to_string(), data);
            }
        }

        if !detectors.is_empty() {
            rsr.insert(band_name.clone(), detectors);
        }
    }

    Ok(rsr)
}

fn get_number_of_detectors(group: &hdf5_pure::Group) -> usize {
    attr_to_int(&group.attrs().unwrap_or_default(), "number_of_detectors").unwrap_or(1)
}

fn read_detector_data(group: &hdf5_pure::Group) -> Option<RsrData> {
    let dataset_name = if group.dataset("wavelength").is_ok() {
        "wavelength"
    } else if group.dataset("wavenumber").is_ok() {
        "wavenumber"
    } else {
        return None;
    };

    let wvl_ds = group.dataset(dataset_name).ok()?;
    let resp_ds = group.dataset("response").ok()?;

    let wavelength_raw: Vec<f64> = wvl_ds.read_f64().ok()?;
    let response: Vec<f64> = resp_ds.read_f64().ok()?;

    let band_attrs = group.attrs().unwrap_or_default();
    let scale: f64 = attr_to_f64(&band_attrs, "scale").unwrap_or(1e-6);

    let wavelength: Vec<f64> = wavelength_raw.iter().map(|&w| w * scale).collect();

    let wvl_arr = Array1::from_vec(wavelength.clone());
    let resp_arr = Array1::from_vec(response.clone());

    let central_wavelength = attr_to_f64(&band_attrs, "central_wavelength").unwrap_or_else(|| {
        get_central_wave(&wvl_arr, &resp_arr, &Array1::from_elem(wvl_arr.len(), 1.0))
    });

    Some(RsrData {
        wavelength: Array1::from_vec(wavelength),
        response: Array1::from_vec(response),
        central_wavelength,
    })
}

fn attr_to_string(attrs: &HashMap<String, hdf5_pure::AttrValue>, key: &str) -> Option<String> {
    attrs.get(key).and_then(|v| match v {
        hdf5_pure::AttrValue::String(s) => Some(s.clone()),
        hdf5_pure::AttrValue::AsciiString(s) => Some(s.clone()),
        _ => None,
    })
}

fn attr_to_int(attrs: &HashMap<String, hdf5_pure::AttrValue>, key: &str) -> Option<usize> {
    attrs.get(key).and_then(|v| match v {
        hdf5_pure::AttrValue::I32(v) => Some(*v as usize),
        hdf5_pure::AttrValue::I64(v) => Some(*v as usize),
        hdf5_pure::AttrValue::U32(v) => Some(*v as usize),
        hdf5_pure::AttrValue::U64(v) => Some(*v as usize),
        hdf5_pure::AttrValue::F64(v) => Some(*v as usize),
        _ => None,
    })
}

fn attr_to_f64(attrs: &HashMap<String, hdf5_pure::AttrValue>, key: &str) -> Option<f64> {
    attrs.get(key).and_then(|v| match v {
        hdf5_pure::AttrValue::F64(v) => Some(*v),
        hdf5_pure::AttrValue::I32(v) => Some(*v as f64),
        hdf5_pure::AttrValue::I64(v) => Some(*v as f64),
        hdf5_pure::AttrValue::U32(v) => Some(*v as f64),
        _ => None,
    })
}

pub fn check_and_download(dest_dir: Option<&Path>, dry_run: bool) {
    let config = get_config(None);
    let version_file = config.rsr_dir.join(RSR_DATA_VERSION_FILENAME);

    let current_version = if version_file.exists() {
        fs::read_to_string(&version_file).unwrap_or_default()
    } else {
        "v0.0.0".to_string()
    };

    if current_version.trim() != RSR_DATA_VERSION {
        info!(
            "RSR data version {} is outdated, downloading {}",
            current_version, RSR_DATA_VERSION
        );
        match crate::download::download_rsr(dest_dir, dry_run) {
            Ok(()) => info!("RSR data downloaded successfully"),
            Err(e) => warn!("Failed to download RSR data: {}", e),
        }
    }
}

fn ensure_rsr_data(filename: &Path) -> Result<(), String> {
    if !filename.exists() {
        let config = get_config(None);
        if config.download_from_internet {
            info!("RSR file not found, attempting download...");
            if let Err(e) = crate::download::download_rsr(None, false) {
                warn!("Failed to auto-download RSR data: {}", e);
            }
        }
        if !filename.exists() {
            return Err(format!(
                "RSR file not found at {} and auto-download failed. Run download_rsr first.",
                filename.display()
            ));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_requires_args() {
        let result = RelativeSpectralResponse::new(None, None, None);
        assert!(result.is_err());
    }

    #[test]
    fn test_sanitize_filename_only() {
        let (_, pn, instr) =
            RelativeSpectralResponse::sanitize_inputs(Some(Path::new("/test/rsr.h5")), None, None)
                .unwrap();
        assert!(pn.is_none());
        assert!(instr.is_none());
    }

    #[test]
    fn test_sanitize_requires_either() {
        let result = RelativeSpectralResponse::sanitize_inputs(None, None, Some("modis"));
        assert!(result.is_err());
    }
}
