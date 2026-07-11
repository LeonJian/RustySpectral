use std::collections::HashMap;
use std::path::PathBuf;

use crate::config::get_config;
use crate::utils::check_and_adjust_instrument_name;
use log::{debug, warn};

pub struct InstrumentRSR {
    pub platform_name: String,
    pub instrument: String,
    pub bandname: String,
    pub bandnames: Vec<String>,
    pub filenames: HashMap<String, Option<String>>,
    pub rsr: Option<()>,
    pub output_dir: PathBuf,
    pub path: Option<PathBuf>,
    pub filename: Option<String>,
}

impl InstrumentRSR {
    pub fn new(bandname: &str, platform_name: &str, bandnames: &[String]) -> Self {
        let instrument = check_and_adjust_instrument_name(platform_name, "");
        let mut filenames = HashMap::new();
        for band in bandnames {
            filenames.insert(band.clone(), None);
        }

        InstrumentRSR {
            platform_name: platform_name.to_string(),
            instrument,
            bandname: bandname.to_string(),
            bandnames: bandnames.to_vec(),
            filenames,
            rsr: None,
            output_dir: PathBuf::from("./"),
            path: None,
            filename: None,
        }
    }

    pub fn get_options_from_config(&mut self) {
        let config = get_config(None);
        self.output_dir = config.rsr_dir.clone();

        let key = format!("{}-{}", self.platform_name, self.instrument);
        if let Some(val) = config.raw.get(&key) {
            if let Some(path_str) = val.get("path").and_then(|v| v.as_str()) {
                self.path = Some(PathBuf::from(path_str));
            }
            if let Some(fn_str) = val.get("filename").and_then(|v| v.as_str()) {
                self.filename = Some(fn_str.to_string());
            }
        }

        debug!(
            "RSR config: output_dir={:?}, path={:?}",
            self.output_dir, self.path
        );
    }

    pub fn get_bandfilenames(&mut self) {
        let config = get_config(None);

        if let Some(path) = &self.path {
            for band in &self.bandnames.clone() {
                let lookup_key = format!("{}-{}", self.platform_name, self.instrument);
                if let Some(platform_cfg) = config.raw.get(&lookup_key) {
                    if let Some(band_cfg) = platform_cfg.get(band) {
                        if let Some(fn_str) = band_cfg.as_str() {
                            let full_path = path.join(fn_str);
                            if full_path.exists() {
                                self.filenames.insert(
                                    band.clone(),
                                    Some(full_path.to_string_lossy().to_string()),
                                );
                            } else {
                                warn!(
                                    "Couldn't find an existing file for this band: {}",
                                    full_path.display()
                                );
                                self.filenames.insert(band.clone(), None);
                            }
                        }
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_instrument_rsr() {
        let bandnames = vec!["VIS006".to_string(), "IR_108".to_string()];
        let rsr = InstrumentRSR::new("VIS006", "Meteosat-10", &bandnames);
        assert_eq!(rsr.bandname, "VIS006");
        assert_eq!(rsr.bandnames.len(), 2);
    }

    #[test]
    fn test_get_options_from_config() {
        let bandnames = vec!["VIS006".to_string()];
        let mut rsr = InstrumentRSR::new("VIS006", "Meteosat-10", &bandnames);
        rsr.get_options_from_config();
        assert!(rsr.output_dir.to_string_lossy().len() > 0);
    }
}
