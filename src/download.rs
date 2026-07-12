use std::fs;
use std::io::{self, Read};
use std::path::{Path, PathBuf};

use log::{debug, info};
use ureq::Agent;

use crate::config::{get_config, Config};
use crate::utils::{
    AEROSOL_TYPES, ATM_CORRECTION_LUT_VERSION, HTTPS_RAYLEIGH_LUTS, HTTP_PYSPECTRAL_RSR,
    RSR_DATA_VERSION, RSR_DATA_VERSION_FILENAME,
};

pub fn download_rsr(dest_dir: Option<&Path>, dry_run: bool) -> io::Result<()> {
    let config = get_config(None);
    let dest = dest_dir.unwrap_or(&config.rsr_dir);
    let dest = dest.to_path_buf();

    fs::create_dir_all(&dest)?;

    let tarball_path = dest.join("pyspectral_rsr_data.tgz");

    if dry_run {
        info!(
            "Dry run: would download {} to {}",
            HTTP_PYSPECTRAL_RSR,
            tarball_path.display()
        );
        return Ok(());
    }

    if !config.download_from_internet {
        info!("Download from internet disabled in config");
        return Ok(());
    }

    info!("Downloading RSR data from {}", HTTP_PYSPECTRAL_RSR);
    download_file(HTTP_PYSPECTRAL_RSR, &tarball_path)?;

    extract_tarball(&tarball_path, &dest)?;

    let version_file = dest.join(RSR_DATA_VERSION_FILENAME);
    fs::write(&version_file, RSR_DATA_VERSION)?;

    info!("RSR data downloaded to {}", dest.display());
    Ok(())
}

pub fn download_luts(aerosol_types: Option<&[String]>, dry_run: bool) -> io::Result<()> {
    let config = get_config(None);
    let types: Vec<String> = aerosol_types
        .map(|v| v.to_vec())
        .unwrap_or_else(|| AEROSOL_TYPES.iter().map(|s| s.to_string()).collect());

    for aerosol_type in &types {
        let url = match HTTPS_RAYLEIGH_LUTS.get(aerosol_type.as_str()) {
            Some(u) => u,
            None => {
                debug!("No LUT URL for aerosol type: {}", aerosol_type);
                continue;
            }
        };

        let lut_dir = get_rayleigh_lut_dir(&config, aerosol_type);
        if !dry_run {
            fs::create_dir_all(&lut_dir)?;
        }

        let tarball_path = lut_dir.join("pyspectral_rayleigh_correction_luts.tgz");

        if dry_run {
            info!(
                "Dry run: would download {} to {}",
                url,
                tarball_path.display()
            );
            continue;
        }

        info!("Downloading LUT data for {} from {}", aerosol_type, url);
        download_file(url, &tarball_path)?;

        extract_tarball(&tarball_path, &lut_dir)?;
        fs::remove_file(&tarball_path)?;

        if let Some(ver_info) = ATM_CORRECTION_LUT_VERSION.get(aerosol_type.as_str()) {
            let version_file = lut_dir.join(ver_info.filename);
            fs::write(&version_file, ver_info.version)?;
        }
    }

    Ok(())
}

pub fn get_rayleigh_lut_dir(config: &Config, aerosol_type: &str) -> PathBuf {
    config.rayleigh_dir.join(aerosol_type)
}

fn download_file(url: &str, dest: &Path) -> io::Result<()> {
    let agent = Agent::new_with_defaults();
    let resp = agent
        .get(url)
        .call()
        .map_err(|e| io::Error::other(format!("HTTP request failed: {}", e)))?;

    let mut body = Vec::new();
    resp.into_body().as_reader().read_to_end(&mut body)?;

    fs::write(dest, &body)?;
    debug!("Downloaded {} to {}", url, dest.display());
    Ok(())
}

fn extract_tarball(tarball_path: &Path, dest_dir: &Path) -> io::Result<()> {
    let file = fs::File::open(tarball_path)?;
    let decoder = flate2::read::GzDecoder::new(file);
    let mut archive = tar::Archive::new(decoder);
    archive.unpack(dest_dir)?;
    debug!(
        "Extracted {} to {}",
        tarball_path.display(),
        dest_dir.display()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn test_download_rsr_dry_run() {
        let dir = tempdir().unwrap();
        let result = download_rsr(Some(dir.path()), true);
        assert!(result.is_ok());
    }

    #[test]
    fn test_download_luts_dry_run() {
        let result = download_luts(None, true);
        assert!(result.is_ok());
    }

    #[test]
    fn test_download_luts_specific_type_dry_run() {
        let types: Vec<String> = vec!["marine_clean_aerosol".to_string()];
        let result = download_luts(Some(&types), true);
        assert!(result.is_ok());
    }

    #[test]
    fn test_lut_urls_are_valid_http() {
        for url in HTTPS_RAYLEIGH_LUTS.values() {
            assert!(
                url.starts_with("https://"),
                "URL missing https:// scheme: {}",
                url
            );
            assert!(
                url.contains("zenodo.org/records/"),
                "URL missing zenodo.org host: {}",
                url
            );
            assert!(url.ends_with(".tgz"), "URL should end with .tgz: {}", url);
        }
    }

    #[test]
    fn test_marine_clean_aerosol_url() {
        let url = HTTPS_RAYLEIGH_LUTS
            .get("marine_clean_aerosol")
            .expect("marine_clean_aerosol should exist");
        assert_eq!(
            url.as_str(),
            "https://zenodo.org/records/1288441/files/pyspectral_atm_correction_luts_marine_clean_aerosol.tgz"
        );
    }

    #[test]
    fn test_all_aerosol_types_have_urls() {
        for aerosol_type in AEROSOL_TYPES {
            assert!(
                HTTPS_RAYLEIGH_LUTS.contains_key(*aerosol_type),
                "AEROSOL_TYPES contains '{}' but HTTPS_RAYLEIGH_LUTS does not",
                aerosol_type
            );
        }
    }

    #[test]
    fn test_get_rayleigh_lut_dir() {
        let config = Config {
            rayleigh_dir: PathBuf::from("/test/rayleigh"),
            ..Default::default()
        };
        let dir = get_rayleigh_lut_dir(&config, "marine_clean_aerosol");
        assert_eq!(dir, PathBuf::from("/test/rayleigh/marine_clean_aerosol"));
    }
}
