use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

use log::debug;
use serde_yaml::Value;

#[allow(dead_code)]
const BUILTIN_RSR_DIR: &str = "~/.local/share/pyspectral/rsr";
#[allow(dead_code)]
const BUILTIN_RAYLEIGH_DIR: &str = "~/.local/share/pyspectral/rayleigh";

#[derive(Debug, Clone)]
pub struct Config {
    pub rsr_dir: PathBuf,
    pub rayleigh_dir: PathBuf,
    pub download_from_internet: bool,
    pub raw: HashMap<String, Value>,
}

impl Default for Config {
    fn default() -> Self {
        let data_dir = default_data_dir();
        Config {
            rsr_dir: data_dir.clone(),
            rayleigh_dir: data_dir,
            download_from_internet: true,
            raw: HashMap::new(),
        }
    }
}

fn default_data_dir() -> PathBuf {
    if let Some(data_dir) = dirs::data_local_dir() {
        data_dir.join("pyspectral")
    } else {
        PathBuf::from("/tmp/pyspectral")
    }
}

fn expand_home(path: &str) -> PathBuf {
    if path.starts_with("~/") {
        if let Some(home) = dirs::home_dir() {
            return home.join(&path[2..]);
        }
    }
    PathBuf::from(path)
}

pub fn get_config(config_file: Option<&Path>) -> Config {
    let cf = match config_file {
        Some(p) => p.to_path_buf(),
        None => {
            if let Ok(env_path) = std::env::var("PSP_CONFIG_FILE") {
                let p = PathBuf::from(&env_path);
                if !p.is_file() {
                    panic!("{} pointed to by PSP_CONFIG_FILE does not exist", env_path);
                }
                p
            } else {
                return Config::default();
            }
        }
    };

    let content = fs::read_to_string(&cf)
        .unwrap_or_else(|e| panic!("Failed to read config file {}: {}", cf.display(), e));

    let raw: Value = serde_yaml::from_str(&content)
        .unwrap_or_else(|e| panic!("Failed to parse config file: {}", e));

    let mapping = raw.as_mapping().cloned().unwrap_or_default();

    let download_from_internet = mapping
        .get(&Value::String("download_from_internet".into()))
        .and_then(|v| v.as_bool())
        .unwrap_or(true);

    let rsr_dir = mapping
        .get(&Value::String("rsr_dir".into()))
        .and_then(|v| v.as_str())
        .map(expand_home)
        .unwrap_or_else(default_data_dir);

    let rayleigh_dir = mapping
        .get(&Value::String("rayleigh_dir".into()))
        .and_then(|v| v.as_str())
        .map(expand_home)
        .unwrap_or_else(default_data_dir);

    fs::create_dir_all(&rsr_dir).ok();
    fs::create_dir_all(&rayleigh_dir).ok();

    let raw_map: HashMap<String, Value> = mapping
        .iter()
        .map(|(k, v)| {
            let key = k.as_str().unwrap_or("").to_string();
            (key, v.clone())
        })
        .collect();

    debug!(
        "Config loaded: rsr_dir={:?}, rayleigh_dir={:?}",
        rsr_dir, rayleigh_dir
    );

    Config {
        rsr_dir,
        rayleigh_dir,
        download_from_internet,
        raw: raw_map,
    }
}

pub fn recursive_dict_update(base: Value, update: &Value) -> Value {
    match (base, update) {
        (Value::Mapping(mut b), Value::Mapping(u)) => {
            for (k, v) in u {
                let entry = b.entry(k.clone());
                match entry {
                    serde_yaml::mapping::Entry::Occupied(mut o) => {
                        let updated = recursive_dict_update(o.get().clone(), v);
                        o.insert(updated);
                    }
                    serde_yaml::mapping::Entry::Vacant(va) => {
                        va.insert(v.clone());
                    }
                }
            }
            Value::Mapping(b)
        }
        (_, u) => u.clone(),
    }
}
