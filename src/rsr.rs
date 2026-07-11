use ndarray::Array1;

use crate::utils::{convert2wavenumber_rsr, get_central_wave};

pub struct Rsr {
    pub wavelength: Array1<f64>,
    pub response: Array1<f64>,
    pub central_wavelength: f64,
}

pub struct DetectorRsr {
    pub detectors: Vec<PerDetectorRsr>,
}

pub struct PerDetectorRsr {
    pub name: String,
    pub wavelength: Array1<f64>,
    pub response: Array1<f64>,
    pub central_wavelength: f64,
    pub wavenumber: Option<Array1<f64>>,
}

impl PerDetectorRsr {
    pub fn new(name: String, wavelength: Array1<f64>, response: Array1<f64>) -> Self {
        let central_wavelength = get_central_wave(&wavelength, &response, 1.0);
        PerDetectorRsr {
            name,
            wavelength,
            response,
            central_wavelength,
            wavenumber: None,
        }
    }
}

impl Rsr {
    pub fn new(wavelength: Array1<f64>, response: Array1<f64>) -> Self {
        let central_wavelength = get_central_wave(&wavelength, &response, 1.0);
        Rsr {
            wavelength,
            response,
            central_wavelength,
        }
    }

    pub fn from_detectors(detectors: Vec<PerDetectorRsr>) -> DetectorRsr {
        DetectorRsr { detectors }
    }

    pub fn to_wavenumber(&self) -> (Array1<f64>, Array1<f64>) {
        convert2wavenumber_rsr(&self.wavelength, &self.response)
    }
}
