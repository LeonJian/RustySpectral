use ndarray::Array2;

pub const H_PLANCK: f64 = 6.626_069_57e-34;
pub const K_BOLTZMANN: f64 = 1.380_648_8e-23;
pub const C_SPEED: f64 = 2.997_924_58e8;

const PLANCK_C1: f64 = H_PLANCK * C_SPEED / K_BOLTZMANN;
const PLANCK_C2: f64 = 2.0 * H_PLANCK * C_SPEED * C_SPEED;
const EPSILON: f64 = 0.000_001;

pub fn blackbody(wave: f64, temperature: f64) -> f64 {
    planck(wave, temperature)
}

pub fn blackbody_wn(wavenumber: f64, temperature: f64) -> f64 {
    planck_wn(wavenumber, temperature)
}

pub fn planck(wave: f64, temperature: f64) -> f64 {
    if temperature.abs() < EPSILON {
        return f64::NAN;
    }
    let nom = PLANCK_C2 / wave.powi(5);
    let arg1 = PLANCK_C1 / wave;
    let exp_arg = arg1 / temperature;
    if exp_arg.is_infinite() || exp_arg < 0.0 {
        return f64::NAN;
    }
    nom / (exp_arg.exp() - 1.0)
}

pub fn planck_wn(wavenumber: f64, temperature: f64) -> f64 {
    if temperature.abs() < EPSILON {
        return f64::NAN;
    }
    let nom = PLANCK_C2 * wavenumber.powi(3);
    let arg1 = PLANCK_C1 * wavenumber;
    let exp_arg = arg1 / temperature;
    if exp_arg.is_infinite() || exp_arg < 0.0 {
        return f64::NAN;
    }
    nom / (exp_arg.exp() - 1.0)
}

pub fn planck_array_wavelength(wave: f64, temperature: &Array2<f64>) -> Array2<f64> {
    temperature.mapv(|t| planck(wave, t))
}

pub fn planck_array_wn(wavenumber: f64, temperature: &Array2<f64>) -> Array2<f64> {
    temperature.mapv(|t| planck_wn(wavenumber, t))
}

pub fn blackbody_rad2temp(wavelength: f64, radiance: f64) -> f64 {
    if radiance <= 0.0 {
        return f64::NAN;
    }
    let arg = PLANCK_C2 / (radiance * wavelength.powi(5)) + 1.0;
    if arg <= 0.0 {
        return f64::NAN;
    }
    PLANCK_C1 / (wavelength * arg.ln())
}

pub fn blackbody_wn_rad2temp(wavenumber: f64, radiance: f64) -> f64 {
    if radiance <= 0.0 {
        return f64::NAN;
    }
    let arg = (PLANCK_C2 * wavenumber.powi(3)) / radiance + 1.0;
    if arg <= 0.0 {
        return f64::NAN;
    }
    PLANCK_C1 * wavenumber / arg.ln()
}

pub fn blackbody_rad2temp_array(wavelength: f64, radiance: &Array2<f64>) -> Array2<f64> {
    radiance.mapv(|r| blackbody_rad2temp(wavelength, r))
}

pub fn blackbody_wn_rad2temp_array(wavenumber: f64, radiance: &Array2<f64>) -> Array2<f64> {
    radiance.mapv(|r| blackbody_wn_rad2temp(wavenumber, r))
}
