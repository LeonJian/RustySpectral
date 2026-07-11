use ndarray::Array1;

use crate::blackbody::{blackbody, blackbody_rad2temp, C_SPEED, H_PLANCK, K_BOLTZMANN};
use crate::utils::trapezoid;

pub fn radiance2tb(radiance: f64, wavelength: f64) -> f64 {
    blackbody_rad2temp(wavelength, radiance)
}

pub fn tb2radiance_simple(tb: f64, wavelength: &Array1<f64>, response: &Array1<f64>) -> f64 {
    let n = wavelength.len();
    let mut planck_vals = Array1::zeros(n);
    for i in 0..n {
        planck_vals[i] = blackbody(wavelength[i], tb);
    }
    let product = &planck_vals * response;
    trapezoid(&product, wavelength)
}

pub fn tb2radiance_array(
    tb: &Array1<f64>,
    wavelength: &Array1<f64>,
    response: &Array1<f64>,
) -> Array1<f64> {
    tb.mapv(|t| tb2radiance_simple(t, wavelength, response))
}

pub fn tb2radiance_normalized(tb: f64, wavelength: &Array1<f64>, response: &Array1<f64>) -> f64 {
    let integrated = tb2radiance_simple(tb, wavelength, response);
    let rsr_integral = trapezoid(response, wavelength);
    if rsr_integral == 0.0 {
        return 0.0;
    }
    integrated / rsr_integral
}

pub fn make_tb2rad_lut(
    wavelength: &Array1<f64>,
    response: &Array1<f64>,
    tb_resolution: f64,
) -> (Array1<f64>, Array1<f64>) {
    let tb_min = 150.0;
    let tb_max = 360.0;
    let n = ((tb_max - tb_min) / tb_resolution).round() as usize + 1;
    let mut lut_tb = Array1::zeros(n);
    let mut lut_rad = Array1::zeros(n);

    for i in 0..n {
        let tb = tb_min + i as f64 * tb_resolution;
        lut_tb[i] = tb;
        lut_rad[i] = tb2radiance_normalized(tb, wavelength, response);
    }

    (lut_tb, lut_rad)
}

pub fn seviri_radiance2tb(radiance: f64, central_wavenumber: f64, alpha: f64, beta: f64) -> f64 {
    let c1 = 2.0 * H_PLANCK * C_SPEED * C_SPEED;
    let c2 = H_PLANCK * C_SPEED / K_BOLTZMANN;

    let vc = central_wavenumber;
    let arg = c1 * vc.powi(3) / radiance + 1.0;

    c2 * vc / (alpha * arg.ln()) - beta / alpha
}

pub fn seviri_tb2radiance(tb: f64, central_wavenumber: f64, alpha: f64, beta: f64) -> f64 {
    let c1 = 2.0 * H_PLANCK * C_SPEED * C_SPEED;
    let c2 = H_PLANCK * C_SPEED / K_BOLTZMANN;

    let vc = central_wavenumber;
    c1 * vc.powi(3) / ((c2 * vc / (alpha * tb + beta)).exp() - 1.0)
}
