use ndarray::Array1;
use std::collections::HashMap;

use crate::blackbody::{blackbody, blackbody_rad2temp, C_SPEED, H_PLANCK, K_BOLTZMANN};
use crate::utils::trapezoid;

pub const TB_MIN: f64 = 150.0;
pub const TB_MAX: f64 = 360.0;
pub const EPSILON: f64 = 0.01;

pub type SeviriParams = (f64, f64, f64);

pub fn get_seviri_params() -> HashMap<&'static str, HashMap<&'static str, SeviriParams>> {
    let mut m = HashMap::new();
    macro_rules! b { ($band:expr, $($plat:expr => ($vc:expr, $alpha:expr, $beta:expr)),+) => {
        let mut p = HashMap::new();
        $( p.insert($plat, ($vc, $alpha, $beta)); )+
        m.insert($band, p);
    }; }
    b!("IR3.9",
        "Meteosat-8" => (2567.330, 0.9956, 3.410),
        "Meteosat-9" => (2568.832, 0.9954, 3.438)
    );
    b!("WV6.2",
        "Meteosat-8" => (1598.103, 0.9962, 2.218),
        "Meteosat-9" => (1600.548, 0.9963, 2.185)
    );
    b!("WV7.3",
        "Meteosat-8" => (1362.081, 0.9991, 0.478),
        "Meteosat-9" => (1360.330, 0.9991, 0.470)
    );
    b!("IR8.7",
        "Meteosat-8" => (1149.069, 0.9996, 0.179),
        "Meteosat-9" => (1148.620, 0.9996, 0.179)
    );
    b!("IR9.7",
        "Meteosat-8" => (1034.343, 0.9999, 0.060),
        "Meteosat-9" => (1035.289, 0.9999, 0.056)
    );
    b!("IR10.8",
        "Meteosat-8" => (930.647, 0.9983, 0.625),
        "Meteosat-9" => (931.700, 0.9983, 0.640)
    );
    b!("IR12.0",
        "Meteosat-8" => (839.660, 0.9988, 0.397),
        "Meteosat-9" => (836.445, 0.9988, 0.408)
    );
    b!("IR13.4",
        "Meteosat-8" => (752.387, 0.9981, 0.578),
        "Meteosat-9" => (751.792, 0.9981, 0.561)
    );
    m
}

pub static SEVIRI: once_cell::sync::Lazy<
    HashMap<&'static str, HashMap<&'static str, SeviriParams>>,
> = once_cell::sync::Lazy::new(get_seviri_params);

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

pub struct RadTbConverter {
    pub platform_name: String,
    pub instrument: String,
    pub band: String,
    pub wavelength: Array1<f64>,
    pub response: Array1<f64>,
    pub central_wavelength: f64,
    pub rsr_integral: f64,
    pub detector: String,
}

impl RadTbConverter {
    pub fn new(
        platform_name: &str,
        instrument: &str,
        band: &str,
        wavelength: Array1<f64>,
        response: Array1<f64>,
    ) -> Self {
        let central_wavelength = crate::utils::get_central_wave(
            &wavelength,
            &response,
            &Array1::from_elem(wavelength.len(), 1.0),
        );
        let rsr_integral = trapezoid(&response, &wavelength);
        RadTbConverter {
            platform_name: platform_name.to_string(),
            instrument: instrument.to_string(),
            band: band.to_string(),
            wavelength,
            response,
            central_wavelength,
            rsr_integral,
            detector: "det-1".to_string(),
        }
    }

    pub fn with_detector(mut self, detector: &str) -> Self {
        self.detector = detector.to_string();
        self
    }

    pub fn tb2radiance(&self, tb: &Array1<f64>, normalized: bool) -> Array1<f64> {
        tb.mapv(|t| {
            let planck_vals: Array1<f64> = self.wavelength.mapv(|w| blackbody(w, t));
            let product = &planck_vals * &self.response;
            let integrated = trapezoid(&product, &self.wavelength);
            if normalized && self.rsr_integral > 0.0 {
                integrated / self.rsr_integral
            } else {
                integrated
            }
        })
    }

    pub fn radiance2tb(&self, rad: &Array1<f64>) -> Array1<f64> {
        rad.mapv(|r| radiance2tb(r, self.central_wavelength * 1e-6))
    }

    pub fn make_tb2rad_lut(
        &self,
        tb_resolution: f64,
        normalized: bool,
    ) -> (Array1<f64>, Array1<f64>) {
        let n = ((TB_MAX - TB_MIN) / tb_resolution).round() as usize + 1;
        let mut lut_tb = Array1::zeros(n);
        let mut lut_rad = Array1::zeros(n);
        for i in 0..n {
            let tb = TB_MIN + i as f64 * tb_resolution;
            lut_tb[i] = tb;
            lut_rad[i] = self.tb2radiance(&Array1::from_elem(1, tb), normalized)[0];
        }
        (lut_tb, lut_rad)
    }
}

pub struct SeviriRadTbConverter {
    pub platform_name: String,
    pub band: String,
    pub vc: f64,
    pub alpha: f64,
    pub beta: f64,
}

impl SeviriRadTbConverter {
    pub fn new(platform_name: &str, band: &str) -> Option<Self> {
        let seviri = SEVIRI.get(band)?;
        let params = seviri.get(platform_name)?;
        Some(SeviriRadTbConverter {
            platform_name: platform_name.to_string(),
            band: band.to_string(),
            vc: params.0 * 100.0,
            alpha: params.1,
            beta: params.2,
        })
    }

    pub fn radiance2tb(&self, rad: f64) -> f64 {
        seviri_radiance2tb(rad, self.vc, self.alpha, self.beta)
    }

    pub fn tb2radiance(&self, tb: f64) -> f64 {
        seviri_tb2radiance(tb, self.vc, self.alpha, self.beta)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::{assert_abs_diff_eq, assert_relative_eq};
    use ndarray::array;

    const TEST_WVL: &[f64] = &[
        3.6123999, 3.6163599, 3.6264927, 3.6363862, 3.646468, 3.6564937, 3.6664478, 3.6765388,
        3.6865413, 3.6964585, 3.7065142, 3.716509, 3.7264658, 3.7364102, 3.7463682, 3.7563652,
        3.7664226, 3.7763396, 3.7863384, 3.7964207, 3.8063589, 3.8163606, 3.8264089, 3.8364836,
        3.8463381, 3.8563975, 3.8664163, 3.8763755, 3.8864797, 3.8964978, 3.9064275, 3.9164873,
        3.9264729, 3.9364026, 3.9465107, 3.9535347,
    ];

    const TEST_RESP: &[f64] = &[
        0.01, 0.0118, 0.01987, 0.03226, 0.05028, 0.0849, 0.16645, 0.33792, 0.59106, 0.81815,
        0.96077, 0.92855, 0.86008, 0.8661, 0.87697, 0.85412, 0.88922, 0.9541, 0.95687, 0.91037,
        0.91058, 0.94256, 0.94719, 0.94808, 1.0, 0.92676, 0.67429, 0.44715, 0.27762, 0.14852,
        0.07141, 0.04151, 0.02925, 0.02085, 0.01414, 0.01,
    ];

    fn test_wavelength() -> ndarray::Array1<f64> {
        ndarray::Array1::from_vec(TEST_WVL.iter().map(|v| v * 1e-6).collect())
    }

    fn test_response() -> ndarray::Array1<f64> {
        ndarray::Array1::from_vec(TEST_RESP.to_vec())
    }

    #[test]
    fn test_radiance2tb_simple() {
        let wvl = test_wavelength();
        let resp = test_response();
        let central_wvl = 3.780_281_935e-6;
        let rad = tb2radiance_normalized(300.0, &wvl, &resp);
        let tb = radiance2tb(rad, central_wvl);
        assert_abs_diff_eq!(tb, 300.0, epsilon = 0.5);
    }

    #[test]
    fn test_tb2radiance_single() {
        let wvl = test_wavelength();
        let resp = test_response();
        let rad = tb2radiance_simple(300.0, &wvl, &resp);
        assert!(rad > 0.0);
        let rad2 = tb2radiance_simple(350.0, &wvl, &resp);
        assert!(rad2 > rad);
    }

    #[test]
    fn test_tb2radiance_array() {
        let wvl = test_wavelength();
        let resp = test_response();
        let tbs = array![200.0, 270.0, 300.0, 302.0, 350.0_f64];
        let rads = tb2radiance_array(&tbs, &wvl, &resp);
        assert_eq!(rads.len(), 5);
        assert!(rads[0] > 0.0);
        assert!(rads[0] < rads[1]);
        assert!(rads[3] < rads[4]);
    }

    #[test]
    fn test_tb2radiance_normalized() {
        let wvl = test_wavelength();
        let resp = test_response();
        let rad = tb2radiance_normalized(300.0, &wvl, &resp);
        assert!(rad > 0.0);
    }

    #[test]
    fn test_make_tb2rad_lut() {
        let wvl = test_wavelength();
        let resp = test_response();
        let (lut_tb, lut_rad) = make_tb2rad_lut(&wvl, &resp, 0.1);
        assert!(!lut_tb.is_empty());
        assert_eq!(lut_tb.len(), lut_rad.len());
        for i in 1..lut_rad.len() {
            assert!(lut_rad[i] > lut_rad[i - 1]);
        }
        assert!(lut_tb[0] >= 150.0);
        assert!(lut_tb[lut_tb.len() - 1] <= 360.0);
    }

    #[test]
    fn test_seviri_tb2radiance() {
        let vc = 2568.832 * 100.0;
        let alpha = 0.9954;
        let beta = 3.438;
        let rad = seviri_tb2radiance(300.0, vc, alpha, beta);
        assert!(rad > 0.0);
        assert_relative_eq!(rad, 9.797_091e-6, epsilon = 1e-8);
    }

    #[test]
    fn test_seviri_radiance2tb() {
        let vc = 2568.832 * 100.0;
        let alpha = 0.9954;
        let beta = 3.438;
        let rad = seviri_tb2radiance(300.0, vc, alpha, beta);
        let tb = seviri_radiance2tb(rad, vc, alpha, beta);
        assert_relative_eq!(tb, 300.0, epsilon = 1e-4);
    }

    #[test]
    fn test_seviri_multiple_temperatures() {
        let vc = 2568.832 * 100.0;
        let alpha = 0.9954;
        let beta = 3.438;
        for tb in [200.0, 250.0, 300.0, 350.0].iter() {
            let rad = seviri_tb2radiance(*tb, vc, alpha, beta);
            let tb_back = seviri_radiance2tb(rad, vc, alpha, beta);
            assert_relative_eq!(tb_back, *tb, epsilon = 1e-4);
        }
    }
}
