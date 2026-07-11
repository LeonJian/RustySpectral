use ndarray::Array1;
use std::f64::consts::PI;

use crate::blackbody::blackbody;
use crate::utils::trapezoid;

pub const TERMINATOR_LIMIT: f64 = 85.0;
const EPSILON: f64 = 0.005;

pub struct ReflectanceCalculator {
    pub platform_name: String,
    pub instrument: String,
    pub solar_flux: f64,
    pub central_wavelength: f64,
    wavelength: Array1<f64>,
    response: Array1<f64>,
    rsr_integral: f64,
    pub sunz_threshold: f64,
    pub masking_limit: Option<f64>,
}

impl ReflectanceCalculator {
    pub fn new(platform_name: &str, instrument: &str) -> Self {
        ReflectanceCalculator {
            platform_name: platform_name.to_string(),
            instrument: instrument.to_string(),
            solar_flux: 0.0,
            central_wavelength: 3.78e-6,
            wavelength: Array1::zeros(0),
            response: Array1::zeros(0),
            rsr_integral: 0.0,
            sunz_threshold: TERMINATOR_LIMIT,
            masking_limit: Some(TERMINATOR_LIMIT),
        }
    }

    pub fn with_rsr(mut self, wavelength: Array1<f64>, response: Array1<f64>) -> Self {
        self.rsr_integral = trapezoid(&response, &wavelength);
        self.central_wavelength = crate::utils::get_central_wave(&wavelength, &response, 1.0);
        self.wavelength = wavelength;
        self.response = response;
        self
    }

    pub fn with_solar_flux(mut self, solar_flux: f64) -> Self {
        self.solar_flux = solar_flux;
        self
    }

    pub fn with_central_wavelength(mut self, cw: f64) -> Self {
        self.central_wavelength = cw;
        self
    }

    pub fn with_sunz_threshold(mut self, threshold: f64) -> Self {
        self.sunz_threshold = threshold;
        self
    }

    pub fn with_masking_limit(mut self, limit: Option<f64>) -> Self {
        self.masking_limit = limit;
        self
    }

    pub fn solar_radiance(&self, sun_zenith: &Array1<f64>) -> Array1<f64> {
        let sunz = sun_zenith.mapv(|sz| sz.clamp(0.0, self.sunz_threshold));
        let mu0 = sunz.mapv(|sz| sz.to_radians().cos());
        mu0.mapv(|mu| self.solar_flux * mu / PI)
    }

    pub fn tb2radiance(&self, tb: &Array1<f64>) -> Array1<f64> {
        if self.wavelength.len() > 0 && self.rsr_integral > 0.0 {
            tb.mapv(|t| {
                let planck_vals: Array1<f64> = self.wavelength.mapv(|w| blackbody(w, t));
                let product = &planck_vals * &self.response;
                trapezoid(&product, &self.wavelength) / self.rsr_integral
            })
        } else {
            tb.mapv(|t| blackbody(self.central_wavelength, t))
        }
    }

    pub fn derive_rad39_corr(&self, bt11: &Array1<f64>, bt13: &Array1<f64>) -> Array1<f64> {
        let n = bt11.len();
        let mut corr = Array1::ones(n);
        for i in 0..n {
            let b11 = bt11[i];
            let b13 = bt13[i];
            if b11 > 0.0 {
                corr[i] = (b11 - 0.25 * (b11 - b13)).powi(4) / b11.powi(4);
            }
        }
        corr
    }

    pub fn reflectance_from_tbs(
        &self,
        sun_zenith: &Array1<f64>,
        tb_near_ir: &Array1<f64>,
        tb_thermal: &Array1<f64>,
        tb_ir_co2: Option<&Array1<f64>>,
    ) -> Array1<f64> {
        let n = tb_near_ir.len();
        assert_eq!(sun_zenith.len(), n);
        assert_eq!(tb_thermal.len(), n);

        let rad3x_t11 = self.tb2radiance(tb_thermal);
        let rad3x = self.tb2radiance(tb_near_ir);

        let solar_rad = self.solar_radiance(sun_zenith);

        let rad3x_correction = if let Some(co2_tb) = tb_ir_co2 {
            self.derive_rad39_corr(tb_thermal, co2_tb)
        } else {
            Array1::ones(n)
        };

        let corrected_thermal = if self.rsr_integral > 0.0 {
            &rad3x_t11 * self.rsr_integral * &rad3x_correction
        } else {
            &rad3x_t11 * &rad3x_correction
        };

        let l_nir = if self.rsr_integral > 0.0 {
            &rad3x * self.rsr_integral
        } else {
            rad3x.clone()
        };

        let nomin = &l_nir - &corrected_thermal;
        let denom = &solar_rad - &corrected_thermal;

        let mut result = Array1::zeros(n);
        for i in 0..n {
            let dn = denom[i];
            if dn < EPSILON {
                result[i] = f64::NAN;
                continue;
            }
            let sz = sun_zenith[i];
            if let Some(limit) = self.masking_limit {
                if sz < 0.0 || sz > limit {
                    result[i] = f64::NAN;
                    continue;
                }
            }
            if tb_near_ir[i].is_nan() {
                result[i] = f64::NAN;
                continue;
            }
            let val = nomin[i] / dn;
            result[i] = val.clamp(0.0, 1.0);
        }
        result
    }

    pub fn emissive_part(&self, tb_nir: &Array1<f64>) -> Array1<f64> {
        self.tb2radiance(tb_nir)
    }
}

pub fn get_as_array(value: f64, n: usize) -> Array1<f64> {
    Array1::from_elem(n, value)
}

pub fn from_scalar(a: f64, n: usize) -> Array1<f64> {
    Array1::from_elem(n, a)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use ndarray::arr1;

    #[test]
    fn test_terminator_limit_value() {
        assert_relative_eq!(TERMINATOR_LIMIT, 85.0);
    }

    #[test]
    fn test_reflectance_terminator_limit() {
        let calc = ReflectanceCalculator::new("test", "test")
            .with_central_wavelength(3.78e-6)
            .with_solar_flux(2.0);
        let sunz = arr1(&[TERMINATOR_LIMIT]);
        let tb_nir = arr1(&[290.0]);
        let tb_thermal = arr1(&[280.0]);
        let refl = calc.reflectance_from_tbs(&sunz, &tb_nir, &tb_thermal, None);
        assert!(refl[0].is_nan());
    }

    #[test]
    fn test_reflectance_baseline() {
        let wvl = Array1::linspace(3.6e-6, 3.95e-6, 36);
        let resp = Array1::ones(36);
        let calc = ReflectanceCalculator::new("test", "test")
            .with_rsr(wvl, resp)
            .with_solar_flux(2.0);
        let sunz = arr1(&[30.0]);
        let tb_nir = arr1(&[300.0]);
        let tb_thermal = arr1(&[290.0]);
        let refl = calc.reflectance_from_tbs(&sunz, &tb_nir, &tb_thermal, None);
        assert!(refl[0] > 0.0);
        assert!(refl[0] < 1.0);
    }

    #[test]
    fn test_reflectance_clamped() {
        let wvl = Array1::linspace(3.6e-6, 3.95e-6, 36);
        let resp = Array1::ones(36);
        let calc = ReflectanceCalculator::new("test", "test")
            .with_rsr(wvl, resp)
            .with_solar_flux(100.0);
        let sunz = arr1(&[0.0]);
        let tb_nir = arr1(&[350.0]);
        let tb_thermal = arr1(&[200.0]);
        let refl = calc.reflectance_from_tbs(&sunz, &tb_nir, &tb_thermal, None);
        assert!(refl[0] <= 1.0);
    }

    #[test]
    fn test_solar_radiance() {
        let calc = ReflectanceCalculator::new("test", "test").with_solar_flux(1400.0);
        let sunz = arr1(&[0.0]);
        let sr = calc.solar_radiance(&sunz);
        assert_relative_eq!(sr[0], 1400.0 / PI, epsilon = 1e-10);
    }

    #[test]
    fn test_solar_radiance_sunz60() {
        let calc = ReflectanceCalculator::new("test", "test").with_solar_flux(1400.0);
        let sunz = arr1(&[60.0]);
        let sr = calc.solar_radiance(&sunz);
        assert_relative_eq!(sr[0], 1400.0 * 0.5 / PI, epsilon = 1e-10);
    }

    #[test]
    fn test_derive_rad39_corr() {
        let calc = ReflectanceCalculator::new("test", "test");
        let bt11 = arr1(&[290.0]);
        let bt13 = arr1(&[270.0]);
        let corr = calc.derive_rad39_corr(&bt11, &bt13);
        assert!(corr[0] > 0.0);
    }

    #[test]
    fn test_with_rsr() {
        let wvl = Array1::linspace(3.6e-6, 3.95e-6, 10);
        let resp = Array1::ones(10);
        let calc = ReflectanceCalculator::new("test", "test")
            .with_rsr(wvl, resp)
            .with_solar_flux(2.0);
        assert!(calc.rsr_integral > 0.0);
    }
}
