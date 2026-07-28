use ndarray::Array1;

use crate::rsr::Rsr;
use crate::utils::trapezoid;

pub struct SolarIrradianceSpectrum {
    pub wavelength: Array1<f64>,
    pub irradiance: Array1<f64>,
    pub ipol_wavelength: Option<Array1<f64>>,
    pub ipol_irradiance: Option<Array1<f64>>,
    pub wavenumber: Option<Array1<f64>>,
    wavespace: WaveSpace,
    dlambda: f64,
}

enum WaveSpace {
    Wavelength,
    Wavenumber,
}

impl SolarIrradianceSpectrum {
    pub fn new(filename: &str, dlambda: f64) -> Self {
        let mut wavelength = Vec::new();
        let mut irradiance = Vec::new();

        let content =
            std::fs::read_to_string(filename).expect("Failed to read solar spectrum file");
        for line in content.lines() {
            let line = line.trim();
            if line.is_empty() {
                continue;
            }
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 2 {
                if let (Ok(w), Ok(i)) = (parts[0].parse::<f64>(), parts[1].parse::<f64>()) {
                    wavelength.push(w);
                    irradiance.push(i);
                }
            }
        }

        SolarIrradianceSpectrum {
            wavelength: Array1::from_vec(wavelength),
            irradiance: Array1::from_vec(irradiance),
            ipol_wavelength: None,
            ipol_irradiance: None,
            wavenumber: None,
            wavespace: WaveSpace::Wavelength,
            dlambda,
        }
    }

    pub fn solar_constant(&self) -> f64 {
        match self.wavespace {
            WaveSpace::Wavelength => trapezoid(&self.irradiance, &self.wavelength),
            WaveSpace::Wavenumber => {
                if let Some(ref wn) = self.wavenumber {
                    trapezoid(&self.irradiance, wn)
                } else {
                    0.0
                }
            }
        }
    }

    pub fn inband_solarflux(&mut self, rsr: &Rsr, scale: f64) -> f64 {
        self._band_calculations(rsr, true, scale)
    }

    pub fn inband_solarirradiance(&mut self, rsr: &Rsr, scale: f64) -> f64 {
        self._band_calculations(rsr, false, scale)
    }

    pub fn interpolate(&mut self, dlambda: f64, ival_wavelength: Option<(f64, f64)>) {
        self.dlambda = dlambda;

        let (start, end) = match ival_wavelength {
            Some((s, e)) => (s, e),
            None => match self.wavespace {
                WaveSpace::Wavelength => (
                    self.wavelength[0],
                    self.wavelength[self.wavelength.len() - 1],
                ),
                WaveSpace::Wavenumber => {
                    if let Some(ref wn) = self.wavenumber {
                        (wn[0], wn[wn.len() - 1])
                    } else {
                        return;
                    }
                }
            },
        };

        let n = ((end - start) / dlambda).round() as usize + 1;
        let xspl: Vec<f64> = (0..n).map(|i| start + i as f64 * dlambda).collect();
        let xspl = Array1::from_vec(xspl);

        let (src_x, src_y) = match self.wavespace {
            WaveSpace::Wavelength => (&self.wavelength, &self.irradiance),
            WaveSpace::Wavenumber => {
                if let Some(ref wn) = self.wavenumber {
                    (wn, &self.irradiance)
                } else {
                    return;
                }
            }
        };

        let yspl = Array1::from_vec(cubic_spline_interpolate(&xspl, src_x, src_y));

        self.ipol_wavelength = Some(xspl);
        self.ipol_irradiance = Some(yspl);
    }

    pub fn set_wavespace_wavenumber(&mut self) {
        let n = self.wavelength.len();
        let mut wavenumber = Array1::zeros(n);
        let mut irradiance = Array1::zeros(n);

        for i in 0..n {
            let j = n - 1 - i;
            wavenumber[i] = 1.0 / (1e-4 * self.wavelength[j]);
            irradiance[i] = self.irradiance[j] * self.wavelength[j] * self.wavelength[j] * 0.1;
        }

        self.wavenumber = Some(wavenumber);
        self.irradiance = irradiance;
        self.wavespace = WaveSpace::Wavenumber;
    }

    fn _band_calculations(&mut self, rsr: &Rsr, flux: bool, scale: f64) -> f64 {
        let wvl = &rsr.wavelength * scale;
        let resp = &rsr.response;

        let start = wvl[0];
        let end = wvl[wvl.len() - 1];

        self.interpolate(self.dlambda, Some((start, end)));

        let ipol_w = self.ipol_wavelength.as_ref().unwrap();
        let ipol_i = self.ipol_irradiance.as_ref().unwrap();

        let n_ipol = ipol_w.len();
        let n_expected = ((end - start) / self.dlambda).round() as usize + 1;
        let capacity = n_expected.min(n_ipol);

        let mut masked_wvl = Vec::with_capacity(capacity);
        let mut masked_irr = Vec::with_capacity(capacity);
        let mut masked_resp = Vec::with_capacity(capacity);

        let resp_ipol = cubic_spline_interpolate(ipol_w, &wvl, resp);

        for i in 0..n_ipol {
            let w = ipol_w[i];
            if w >= start && w <= end {
                masked_wvl.push(w);
                masked_irr.push(ipol_i[i]);
                masked_resp.push(resp_ipol[i]);
            }
        }

        let mw = Array1::from_vec(masked_wvl);
        let mi = Array1::from_vec(masked_irr);
        let mr = Array1::from_vec(masked_resp);

        let product = &mi * &mr;
        let integrated = trapezoid(&product, &mw);

        if flux {
            integrated
        } else {
            let resp_sum = trapezoid(&mr, &mw);
            integrated / resp_sum
        }
    }
}

fn linear_interpolate(x_new: &Array1<f64>, x_vals: &Array1<f64>, y_vals: &Array1<f64>) -> Vec<f64> {
    if x_vals.len() < 2 {
        return vec![y_vals[0]; x_new.len()];
    }

    let x_slice = x_vals.as_slice().unwrap();
    let y_slice = y_vals.as_slice().unwrap();

    x_new
        .iter()
        .map(|&x| {
            if x <= x_slice[0] {
                return y_slice[0];
            }
            if x >= x_slice[x_slice.len() - 1] {
                return y_slice[x_slice.len() - 1];
            }

            let idx = match x_slice.binary_search_by(|&v| v.partial_cmp(&x).unwrap()) {
                Ok(i) => i,
                Err(i) => i.saturating_sub(1),
            };

            let i = idx.min(x_slice.len() - 2);
            let x0 = x_slice[i];
            let x1 = x_slice[i + 1];
            let y0 = y_slice[i];
            let y1 = y_slice[i + 1];

            if (x1 - x0).abs() < 1e-15 {
                return y0;
            }

            y0 + (y1 - y0) * (x - x0) / (x1 - x0)
        })
        .collect()
}

fn cubic_spline_interpolate(
    x_new: &Array1<f64>,
    x_vals: &Array1<f64>,
    y_vals: &Array1<f64>,
) -> Vec<f64> {
    let n = x_vals.len();
    if n < 2 {
        return vec![y_vals[0]; x_new.len()];
    }
    if n == 2 {
        return linear_interpolate(x_new, x_vals, y_vals);
    }

    let x_slice = x_vals.as_slice().unwrap();
    let y_slice = y_vals.as_slice().unwrap();

    let m = n - 2;
    let mut h = vec![0.0; n - 1];
    for i in 0..n - 1 {
        h[i] = x_slice[i + 1] - x_slice[i];
        if h[i].abs() < 1e-15 {
            return linear_interpolate(x_new, x_vals, y_vals);
        }
    }

    let mut a = vec![0.0; m]; // subdiagonal
    let mut b = vec![0.0; m]; // main diagonal
    let mut c = vec![0.0; m]; // superdiagonal
    let mut d = vec![0.0; m]; // RHS

    // Natural spline: M_0 = 0, M_{n-1} = 0
    // m = n-2 equations, unknowns M_1..M_{n-2}
    for k in 0..m {
        // Equation for M_{k+1} (i = k+1)
        let i = k + 1;
        // h_{i-1} * M_{i-1} + 2*(h_{i-1}+h_i) * M_i + h_i * M_{i+1} = 6*(dy_i/h_i - dy_{i-1}/h_{i-1})
        a[k] = h[i - 1];
        b[k] = 2.0 * (h[i - 1] + h[i]);
        c[k] = h[i];
        d[k] =
            6.0 * ((y_slice[i + 1] - y_slice[i]) / h[i] - (y_slice[i] - y_slice[i - 1]) / h[i - 1]);
    }

    // Thomas algorithm (tridiagonal solver)
    // Forward sweep: eliminate lower diagonal
    for k in 1..m {
        let w = a[k] / b[k - 1];
        b[k] -= w * c[k - 1];
        d[k] -= w * d[k - 1];
    }

    // Back substitution
    let mut mm = vec![0.0; n]; // M values, index 0..n-1, M_0=0, M_{n-1}=0 by natural BC
    mm[m - 1] = d[m - 1] / b[m - 1];
    for k in (0..m - 1).rev() {
        mm[k + 1] = (d[k] - c[k] * mm[k + 2]) / b[k];
    }

    // Evaluate spline at each output point
    x_new
        .iter()
        .map(|&xq| {
            if xq <= x_slice[0] {
                return y_slice[0];
            }
            if xq >= x_slice[n - 1] {
                return y_slice[n - 1];
            }

            let idx = match x_slice.binary_search_by(|&v| v.partial_cmp(&xq).unwrap()) {
                Ok(i) => i,
                Err(i) => i.saturating_sub(1),
            };
            let i = idx.min(n - 2);
            let hi = x_slice[i + 1] - x_slice[i];
            let dx1 = x_slice[i + 1] - xq;
            let dx0 = xq - x_slice[i];

            let a_term = mm[i] * dx1.powi(3) / (6.0 * hi);
            let b_term = mm[i + 1] * dx0.powi(3) / (6.0 * hi);
            let c_term = (y_slice[i] - mm[i] * hi.powi(2) / 6.0) * dx1 / hi;
            let d_term = (y_slice[i + 1] - mm[i + 1] * hi.powi(2) / 6.0) * dx0 / hi;

            a_term + b_term + c_term + d_term
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rsr::Rsr;
    use approx::{assert_abs_diff_eq, assert_relative_eq};
    use ndarray::array;

    const SOLAR_SPECTRUM_FILE: &str = "data/e490_00a.dat";

    fn test_rsr() -> Rsr {
        let wvl = array![
            3.6123999,
            3.6163599,
            3.6264927,
            3.6363862,
            3.646468,
            3.6564937,
            3.6664478,
            3.6765388,
            3.6865413,
            3.6964585,
            3.7065142,
            3.716509,
            3.7264658,
            3.7364102,
            3.7463682,
            3.7563652,
            3.7664226,
            3.7763396,
            3.7863384,
            3.7964207,
            3.8063589,
            3.8163606,
            3.8264089,
            3.8364836,
            3.8463381,
            3.8563975,
            3.8664163,
            3.8763755,
            3.8864797,
            3.8964978,
            3.9064275,
            3.9164873,
            3.9264729,
            3.9364026,
            3.9465107,
            3.9535347_f64,
        ];
        let resp = array![
            0.01, 0.0118, 0.01987, 0.03226, 0.05028, 0.0849, 0.16645, 0.33792, 0.59106, 0.81815,
            0.96077, 0.92855, 0.86008, 0.8661, 0.87697, 0.85412, 0.88922, 0.9541, 0.95687, 0.91037,
            0.91058, 0.94256, 0.94719, 0.94808, 1.0, 0.92676, 0.67429, 0.44715, 0.27762, 0.14852,
            0.07141, 0.04151, 0.02925, 0.02085, 0.01414, 0.01_f64,
        ];
        Rsr::new(wvl, resp)
    }

    #[test]
    fn test_load_solar_spectrum() {
        let solar = SolarIrradianceSpectrum::new(SOLAR_SPECTRUM_FILE, 0.005);
        assert_eq!(solar.wavelength.len(), 1697);
        assert_eq!(solar.irradiance.len(), 1697);
        assert!(solar.wavelength[0] < solar.wavelength[1]);
    }

    #[test]
    fn test_solar_constant() {
        let solar = SolarIrradianceSpectrum::new(SOLAR_SPECTRUM_FILE, 0.005);
        let sc = solar.solar_constant();
        assert!(sc > 1364.0 && sc < 1368.0, "sc = {}", sc);
    }

    #[test]
    fn test_interpolate() {
        let mut solar = SolarIrradianceSpectrum::new(SOLAR_SPECTRUM_FILE, 0.005);
        solar.interpolate(0.001, Some((0.200, 0.240)));
        let ipol = solar.ipol_wavelength.as_ref().unwrap();
        assert!(!ipol.is_empty());
        assert_relative_eq!(ipol[0], 0.200, epsilon = 1e-6);
        assert_relative_eq!(ipol[ipol.len() - 1], 0.240, epsilon = 1e-6);
        let expected_n = f64::round((0.240 - 0.200) / 0.001) as usize + 1;
        assert_eq!(ipol.len(), expected_n);
    }

    #[test]
    fn test_inband_solarflux() {
        let mut solar = SolarIrradianceSpectrum::new(SOLAR_SPECTRUM_FILE, 0.005);
        let rsr = test_rsr();
        let flux = solar.inband_solarflux(&rsr, 1.0);
        assert_abs_diff_eq!(flux, 2.002_927_627, epsilon = 1e-3);
    }

    #[test]
    fn test_inband_solarflux_wavenumber() {
        let mut solar = SolarIrradianceSpectrum::new(SOLAR_SPECTRUM_FILE, 0.005);
        solar.set_wavespace_wavenumber();
        let rsr = test_rsr();
        let flux = solar.inband_solarflux(&rsr, 1.0);
        assert!(flux > 0.0);
    }

    #[test]
    fn test_inband_solarirradiance() {
        let mut solar = SolarIrradianceSpectrum::new(SOLAR_SPECTRUM_FILE, 0.005);
        let rsr = test_rsr();
        let irradiance = solar.inband_solarirradiance(&rsr, 1.0);
        assert!(irradiance > 0.0);
        let flux = solar.inband_solarflux(&rsr, 1.0);
        assert!(irradiance != flux);
    }

    #[test]
    fn test_solar_constant_wavenumber() {
        let mut solar = SolarIrradianceSpectrum::new(SOLAR_SPECTRUM_FILE, 0.005);
        solar.set_wavespace_wavenumber();
        let sc = solar.solar_constant();
        assert!(sc > 0.0);
    }
}
