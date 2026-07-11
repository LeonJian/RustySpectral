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

        let content = std::fs::read_to_string(filename).expect("Failed to read solar spectrum file");
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
                WaveSpace::Wavelength => (self.wavelength[0], self.wavelength[self.wavelength.len() - 1]),
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

        let yspl = Array1::from_vec(linear_interpolate(&xspl, src_x, src_y));

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
        let n = ((end - start) / self.dlambda).round() as usize + 1;
        let xspl: Vec<f64> = (0..n).map(|i| start + i as f64 * self.dlambda).collect();
        let xspl_a = Array1::from_vec(xspl);

        let resp_ipol = Array1::from_vec(linear_interpolate(&xspl_a, &wvl, resp));

        self.interpolate(self.dlambda, Some((start, end)));

        let ipol_w = self.ipol_wavelength.as_ref().unwrap();
        let ipol_i = self.ipol_irradiance.as_ref().unwrap();

        let mask: Vec<bool> = ipol_w
            .iter()
            .map(|&x| x >= start && x <= end)
            .collect();

        let mut masked_wvl = Vec::new();
        let mut masked_irr = Vec::new();
        let mut masked_resp = Vec::new();

        let mut resp_idx = 0usize;
        for (i, &m) in mask.iter().enumerate() {
            if m {
                masked_wvl.push(ipol_w[i]);
                masked_irr.push(ipol_i[i]);
                // Map back to response grid using nearest-neighbor or linear
                if i < resp_ipol.len() {
                    masked_resp.push(resp_ipol[i]);
                } else {
                    masked_resp.push(resp_ipol[resp_idx.min(resp_ipol.len() - 1)]);
                }
                resp_idx += 1;
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
