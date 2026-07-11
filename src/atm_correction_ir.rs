use ndarray::Array2;

pub fn viewzen_corr(data: &Array2<f64>, view_zen: &Array2<f64>) -> Array2<f64> {
    let mut result = data.clone();

    fn ratio(value: f64, v_null: f64, v_ref: f64) -> f64 {
        (value - v_null) / (v_ref - v_null)
    }

    fn tau0(t: f64) -> f64 {
        let t0: f64 = 210.0;
        let t_ref: f64 = 320.0;
        let tau_ref: f64 = 9.85;
        (1.0_f64 + tau_ref).powf(ratio(t, t0, t_ref)) - 1.0
    }

    fn tau(t: f64) -> f64 {
        let t0 = 170.0;
        let t_ref = 295.0;
        let tau_ref = 1.0;
        let m = 4;
        tau_ref * ratio(t, t0, t_ref).powi(m)
    }

    fn delta(z: f64) -> f64 {
        let z0: f64 = 0.0;
        let z_ref: f64 = 70.0;
        let delta_ref: f64 = 6.2;
        (1.0_f64 + delta_ref).powf(ratio(z, z0, z_ref)) - 1.0
    }

    for ((i, j), &z) in view_zen.indexed_iter() {
        if z == 0.0 {
            result[(i, j)] += tau0(data[(i, j)]);
        } else if z > 0.0 && z < 90.0 {
            result[(i, j)] += tau(data[(i, j)]) * delta(z);
        }
    }

    result
}

pub struct AtmosphericalCorrection {
    pub platform_name: String,
    pub sensor: String,
}

impl AtmosphericalCorrection {
    pub fn new(platform_name: &str, sensor: &str) -> Self {
        AtmosphericalCorrection {
            platform_name: platform_name.to_string(),
            sensor: sensor.to_string(),
        }
    }

    pub fn get_correction(&self, sat_zenith: &Array2<f64>, data: &Array2<f64>) -> Array2<f64> {
        viewzen_corr(data, sat_zenith)
    }
}
