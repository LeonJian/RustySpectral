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

#[inline]
pub fn planck(wave: f64, temperature: f64) -> f64 {
    if temperature.abs() <= EPSILON {
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

#[inline]
pub fn planck_wn(wavenumber: f64, temperature: f64) -> f64 {
    if temperature.abs() <= EPSILON {
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

#[inline]
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

#[inline]
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

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use ndarray::arr2;

    const RAD_11MICRON_300K: f64 = 9573176.935507433;
    const RAD_11MICRON_301K: f64 = 9714686.576498277;
    const WN_RAD_11MICRON_300K: f64 = 0.00115835441353;
    const WN_RAD_11MICRON_301K: f64 = 0.00117547716523;

    const H: f64 = 6.626_069_57e-34;
    const K: f64 = 1.380_648_8e-23;
    const C: f64 = 2.997_924_58e8;

    #[test]
    fn test_physical_constants() {
        assert_abs_diff_eq!(H_PLANCK, H, epsilon = 1e-40);
        assert_abs_diff_eq!(K_BOLTZMANN, K, epsilon = 1e-30);
        assert_abs_diff_eq!(C_SPEED, C, epsilon = 1e-1);
    }

    #[test]
    fn test_planck_wavelength_300k() {
        let rad = planck(11e-6, 300.0);
        assert_abs_diff_eq!(rad, RAD_11MICRON_300K, epsilon = 1e-4);
    }

    #[test]
    fn test_planck_wavelength_301k() {
        let rad = planck(11e-6, 301.0);
        assert_abs_diff_eq!(rad, RAD_11MICRON_301K, epsilon = 1e-4);
    }

    #[test]
    fn test_planck_wn_300k() {
        let rad = planck_wn(90909.1, 300.0);
        assert_abs_diff_eq!(rad, WN_RAD_11MICRON_300K, epsilon = 3e-10);
    }

    #[test]
    fn test_planck_wn_301k() {
        let rad = planck_wn(90909.1, 301.0);
        assert_abs_diff_eq!(rad, WN_RAD_11MICRON_301K, epsilon = 3e-10);
    }

    #[test]
    fn test_rad2temp_roundtrip() {
        let wavel = 11e-6;
        let rad = planck(wavel, 300.0);
        let t = blackbody_rad2temp(wavel, rad);
        assert_abs_diff_eq!(t, 300.0, epsilon = 1e-8);

        let rad = planck(wavel, 301.0);
        let t = blackbody_rad2temp(wavel, rad);
        assert_abs_diff_eq!(t, 301.0, epsilon = 1e-8);
    }

    #[test]
    fn test_rad2temp_wn_roundtrip() {
        let wn = 90909.1;
        let rad = planck_wn(wn, 300.0);
        let t = blackbody_wn_rad2temp(wn, rad);
        assert_abs_diff_eq!(t, 300.0, epsilon = 1e-8);

        let rad = planck_wn(wn, 301.0);
        let t = blackbody_wn_rad2temp(wn, rad);
        assert_abs_diff_eq!(t, 301.0, epsilon = 1e-8);
    }

    #[test]
    fn test_planck_wavelength_array() {
        let temps = arr2(&[[300.0, 301.0], [299.0, 298.0], [279.0, 286.0]]);
        let result = planck_array_wavelength(10e-6, &temps);
        assert_eq!(result.shape(), &[3, 2]);
        assert!(result[[0, 0]] > result[[1, 1]]);
        assert!(result[[0, 1]] > result[[0, 0]]);
    }

    #[test]
    fn test_rad2temp_array() {
        let radiances = arr2(&[[0.001, 0.0009], [0.0012, 0.0018]]);
        let temps = blackbody_wn_rad2temp_array(90909.1, &radiances);
        let expected = arr2(&[
            [290.327_691_6, 283.761_154_41],
            [302.418_133_0, 333.141_416_4],
        ]);
        assert_abs_diff_eq!(temps[[0, 0]], expected[[0, 0]], epsilon = 1e-5);
        assert_abs_diff_eq!(temps[[1, 1]], expected[[1, 1]], epsilon = 1e-5);
    }

    #[test]
    fn test_zero_radiance_returns_nan() {
        assert!(blackbody_rad2temp(11e-6, 0.0).is_nan());
    }

    #[test]
    fn test_zero_temperature_returns_nan() {
        assert!(planck(11e-6, 0.0).is_nan());
    }

    #[test]
    fn test_planck_wn_zero_temperature_returns_nan() {
        assert!(planck_wn(90909.1, 0.0).is_nan());
    }

    #[test]
    fn test_wn_rad2temp_zero_returns_nan() {
        assert!(blackbody_wn_rad2temp(90909.1, 0.0).is_nan());
    }
}
