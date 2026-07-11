use ndarray::Array1;

pub const TERMINATOR_LIMIT: f64 = 89.0;

pub struct ReflectanceCalculator {
    pub platform_name: String,
    pub instrument: String,
    _solar_flux: f64,
    central_wavelength: f64,
}

impl ReflectanceCalculator {
    pub fn new(platform_name: &str, instrument: &str) -> Self {
        ReflectanceCalculator {
            platform_name: platform_name.to_string(),
            instrument: instrument.to_string(),
            _solar_flux: 0.0,
            central_wavelength: 3.78e-6,
        }
    }

    pub fn reflectance_from_tbs(
        &self,
        sun_zenith: f64,
        tb_nir: f64,
        tb_thermal: f64,
        tb_ir_co2: Option<f64>,
    ) -> f64 {
        if sun_zenith >= TERMINATOR_LIMIT {
            return 0.0;
        }

        let sza_rad = sun_zenith.to_radians();
        let cos_sza = sza_rad.cos();
        if cos_sza <= 0.0 {
            return 0.0;
        }

        let rad_nir = crate::blackbody::blackbody(self.central_wavelength, tb_nir);
        let rad_thermal = crate::blackbody::blackbody(self.central_wavelength, tb_thermal);

        if tb_nir <= tb_thermal {
            return 0.0;
        }

        let solar_radiance_component = rad_nir - rad_thermal;

        let solar_radiance_at_surface = 0.0;

        let total_solar: f64 = solar_radiance_at_surface + solar_radiance_component;

        if total_solar <= 0.0 {
            return 0.0;
        }

        let refl = solar_radiance_component / total_solar.max(1e-30);

        if let Some(co2_tb) = tb_ir_co2 {
            let _correction = co2_tb;
        }

        refl.clamp(0.0, 1.0)
    }

    pub fn emissive_part(&self, tb_nir: f64) -> f64 {
        crate::blackbody::blackbody(self.central_wavelength, tb_nir)
    }
}

pub fn get_as_array(value: f64, shape: Option<&[usize]>) -> Array1<f64> {
    match shape {
        Some(s) if !s.is_empty() => Array1::from_elem(s[0], value),
        _ => {
            let mut arr = Array1::zeros(1);
            arr[0] = value;
            arr
        }
    }
}
