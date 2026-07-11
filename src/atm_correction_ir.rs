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

    pub fn get_correction(
        &self,
        sat_zenith: &Array2<f64>,
        _bandname: &str,
        data: &Array2<f64>,
    ) -> Array2<f64> {
        viewzen_corr(data, sat_zenith)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    fn satz_data() -> Array2<f64> {
        let mut rows = Vec::new();
        for i in 0..10 {
            let base = 48.03 + i as f64 * 0.06;
            let row: Vec<f64> = (0..10).map(|j| base + j as f64 * 0.00002).collect();
            rows.extend(row);
        }
        Array2::from_shape_vec((10, 10), rows).unwrap()
    }

    fn tbs_data() -> Array2<f64> {
        let mut rows = Vec::new();
        for i in 0..10 {
            let base = 284.04 + i as f64 * 0.08;
            let row: Vec<f64> = (0..10).map(|j| base + j as f64 * 0.000_026_67).collect();
            rows.extend(row);
        }
        Array2::from_shape_vec((10, 10), rows).unwrap()
    }

    fn expected_result() -> Array2<f64> {
        let rows: Vec<f64> = vec![
            286.031_594_12, 286.031_624_17, 286.031_654_21, 286.031_684_26, 286.031_714_30,
            286.031_744_34, 286.031_774_39, 286.031_804_43, 286.031_834_47, 286.031_864_52,
            286.121_747_23, 286.121_777_29, 286.121_807_35, 286.121_837_41, 286.121_867_47,
            286.121_897_52, 286.121_927_58, 286.121_957_64, 286.121_987_70, 286.122_017_76,
            286.211_945_45, 286.211_975_52, 286.212_005_60, 286.212_035_67, 286.212_065_74,
            286.212_095_82, 286.212_125_89, 286.212_155_97, 286.212_186_04, 286.212_216_11,
            286.302_188_96, 286.302_219_05, 286.302_249_13, 286.302_279_22, 286.302_309_31,
            286.302_339_40, 286.302_369_49, 286.302_399_58, 286.302_429_67, 286.302_459_76,
            286.392_477_93, 286.392_508_03, 286.392_538_14, 286.392_568_24, 286.392_598_34,
            286.392_628_45, 286.392_658_55, 286.392_688_66, 286.392_718_76, 286.392_748_86,
            286.482_812_54, 286.482_842_66, 286.482_872_78, 286.482_902_90, 286.482_933_02,
            286.482_963_14, 286.482_993_25, 286.483_023_37, 286.483_053_49, 286.483_083_61,
            286.573_192_97, 286.573_223_10, 286.573_253_24, 286.573_283_37, 286.573_313_51,
            286.573_343_64, 286.573_373_78, 286.573_403_91, 286.573_434_05, 286.573_464_18,
            286.663_619_40, 286.663_649_55, 286.663_679_70, 286.663_709_85, 286.663_740_00,
            286.663_770_15, 286.663_800_30, 286.663_830_45, 286.663_860_60, 286.663_890_75,
            286.754_092_00, 286.754_122_16, 286.754_152_33, 286.754_182_49, 286.754_212_66,
            286.754_242_83, 286.754_272_99, 286.754_303_16, 286.754_333_32, 286.754_363_49,
            286.844_610_96, 286.844_641_14, 286.844_671_32, 286.844_701_50, 286.844_731_68,
            286.844_761_86, 286.844_792_04, 286.844_822_22, 286.844_852_40, 286.844_882_58,
        ];
        Array2::from_shape_vec((10, 10), rows).unwrap()
    }

    #[test]
    fn test_viewzen_corr_accuracy() {
        let satz = satz_data();
        let tbs = tbs_data();
        let result = viewzen_corr(&tbs, &satz);
        let expected = expected_result();
        for i in 0..10 {
            for j in 0..10 {
                assert_abs_diff_eq!(result[(i, j)], expected[(i, j)], epsilon = 1e-4);
            }
        }
    }

    #[test]
    fn test_viewzen_corr_matches_python() {
        let satz = satz_data();
        let tbs = tbs_data();
        let result = viewzen_corr(&tbs, &satz);
        let expected = expected_result();
        let diff = (&result - &expected).mapv(|x| x.abs());
        let max_diff = diff.fold(f64::NEG_INFINITY, |a, &b| a.max(b));
        assert!(max_diff < 1e-3, "max diff = {}", max_diff);
    }
}
