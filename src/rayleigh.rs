use ndarray::{Array1, Array3, Array4, Axis};

pub fn clip_angles_inside_coordinate_range(zenith_angle: &Array1<f64>, zenith_secant_max: f64) -> Array1<f64> {
    let clip_angle = (1.0 / zenith_secant_max).acos().to_degrees();
    zenith_angle.mapv(|z| if z.is_nan() { 0.0 } else { z.clamp(0.0, clip_angle) })
}

pub fn clip_angles_inside_coordinate_range_scalar(zenith_angle: f64, zenith_secant_max: f64) -> f64 {
    if zenith_angle.is_nan() {
        0.0f64
    } else {
        let clip_angle = (1.0 / zenith_secant_max).acos().to_degrees();
        zenith_angle.clamp(0.0, clip_angle)
    }
}

pub fn reduce_rayleigh_highzenith(
    zenith: &Array1<f64>,
    rayref: &Array1<f64>,
    thresh_zen: f64,
    maxzen: f64,
    strength: f64,
) -> Array1<f64> {
    let factor: Array1<f64> = zenith.mapv(|z| {
        if z < thresh_zen { 0.0 } else { (z - thresh_zen) / (maxzen - thresh_zen) }
    });
    let factor = 1.0 - strength * &factor;
    let factor = factor.mapv(|f| f.clamp(0.0, 1.0));
    rayref * &factor
}

pub fn get_wavelength_index_and_factor(wvl_coord: &Array1<f64>, wvl: f64) -> (usize, f64) {
    let idx = match wvl_coord.iter().position(|&v| v > wvl) {
        Some(i) => i,
        None => wvl_coord.len() - 1,
    };
    let wavelength_index = idx;
    let wvl1 = wvl_coord[wavelength_index - 1];
    let wvl2 = wvl_coord[wavelength_index];
    let wavelength_factor = (wvl2 - wvl) / (wvl2 - wvl1);
    (wavelength_index, wavelength_factor)
}

pub fn get_wavelength_adjusted_lut(
    rayleigh_refl: &Array4<f64>,
    wvl_coord: &Array1<f64>,
    wvl: f64,
) -> Array3<f64> {
    let (wi, wf) = get_wavelength_index_and_factor(wvl_coord, wvl);
    let slice1 = rayleigh_refl.index_axis(Axis(0), wi - 1);
    let slice2 = rayleigh_refl.index_axis(Axis(0), wi);
    wf * &slice1 + (1.0 - wf) * &slice2
}

pub fn trilinear_interpolate(
    grid: &Array3<f64>,
    sunz_sec: f64,
    azidiff_in: f64,
    satz_sec: f64,
    sunz_coord: &Array1<f64>,
    azid_coord: &Array1<f64>,
    satz_coord: &Array1<f64>,
) -> f64 {
    let azidiff = 180.0 - azidiff_in;

    let si = find_interval_index(sunz_coord, sunz_sec);
    let ai = find_interval_index(azid_coord, azidiff);
    let ti = find_interval_index(satz_coord, satz_sec);

    let s0 = sunz_coord[si];
    let s1 = sunz_coord[si + 1];
    let a0 = azid_coord[ai];
    let a1 = azid_coord[ai + 1];
    let t0 = satz_coord[ti];
    let t1 = satz_coord[ti + 1];

    let sd = (sunz_sec - s0) / (s1 - s0);
    let ad = (azidiff - a0) / (a1 - a0);
    let td = (satz_sec - t0) / (t1 - t0);

    let c000 = grid[(si, ai, ti)];
    let c001 = grid[(si, ai, ti + 1)];
    let c010 = grid[(si, ai + 1, ti)];
    let c011 = grid[(si, ai + 1, ti + 1)];
    let c100 = grid[(si + 1, ai, ti)];
    let c101 = grid[(si + 1, ai, ti + 1)];
    let c110 = grid[(si + 1, ai + 1, ti)];
    let c111 = grid[(si + 1, ai + 1, ti + 1)];

    let c00 = c000 * (1.0 - td) + c001 * td;
    let c01 = c010 * (1.0 - td) + c011 * td;
    let c10 = c100 * (1.0 - td) + c101 * td;
    let c11 = c110 * (1.0 - td) + c111 * td;

    let c0 = c00 * (1.0 - ad) + c01 * ad;
    let c1 = c10 * (1.0 - ad) + c11 * ad;

    (c0 * (1.0 - sd) + c1 * sd) * 100.0
}

fn find_interval_index(coords: &Array1<f64>, value: f64) -> usize {
    let idx = coords.iter().position(|&v| v > value).unwrap_or(coords.len() - 1);
    (idx.saturating_sub(1)).min(coords.len() - 2)
}

pub fn rayleigh_interpolate_by_angles(
    sun_zenith: &Array1<f64>,
    sat_zenith: &Array1<f64>,
    azidiff: &Array1<f64>,
    rayleigh_refl: &Array4<f64>,
    wvl_coord: &Array1<f64>,
    wvl: f64,
    sunz_sec_coord: &Array1<f64>,
    satz_sec_coord: &Array1<f64>,
    azid_coord: &Array1<f64>,
) -> Array1<f64> {
    let grid3 = get_wavelength_adjusted_lut(rayleigh_refl, wvl_coord, wvl);

    let n = sun_zenith.len();
    let mut result = Array1::zeros(n);
    for i in 0..n {
        let sz = clip_angles_inside_coordinate_range_scalar(sun_zenith[i], sunz_sec_coord[sunz_sec_coord.len() - 1]);
        let satz = clip_angles_inside_coordinate_range_scalar(sat_zenith[i], satz_sec_coord[satz_sec_coord.len() - 1]);
        let sunzsec = 1.0 / sz.to_radians().cos();
        let satzsec = 1.0 / satz.to_radians().cos();

        result[i] = trilinear_interpolate(
            &grid3, sunzsec, azidiff[i], satzsec,
            sunz_sec_coord, azid_coord, satz_sec_coord,
        );
    }
    result
}

pub fn normalize_sensor(platform_name: &str, sensor: &str) -> String {
    let instruments: std::collections::HashMap<&str, &str> = std::collections::HashMap::from([
        ("Envisat", "aatsr"), ("GOES-16", "abi"), ("GOES-17", "abi"), ("GOES-18", "abi"),
        ("GOES-19", "abi"), ("FY-4A", "agri"), ("FY-4B", "agri"), ("Himawari-8", "ahi"),
        ("Himawari-9", "ahi"), ("GEO-KOMPSAT-2A", "ami"), ("NOAA-10", "avhrr1"),
        ("NOAA-6", "avhrr1"), ("NOAA-8", "avhrr1"), ("TIROS-N", "avhrr1"),
        ("NOAA-11", "avhrr2"), ("NOAA-12", "avhrr2"), ("NOAA-14", "avhrr2"),
        ("NOAA-7", "avhrr2"), ("NOAA-9", "avhrr2"), ("Metop-A", "avhrr3"),
        ("Metop-B", "avhrr3"), ("Metop-C", "avhrr3"), ("NOAA-15", "avhrr3"),
        ("NOAA-16", "avhrr3"), ("NOAA-17", "avhrr3"), ("NOAA-18", "avhrr3"),
        ("NOAA-19", "avhrr3"), ("HY-1C", "cocts"), ("Meteosat-12", "fci"),
        ("MTG-I1", "fci"), ("Metop-SG-A1", "metimage"), ("EOS-Aqua", "modis"),
        ("EOS-Terra", "modis"), ("Aqua", "modis"), ("Terra", "modis"),
        ("Sentinel-2A", "msi"), ("Sentinel-2B", "msi"), ("Sentinel-2C", "msi"),
        ("Arctica-M-N1", "msugsa"), ("Electro-L-N2", "msugs"),
        ("Sentinel-3A", "olci"), ("Sentinel-3B", "olci"),
        ("Landsat-8", "oli_tirs"), ("Landsat-9", "oli_tirs"),
        ("Meteosat-10", "seviri"), ("Meteosat-11", "seviri"),
        ("Meteosat-8", "seviri"), ("Meteosat-9", "seviri"),
        ("NOAA-20", "viirs"), ("NOAA-21", "viirs"), ("Suomi-NPP", "viirs"),
        ("FY-3D", "mersi2"), ("FY-3F", "mersi3"), ("FY-3G", "mersirm"),
        ("DSCOVR", "epic"),
    ]);

    let _instr = instruments.get(platform_name).copied().unwrap_or(sensor);
    sensor.replace('/', "").replace('-', "")
}

pub const AEROSOL_TYPES: &[&str] = &[
    "antarctic_aerosol", "continental_average_aerosol", "continental_clean_aerosol",
    "continental_polluted_aerosol", "desert_aerosol", "marine_clean_aerosol",
    "marine_polluted_aerosol", "marine_tropical_aerosol", "rayleigh_only",
    "rural_aerosol", "urban_aerosol",
];

pub const ATMOSPHERES: &[&str] = &[
    "subarctic_summer", "subarctic_winter", "midlatitude_summer",
    "midlatitude_winter", "tropical", "us_standard",
];
