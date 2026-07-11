use ndarray::Array1;

pub fn trapezoid(y: &Array1<f64>, x: &Array1<f64>) -> f64 {
    let n = y.len();
    if n < 2 {
        return 0.0;
    }
    let mut sum = 0.0;
    for i in 1..n {
        sum += (y[i] + y[i - 1]) * (x[i] - x[i - 1]);
    }
    0.5 * sum
}

pub fn get_central_wave(wav: &Array1<f64>, resp: &Array1<f64>, weight: f64) -> f64 {
    let numerator = trapezoid(&(resp * wav * weight), wav);
    let denominator = trapezoid(&(resp * weight), wav);
    if denominator == 0.0 {
        return f64::NAN;
    }
    numerator / denominator
}

pub fn convert2wavenumber_rsr(wavelength: &Array1<f64>, response: &Array1<f64>) -> (Array1<f64>, Array1<f64>) {
    let n = wavelength.len();
    let mut wavenumber = Array1::zeros(n);
    let mut resp_out = Array1::zeros(n);

    for i in 0..n {
        let j = n - 1 - i;
        wavenumber[i] = 1.0 / (1e-4 * wavelength[j]);
        resp_out[i] = response[j];
    }

    (wavenumber, resp_out)
}

pub fn sort_data(x_vals: &Array1<f64>, y_vals: &Array1<f64>) -> (Array1<f64>, Array1<f64>) {
    let n = x_vals.len();

    let mut indices: Vec<usize> = (0..n).collect();
    indices.sort_by(|&a, &b| {
        x_vals[a].partial_cmp(&x_vals[b]).unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut sorted_x = Array1::zeros(n);
    let mut sorted_y = Array1::zeros(n);
    for (i, &idx) in indices.iter().enumerate() {
        sorted_x[i] = x_vals[idx];
        sorted_y[i] = y_vals[idx];
    }

    let mut keep = vec![true; n];
    for i in 1..n {
        if sorted_x[i] <= sorted_x[i - 1] {
            keep[i] = false;
        }
    }

    let deduped: Vec<_> = (0..n).filter(|&i| keep[i]).collect();
    let m = deduped.len();
    let mut result_x = Array1::zeros(m);
    let mut result_y = Array1::zeros(m);
    for (i, &idx) in deduped.iter().enumerate() {
        result_x[i] = sorted_x[idx];
        result_y[i] = sorted_y[idx];
    }

    (result_x, result_y)
}

pub fn get_fullwidth_halfmax(rsp: &Array1<f64>, wvl: &Array1<f64>) -> f64 {
    let half_max = 0.5;
    let indices: Vec<usize> = rsp
        .iter()
        .enumerate()
        .filter(|(_, &v)| v >= half_max)
        .map(|(i, _)| i)
        .collect();

    if indices.len() < 2 {
        return f64::NAN;
    }

    wvl[indices[indices.len() - 1]] - wvl[indices[0]]
}

pub fn get_bounds_integrated_energy(rsp: &Array1<f64>, wvl: &Array1<f64>, ener_perc_lim: f64) -> (f64, f64) {
    let n = rsp.len();
    let mut crs = Array1::zeros(n);
    crs[0] = rsp[0];
    for i in 1..n {
        crs[i] = crs[i - 1] + rsp[i];
    }
    let max_val = crs[n - 1];
    for i in 0..n {
        crs[i] = crs[i] / max_val * 100.0;
    }

    let low: Vec<usize> = crs
        .iter()
        .enumerate()
        .filter(|(_, &v)| v >= ener_perc_lim)
        .map(|(i, _)| i)
        .collect();

    let high: Vec<usize> = crs
        .iter()
        .enumerate()
        .filter(|(_, &v)| v <= (100.0 - ener_perc_lim))
        .map(|(i, _)| i)
        .collect();

    let min_wvl = wvl[low[0]];
    let max_wvl = wvl[high[high.len() - 1]];

    (min_wvl, max_wvl)
}

pub fn get_wave_range(wvl: &Array1<f64>, resp: &Array1<f64>, threshold: f64) -> (f64, f64, f64) {
    let cwl = get_central_wave(wvl, resp, 1.0);

    let pts: Vec<usize> = resp
        .iter()
        .enumerate()
        .filter(|(_, &v)| v > threshold)
        .map(|(i, _)| i)
        .collect();

    let min_wvl = wvl[pts[0]];
    let max_wvl = wvl[pts[pts.len() - 1]];

    (min_wvl, cwl, max_wvl)
}

pub fn are_instruments_identical(name1: &str, name2: &str) -> bool {
    if name1 == name2 {
        return true;
    }
    let translate = |s: &str| -> String {
        match s {
            "avhrr-1" => "avhrr/1".to_string(),
            "avhrr-2" => "avhrr/2".to_string(),
            "avhrr-3" => "avhrr/3".to_string(),
            _ => s.to_string(),
        }
    };
    translate(name1) == translate(name2)
}

pub fn check_and_adjust_instrument_name(_platform_name: &str, instrument: &str) -> String {
    instrument.to_lowercase().replace('/', "").replace('-', "")
}
