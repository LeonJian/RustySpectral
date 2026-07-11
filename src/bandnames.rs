use std::collections::HashMap;

pub type BandNames = HashMap<&'static str, &'static str>;

pub fn get_bandnames() -> HashMap<&'static str, BandNames> {
    let mut names: HashMap<&'static str, BandNames> = HashMap::new();
    names.insert("generic", generic_bandnames());
    names.insert("modis", modis_bandnames());
    names.insert("seviri", seviri_bandnames());
    names.insert("viirs", viirs_bandnames());
    names.insert("avhrr3", avhrr3_bandnames());
    names.insert("abi", abi_bandnames());
    names.insert("agri", agri_bandnames());
    names.insert("ahi", ahi_bandnames());
    names.insert("ami", ami_bandnames());
    names.insert("fci", fci_bandnames());
    names.insert("slstr", slstr_bandnames());
    names.insert("vii", vii_bandnames());
    names
}

fn generic_bandnames() -> BandNames {
    let mut m = HashMap::from([
        ("VIS006", "VIS0.6"),
        ("VIS008", "VIS0.8"),
        ("IR_016", "NIR1.6"),
        ("IR_039", "IR3.9"),
        ("WV_062", "IR6.2"),
        ("WV_073", "IR7.3"),
        ("IR_087", "IR8.7"),
        ("IR_097", "IR9.7"),
        ("IR_108", "IR10.8"),
        ("IR_120", "IR12.0"),
        ("IR_134", "IR13.4"),
        ("HRV", "HRV"),
        ("I01", "I1"),
        ("I02", "I2"),
        ("I03", "I3"),
        ("I04", "I4"),
        ("I05", "I5"),
        ("M01", "M1"),
        ("M02", "M2"),
        ("M03", "M3"),
        ("M04", "M4"),
        ("M05", "M5"),
        ("M06", "M6"),
        ("M07", "M7"),
        ("M08", "M8"),
        ("M09", "M9"),
        ("C01", "ch1"),
        ("C02", "ch2"),
        ("C03", "ch3"),
        ("C04", "ch4"),
        ("C05", "ch5"),
        ("C06", "ch6"),
        ("C07", "ch7"),
        ("C08", "ch8"),
        ("C09", "ch9"),
        ("C10", "ch10"),
        ("C11", "ch11"),
        ("C12", "ch12"),
        ("C13", "ch13"),
        ("C14", "ch14"),
        ("C15", "ch15"),
        ("C16", "ch16"),
    ]);
    let nums: Vec<String> = (1..=36).map(|n| n.to_string()).collect();
    let labels: Vec<String> = (1..=36).map(|n| format!("ch{n}")).collect();
    for (num, label) in nums.iter().zip(labels.iter()) {
        // Leak to get &str with 'static lifetime
        let k: &'static str = Box::leak(num.clone().into_boxed_str());
        let v: &'static str = Box::leak(label.clone().into_boxed_str());
        m.insert(k, v);
    }
    m
}

fn modis_bandnames() -> BandNames {
    let mut m = HashMap::new();
    let nums: Vec<String> = (1..=36).map(|n| n.to_string()).collect();
    for n in &nums {
        let s: &'static str = Box::leak(n.clone().into_boxed_str());
        let s2: &'static str = Box::leak(n.clone().into_boxed_str());
        m.insert(s, s2);
    }
    m
}

fn seviri_bandnames() -> BandNames {
    HashMap::from([
        ("VIS006", "VIS0.6"),
        ("VIS008", "VIS0.8"),
        ("IR_016", "NIR1.6"),
        ("IR_039", "IR3.9"),
        ("WV_062", "IR6.2"),
        ("WV_073", "IR7.3"),
        ("IR_087", "IR8.7"),
        ("IR_097", "IR9.7"),
        ("IR_108", "IR10.8"),
        ("IR_120", "IR12.0"),
        ("IR_134", "IR13.4"),
        ("HRV", "HRV"),
    ])
}

fn viirs_bandnames() -> BandNames {
    HashMap::from([
        ("I01", "I1"),
        ("I02", "I2"),
        ("I03", "I3"),
        ("I04", "I4"),
        ("I05", "I5"),
        ("M01", "M1"),
        ("M02", "M2"),
        ("M03", "M3"),
        ("M04", "M4"),
        ("M05", "M5"),
        ("M06", "M6"),
        ("M07", "M7"),
        ("M08", "M8"),
        ("M09", "M9"),
    ])
}

fn avhrr3_bandnames() -> BandNames {
    HashMap::from([
        ("1", "ch1"),
        ("2", "ch2"),
        ("3b", "ch3b"),
        ("3a", "ch3a"),
        ("4", "ch4"),
        ("5", "ch5"),
    ])
}

fn abi_bandnames() -> BandNames {
    HashMap::from([
        ("C01", "ch1"),
        ("C02", "ch2"),
        ("C03", "ch3"),
        ("C04", "ch4"),
        ("C05", "ch5"),
        ("C06", "ch6"),
        ("C07", "ch7"),
        ("C08", "ch8"),
        ("C09", "ch9"),
        ("C10", "ch10"),
        ("C11", "ch11"),
        ("C12", "ch12"),
        ("C13", "ch13"),
        ("C14", "ch14"),
        ("C15", "ch15"),
        ("C16", "ch16"),
    ])
}

fn agri_bandnames() -> BandNames {
    HashMap::from([
        ("C01", "ch1"),
        ("C02", "ch2"),
        ("C03", "ch3"),
        ("C04", "ch4"),
        ("C05", "ch5"),
        ("C06", "ch6"),
        ("C07", "ch7"),
        ("C08", "ch8"),
        ("C09", "ch9"),
        ("C10", "ch10"),
        ("C11", "ch11"),
        ("C12", "ch12"),
        ("C13", "ch13"),
        ("C14", "ch14"),
    ])
}

fn ahi_bandnames() -> BandNames {
    HashMap::from([
        ("B01", "ch1"),
        ("B02", "ch2"),
        ("B03", "ch3"),
        ("B04", "ch4"),
        ("B05", "ch5"),
        ("B06", "ch6"),
        ("B07", "ch7"),
        ("B08", "ch8"),
        ("B09", "ch9"),
        ("B10", "ch10"),
        ("B11", "ch11"),
        ("B12", "ch12"),
        ("B13", "ch13"),
        ("B14", "ch14"),
        ("B15", "ch15"),
        ("B16", "ch16"),
    ])
}

fn ami_bandnames() -> BandNames {
    HashMap::from([
        ("VI004", "ch1"),
        ("VI005", "ch2"),
        ("VI006", "ch3"),
        ("VI008", "ch4"),
        ("NR013", "ch5"),
        ("NR016", "ch6"),
        ("SW038", "ch7"),
        ("WV063", "ch8"),
        ("WV069", "ch9"),
        ("WV073", "ch10"),
        ("IR087", "ch11"),
        ("IR096", "ch12"),
        ("IR105", "ch13"),
        ("IR112", "ch14"),
        ("IR123", "ch15"),
        ("IR133", "ch16"),
    ])
}

fn fci_bandnames() -> BandNames {
    HashMap::from([
        ("vis_04", "VIS0.4"),
        ("vis_05", "VIS0.5"),
        ("vis_06", "VIS0.6_HR"),
        ("vis_08", "VIS0.8"),
        ("vis_09", "VIS0.9"),
        ("nir_13", "NIR1.3"),
        ("nir_16", "NIR1.6"),
        ("nir_22", "NIR2.2_HR"),
        ("ir_38", "IR3.8_HR"),
        ("wv_63", "WV6.3"),
        ("wv_73", "WV7.3"),
        ("ir_87", "IR8.7"),
        ("ir_97", "IR9.7"),
        ("ir_105", "IR10.5_HR"),
        ("ir_123", "IR12.3"),
        ("ir_133", "IR13.3"),
    ])
}

fn slstr_bandnames() -> BandNames {
    HashMap::from([
        ("S1", "ch1"),
        ("S2", "ch2"),
        ("S3", "ch3"),
        ("S4", "ch4"),
        ("S5", "ch5"),
        ("S6", "ch6"),
        ("S7", "ch7"),
        ("S8", "ch8"),
        ("S9", "ch9"),
        ("F1", "ch7"),
        ("F2", "ch8"),
    ])
}

fn vii_bandnames() -> BandNames {
    HashMap::from([
        ("vii_443", "vii-4"),
        ("vii_555", "vii-8"),
        ("vii_668", "vii-12"),
        ("vii_752", "vii-15"),
        ("vii_763", "vii-16"),
        ("vii_865", "vii-17"),
        ("vii_914", "vii-20"),
        ("vii_1240", "vii-22"),
        ("vii_1375", "vii-23"),
        ("vii_1630", "vii-24"),
        ("vii_2250", "vii-25"),
        ("vii_3740", "vii-26"),
        ("vii_3959", "vii-28"),
        ("vii_4050", "vii-30"),
        ("vii_6725", "vii-33"),
        ("vii_7325", "vii-34"),
        ("vii_8540", "vii-35"),
        ("vii_10690", "vii-37"),
        ("vii_12020", "vii-39"),
        ("vii_13345", "vii-40"),
    ])
}
