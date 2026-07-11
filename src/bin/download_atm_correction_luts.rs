use std::process;

fn main() {
    env_logger::init();

    let mut dry_run = false;
    let mut aerosol_types: Vec<String> = Vec::new();
    let args: Vec<String> = std::env::args().collect();
    let mut i = 1;

    while i < args.len() {
        match args[i].as_str() {
            "-a" | "--aerosol_types" => {
                i += 1;
                if i < args.len() {
                    aerosol_types = args[i].split(',').map(|s| s.to_string()).collect();
                }
            }
            "-d" | "--dry_run" => dry_run = true,
            "-v" | "--verbose" => {}
            _ => {}
        }
        i += 1;
    }

    let types = if aerosol_types.is_empty() {
        None
    } else {
        Some(&aerosol_types[..])
    };

    match rustyspectral::download::download_luts(types, dry_run) {
        Ok(()) => println!("Atm correction LUTs downloaded successfully."),
        Err(e) => {
            eprintln!("Error: {}", e);
            process::exit(1);
        }
    }
}
