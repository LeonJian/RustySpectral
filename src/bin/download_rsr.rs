use std::path::Path;
use std::process;

fn main() {
    env_logger::init();

    let mut dest_dir: Option<&Path> = None;
    let mut dry_run = false;
    let args: Vec<String> = std::env::args().collect();
    let mut i = 1;

    while i < args.len() {
        match args[i].as_str() {
            "-o" | "--destination" => {
                i += 1;
                if i < args.len() {
                    dest_dir = Some(Path::new(&args[i]));
                }
            }
            "-d" | "--dry_run" => dry_run = true,
            "-v" | "--verbose" => {
                // already handled by env_logger
            }
            _ => {}
        }
        i += 1;
    }

    match rustyspectral::download::download_rsr(dest_dir, dry_run) {
        Ok(()) => println!("RSR data downloaded successfully."),
        Err(e) => {
            eprintln!("Error: {}", e);
            process::exit(1);
        }
    }
}
