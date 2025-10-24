//! Reads Configuration parameters form a specified location, or returns a default one if none is
//! found.

use std::path::PathBuf;

use config_file::FromConfigFile;
use serde::Deserialize;

/// Make it discoverable for both tests and binaries, even if they are run from a deeper directory.
const PATH: [&str; 3] = ["./config.toml", "../config.toml", "../../config.toml"];

#[derive(Deserialize, Debug, Clone)]
/// See `config.toml` for descriptions
pub struct Config {
    pub max_steps: usize,
    pub rkf45_first_step: f64,
    pub energy_rel_tol: f64,
    pub step_rel_tol: f64,
    pub evolution_init_capacity: usize,
    pub map_threshold: f64,
}

/// Default configuration if no file is found
impl Default for Config {
    fn default() -> Self {
        Self {
            max_steps: 1e7 as usize,
            rkf45_first_step: 1e-3,
            energy_rel_tol: 1e-12,
            step_rel_tol: 1e-12,
            evolution_init_capacity: 2000,
            map_threshold: 1e-9,
        }
    }
}

pub fn get_config() -> Config {
    for path_str in PATH {
        let path = PathBuf::from(path_str).canonicalize().unwrap_or_default();

        if path.exists() {
            match Config::from_config_file(&path) {
                // TODO: find a way to open the config file only once, so the messages are printed
                Ok(file) => {
                    // eprintln!("Reading configuration from {}", path.display());
                    return file;
                }
                Err(err) => {
                    eprintln!("Ignoring invalid config file '{}': {}", path.display(), err);
                }
            }
        }
    }
    eprintln!("Could not find config.toml, using default configuration.");
    Config::default()
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_config_file() {
        let config_file = get_config();
        dbg!(&config_file);

        assert_eq!(config_file.max_steps, 1e7 as usize)
    }
}
