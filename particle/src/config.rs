use std::path::PathBuf;

use config_file::FromConfigFile;
use serde::Deserialize;

const PATH: &str = "../config.toml";

#[derive(Deserialize, Debug)]
pub struct Config {
    pub max_steps: usize,
    pub rkf45_first_step: f64,
    pub energy_rel_tol: f64,
    pub step_rel_tol: f64,
}

/// Default configuration if no file is found
impl Default for Config {
    fn default() -> Self {
        Self {
            max_steps: 1e7 as usize,
            rkf45_first_step: 1e-3,
            energy_rel_tol: 1e-12,
            step_rel_tol: 1e-12,
        }
    }
}

pub fn get_config() -> Config {
    let path = PathBuf::from(PATH).canonicalize().unwrap_or_default();

    match path.exists() {
        true => match Config::from_config_file(&path) {
            Ok(file) => {
                eprintln!("Reading configuration from {}", path.display());
                return file;
            }
            Err(err) => {
                eprintln!("Ignoring invalid config file '{}': {}", path.display(), err);
            }
        },
        false => {
            eprintln!("Config file {} doesn't exists", path.display());
        }
    }
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
