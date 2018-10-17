#[macro_use]
extern crate log;
extern crate fern;
#[macro_use]
extern crate clap;
extern crate csv;
extern crate itertools;
extern crate time;
extern crate serde_json;

extern crate libprosic;
extern crate rust_htslib;
extern crate bio;

use std::process;

use clap::App;

pub mod call;
pub mod estimate;

fn main() {
    // parse command line
    let yaml = load_yaml!("cli.yaml");
    let matches = App::from_yaml(yaml)
                      .version(env!("CARGO_PKG_VERSION"))
                      .get_matches();

    // setup logger
    let logger_config = fern::DispatchConfig {
        format: Box::new(|msg: &str, level: &log::LogLevel, _: &log::LogLocation| {
          match level {
              &log::LogLevel::Debug => format!("DEBUG[{}]: {}", time::now().strftime("%H:%M:%S").unwrap(), msg),
              _ => msg.to_owned()
          }
        }),
        output: vec![fern::OutputConfig::stderr()],
        level: log::LogLevelFilter::Trace,
    };
    if let Err(e) = fern::init_global_logger(
        logger_config,
        if matches.is_present("verbose") { log::LogLevelFilter::Debug } else { log::LogLevelFilter::Info }
    ) {
        panic!("Failed to initialize logger: {}", e);
    }


    if let Some(matches) = matches.subcommand_matches("single-cell-bulk") {
        if let Err(e) = call::single_cell_bulk(matches) {
            error!("Error: {}", e);
            process::exit(1);
        }
    } else if let Some(matches) = matches.subcommand_matches("estimate-mutation-rate") {
        if let Err(e) = estimate::effective_mutation_rate(matches) {
            error!("Error: {}", e);
            process::exit(1);
        }
    } else if let Some(matches) = matches.subcommand_matches("control-fdr") {
        if let Err(e) = estimate::fdr(matches) {
            error!("Error: {}", e);
            process::exit(1);
        }
    }

}
