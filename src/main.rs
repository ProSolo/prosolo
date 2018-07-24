#![cfg_attr(feature="flame_it", feature(plugin, custom_attribute))]
#![cfg_attr(feature="flame_it", plugin(flamer))]

#[cfg(feature="flame_it")]
extern crate flame;

#[macro_use]
extern crate log;
extern crate fern;
#[macro_use]
extern crate clap;
extern crate csv;
extern crate itertools;
extern crate time;
#[macro_use]
extern crate serde_json;

extern crate libprosic;
extern crate rust_htslib;
extern crate bio;

#[cfg(feature="flame_it")]
use std::fs::File;
use std::process;
use std::mem;

use clap::App;

pub mod call;
pub mod estimate;

/// Implementation as given by @ebarnard at: https://github.com/TyOverby/flame/issues/33#issuecomment-352312506
#[cfg(feature="flame_it")]
fn write_flamegraph() {
    let mut spans = flame::threads().into_iter().next().unwrap().spans;
    merge_spans(&mut spans);
    flame::dump_html_custom(&mut File::create("flame-graph.html").unwrap(), &spans)
        .unwrap();
}

#[cfg(feature="flame_it")]
fn merge_spans(spans: &mut Vec<flame::Span>) {
    if spans.is_empty() {
        return;
    }

    // Sort so spans to be merged are adjacent and spans with the most children are
    // merged into to minimise allocations.
    spans.sort_unstable_by(|s1, s2| {
        let a = (&s1.name, s1.depth, usize::max_value() - s1.children.len());
        let b = (&s2.name, s2.depth, usize::max_value() - s2.children.len());
        a.cmp(&b)
    });

    // Copy children and sum delta from spans to be merged
    let mut merge_targets = vec![0];
    {
        let mut spans_iter = spans.iter_mut().enumerate();
        let (_, mut current) = spans_iter.next().unwrap();
        for (i, span) in spans_iter {
            if current.name == span.name && current.depth == span.depth {
                current.delta += span.delta;
                let mut children = mem::replace(&mut span.children, Vec::new());
                current.children.extend(children.into_iter());
            } else {
                current = span;
                merge_targets.push(i);
            }
        }
    }

    // Move merged spans to the front of the spans vector
    for (target_i, &current_i) in merge_targets.iter().enumerate() {
        spans.swap(target_i, current_i);
    }

    // Remove duplicate spans
    spans.truncate(merge_targets.len());

    // Merge children of the newly collapsed spans
    for span in spans {
        merge_spans(&mut span.children);
    }
}


#[cfg(not(feature="flame_it"))]
fn write_flamegraph() {
}


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
    } else if let Some(matches) = matches.subcommand_matches("apply-fdr") {
        if let Err(e) = estimate::apply_fdr(matches) {
            error!("Error: {}", e);
            process::exit(1);
        }
    }

    write_flamegraph();
}
