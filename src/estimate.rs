use std::error::Error;
use std::fs::File;
use std::io;

use bio::stats::{LogProb, Prob};
use clap;
use csv;
use serde_json;

use libprosic;
use libprosic::estimation;
use libprosic::model;
use libprosic::model::AlleleFreq;

use call;

pub fn effective_mutation_rate(matches: &clap::ArgMatches) -> Result<(), Box<Error>> {
    let min_af = value_t!(matches, "min-af", f64).unwrap_or(0.12);
    let max_af = value_t!(matches, "max-af", f64).unwrap_or(0.25);
    let mut reader = csv::Reader::from_reader(io::stdin());
    let freqs = try!(reader.deserialize().collect::<Result<Vec<f64>, _>>());
    let estimate = estimation::effective_mutation_rate::estimate(
        freqs
            .into_iter()
            .filter(|&f| f >= min_af && f <= max_af)
            .map(|f| AlleleFreq(f)),
    );

    // print estimated mutation rate to stdout
    println!("{}", estimate.effective_mutation_rate());

    // if --fit is given, print data visualizing model fit
    if let Some(path) = matches.value_of("fit") {
        let mut f = try!(File::create(path));
        serde_json::to_writer(&mut f, &estimate)?;
    }
    Ok(())
}

struct DummyEvent {
    pub name: String,
}

impl libprosic::Event for DummyEvent {
    fn name(&self) -> &str {
        &self.name
    }
}

/// Parse `VariantType` from command line arguments.
pub fn parse_vartype(
    vartype: &str,
    min_len: Option<u32>,
    max_len: Option<u32>,
) -> Result<model::VariantType, Box<Error>> {
    Ok(match (vartype, min_len, max_len) {
        ("SNV", _, _) => model::VariantType::SNV,
        ("INS", Some(min_len), Some(max_len)) => {
            model::VariantType::Insertion(Some(min_len..max_len))
        }
        ("DEL", Some(min_len), Some(max_len)) => {
            model::VariantType::Deletion(Some(min_len..max_len))
        }
        ("INS", _, _) => model::VariantType::Insertion(None),
        ("DEL", _, _) => model::VariantType::Deletion(None),
        _ => {
            return Err(Box::new(clap::Error {
                message: "unsupported variant type (supported: SNV, INS, DEL)".to_owned(),
                kind: clap::ErrorKind::InvalidValue,
                info: None,
            }));
        }
    })
}

pub fn fdr(matches: &clap::ArgMatches) -> Result<(), Box<Error>> {
    let call_bcf = matches.value_of("calls").unwrap();
    let alpha = value_t!(matches, "alpha", f64).unwrap();
    let events_list = matches.values_of("events").unwrap();
    let vartype = matches.value_of("vartype").unwrap();
    let min_len = value_t!(matches, "min-len", u32).ok();
    let max_len = value_t!(matches, "max-len", u32).ok();
    let vartype = parse_vartype(vartype, min_len, max_len)?;

    let events: Vec<DummyEvent> = events_list
        .map(|ev| DummyEvent {
            name: ev.to_owned(),
        }).collect();
    let alpha = LogProb::from(Prob::checked(alpha)?);

    let out = call::path_or_pipe(matches.value_of("output"));

    estimation::fdr::ev::control_fdr::<_, _, &str>(call_bcf, out, &events, &vartype, alpha)?;

    Ok(())
}
