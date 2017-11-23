use std::io;
use std::error::Error;
use std::io::prelude::*;
use std::fs::File;

use rustc_serialize::json;
use clap;
use csv;
use rust_htslib::bcf;

use libprosic;
use libprosic::model::AlleleFreq;
use libprosic::estimation;
use libprosic::model;
use libprosic::utils;

use call;

pub fn effective_mutation_rate(matches: &clap::ArgMatches) -> Result<(), Box<Error>> {
    let min_af = value_t!(matches, "min-af", f64).unwrap_or(0.12);
    let max_af = value_t!(matches, "max-af", f64).unwrap_or(0.25);
    let mut reader = csv::Reader::from_reader(io::stdin());
    let freqs = try!(reader.decode().collect::<Result<Vec<f64>, _>>());
    let estimate = estimation::effective_mutation_rate::estimate(freqs.into_iter().filter(|&f| {
        f >= min_af && f <= max_af
    }).map(|f| AlleleFreq(f)));

    // print estimated mutation rate to stdout
    println!("{}", estimate.effective_mutation_rate());

    // if --fit is given, print data visualizing model fit
    if let Some(path) = matches.value_of("fit") {
        let json = json::encode(&estimate).unwrap();
        let mut f = try!(File::create(path));
        try!(f.write_all(json.as_bytes()));
    }
    Ok(())
}


struct DummyEvent {
    pub name: String
}


impl libprosic::Event for DummyEvent {
    fn name(&self) -> &str {
        &self.name
    }
}


/// Parse `VariantType` from command line arguments.
pub fn parse_vartype(vartype: &str, min_len: Option<u32>, max_len: Option<u32>) -> Result<model::VariantType, Box<Error>> {
    Ok(match (vartype, min_len, max_len) {
        ("SNV", _, _) => model::VariantType::SNV,
        ("INS", Some(min_len), Some(max_len)) => model::VariantType::Insertion(Some(min_len..max_len)),
        ("DEL", Some(min_len), Some(max_len)) => model::VariantType::Deletion(Some(min_len..max_len)),
        ("INS", _, _) => model::VariantType::Insertion(None),
        ("DEL", _, _) => model::VariantType::Deletion(None),
        _ => {
            return Err(Box::new(clap::Error {
                message: "unsupported variant type (supported: SNV, INS, DEL)".to_owned(),
                kind: clap::ErrorKind::InvalidValue,
                info: None
            }));
        }
    })
}


pub fn fdr(matches: &clap::ArgMatches) -> Result<(), Box<Error>> {
    let call_bcf = matches.value_of("calls").unwrap();
    let events_list = matches.values_of("events").unwrap();
    let vartype = matches.value_of("vartype").unwrap();
    let min_len = value_t!(matches, "min-len", u32).ok();
    let max_len = value_t!(matches, "max-len", u32).ok();
    let vartype = parse_vartype(vartype, min_len, max_len)?;
    let method = matches.value_of("method").unwrap();

    let events: Vec<DummyEvent> = events_list.map(|ev| DummyEvent { name: ev.to_owned() }).collect();
    let mut writer = io::stdout();
    let mut call_reader = try!(bcf::Reader::from_path(&call_bcf));

    match method {
        "bh" => {
            let null_bcf = matches.value_of("null-calls").unwrap_or_else(|| panic!("--null-calls required for controlling fdr with bh method") );
            let mut null_reader = try!(bcf::Reader::from_path(&null_bcf));

            estimation::fdr::bh::control_fdr(
                &mut call_reader,
                &mut null_reader,
                &mut writer,
                &events,
                &vartype
            )?;
        },
        "ev" => {
            estimation::fdr::ev::control_fdr(
                &mut call_reader,
                &mut writer,
                &events,
                &vartype
            )?;
        },
        _ => panic!("Undefined method for controlling FDR.")
    }

    Ok(())
}

pub fn apply_fdr(matches: &clap::ArgMatches) -> Result<(), Box<Error>> {
    let call_bcf = matches.value_of("calls").unwrap();
    let mut call_reader = try!(bcf::Reader::from_path(&call_bcf));

    let threshold = value_t!(matches, "threshold", f64).unwrap();

    let events_list = matches.values_of("events").unwrap();
    let events: Vec<DummyEvent> = events_list.map(|ev| DummyEvent { name: ev.to_owned() }).collect();

    let vartype = matches.value_of("vartype").unwrap();
    let min_len = value_t!(matches, "min-len", u32).ok();
    let max_len = value_t!(matches, "max-len", u32).ok();
    let vartype = parse_vartype(vartype, min_len, max_len)?;

    let out = call::path_or_pipe(matches.value_of("output"));
    let header = bcf::Header::with_template(&call_reader.header);
    let mut writer = match out {
        Some(f) => bcf::Writer::from_path(f, &header, false, false)?,
        None    => bcf::Writer::from_stdout(&header, false, false)?
    };

    utils::filter_by_threshold(
        &mut call_reader,
        &threshold,
        &mut writer,
        &events,
        &vartype
    )?;

    Ok(())
}
