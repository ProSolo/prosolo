use std::error::Error;

use clap;
use libprosic;
use libprosic::model::{AlleleFreq, ContinuousAlleleFreqs, DiscreteAlleleFreqs};
use rust_htslib::bam;
use bio::stats::Prob;

use std::collections::BTreeMap;

pub fn path_or_pipe(arg: Option<&str>) -> Option<&str> {
    arg.map_or(None, |f| if f == "-" { None } else { Some(f) })
}

fn alignment_properties(path: &str) -> Result<libprosic::AlignmentProperties, Box<Error>> {
    let mut bam = bam::Reader::from_path(path)?;
    Ok(libprosic::AlignmentProperties::estimate(&mut bam)?)
}

pub fn single_cell_bulk(matches: &clap::ArgMatches) -> Result<(), Box<Error>> {
    // read command line parameters
//    let bulk_heterozygosity = try!(Prob::checked(value_t!(matches, "heterozygosity", f64).unwrap_or(1.25E-4)));
    let ploidy = value_t!(matches, "ploidy", u32).unwrap();
    let n_bulk_min = value_t!(matches, "bulk-min-n", usize).unwrap();
    let n_bulk_max = value_t!(matches, "bulk-max-n", usize).unwrap();
//    let single_effective_mutation_rate = value_t!(matches, "effective-mutation-rate", f64).unwrap_or(2000.0);
//    let deletion_factor = value_t!(matches, "deletion-factor", f64).unwrap();
//    let insertion_factor = value_t!(matches, "insertion-factor", f64).unwrap();
    let pileup_window = value_t!(matches, "pileup-window", u32).unwrap();
    let no_fragment_evidence = matches.is_present("omit-fragment-evidence");
    let omit_snvs = matches.is_present("omit-snvs");
    let omit_indels = matches.is_present("omit-indels");
    let single = matches.value_of("single-cell").unwrap();
    let bulk = matches.value_of("bulk").unwrap();
    let candidates = path_or_pipe(matches.value_of("candidates"));
    let output = path_or_pipe(matches.value_of("output"));
    let reference = matches.value_of("reference").unwrap();
    let observations = matches.value_of("observations");
//    let flat_priors = matches.is_present("flat-priors");
    let exclusive_end = matches.is_present("exclusive-end");
    let indel_haplotype_window = value_t!(matches, "indel-window", u32).unwrap();

    let prob_spurious_ins = Prob::checked(value_t_or_exit!(matches, "prob-spurious-ins", f64))?;
    let prob_spurious_del = Prob::checked(value_t_or_exit!(matches, "prob-spurious-del", f64))?;
    let prob_ins_extend = Prob::checked(value_t_or_exit!(matches, "prob-ins-extend", f64))?;
    let prob_del_extend = Prob::checked(value_t_or_exit!(matches, "prob-del-extend", f64))?;

    let max_indel_len = value_t!(matches, "max-indel-len", u32).unwrap_or(1000);

    let single_alignment_properties = alignment_properties(&single)?;
    let bulk_alignment_properties = alignment_properties(&bulk)?;

    info!("estimated single properties: {:?}", single_alignment_properties);
    info!("estimated bulk properties: {:?}", bulk_alignment_properties);

    let single_bam = bam::IndexedReader::from_path(&single)?;
    let bulk_bam = bam::IndexedReader::from_path(&bulk)?;

    // init bulk sample
    let bulk_sample = libprosic::Sample::new(
        bulk_bam,
        pileup_window,
        !no_fragment_evidence,
        false,
        false,
        bulk_alignment_properties,
        libprosic::likelihood::LatentVariableModel::with_single_sample(),
        prob_spurious_ins,
        prob_spurious_del,
        prob_ins_extend,
        prob_del_extend,
        indel_haplotype_window
    );

    // init single sample
    let single_sample = libprosic::Sample::new(
        single_bam,
        pileup_window,
        !no_fragment_evidence,
        false,
        false,
        single_alignment_properties,
        libprosic::likelihood::LatentVariableModel::with_single_sample(),
        prob_spurious_ins,
        prob_spurious_del,
        prob_ins_extend,
        prob_del_extend,
        indel_haplotype_window
    );

    // setup events: case = single cell; control = bulk
    let events = [
        libprosic::call::pairwise::PairEvent {
            name: "hom_ref".to_owned(),
            af_case: DiscreteAlleleFreqs::absent(),
            af_control: ContinuousAlleleFreqs::right_exclusive( 0.0..0.5 )
        },
        libprosic::call::pairwise::PairEvent {
            name: "ADO_to_ref".to_owned(),
            af_case: DiscreteAlleleFreqs::absent(),
            af_control: ContinuousAlleleFreqs::right_exclusive( 0.5..1.0 )
        },
        libprosic::call::pairwise::PairEvent {
            name: "ADO_to_alt".to_owned(),
            af_case: DiscreteAlleleFreqs::new( vec![AlleleFreq(1.0)] ),
            af_control: ContinuousAlleleFreqs::left_exclusive( 0.0..0.5 )
        },
        libprosic::call::pairwise::PairEvent {
            name: "hom_alt".to_owned(),
            af_case: DiscreteAlleleFreqs::new( vec![AlleleFreq(1.0)] ),
            af_control: ContinuousAlleleFreqs::left_exclusive( 0.5..1.0 )
        },
        libprosic::call::pairwise::PairEvent {
            name: "err_alt".to_owned(),
            af_case: DiscreteAlleleFreqs::feasible(2).not_absent(),
            af_control: ContinuousAlleleFreqs::inclusive( 0.0..0.0 )
        },
        libprosic::call::pairwise::PairEvent {
            name: "het".to_owned(),
            af_case: DiscreteAlleleFreqs::new( vec![AlleleFreq(0.5)] ),
            af_control: ContinuousAlleleFreqs::exclusive( 0.0..1.0 )
        },
        libprosic::call::pairwise::PairEvent {
            name: "err_ref".to_owned(),
            af_case: DiscreteAlleleFreqs::new( vec![AlleleFreq(0.0), AlleleFreq(0.5)] ),
            af_control: ContinuousAlleleFreqs::inclusive( 1.0..1.0 )
        }
    ];

    // test that the `Event`s exactly cover the `Event` space
    {
        let mut event_collection = BTreeMap::new();
        for daf in DiscreteAlleleFreqs::feasible(2).iter() {
            event_collection.insert(*daf, BTreeMap::new());
        }
        for event in events.iter() {
            for daf in event.af_case.iter() {
                let mut start = event.af_control.start;
                // make sure events get added to the BTreeMap in the correct order, putting left_exclusive after inclusive point ContinuousAlleleFreqs
                if event.af_control.left_exclusive { start += 0.0000000001};
                event_collection.get_mut(&daf).unwrap().insert(start, &event.af_control);
            }
        }
        for (daf, cafs) in event_collection.iter() {
            let mut first_caf = cafs.iter();
            if let Some( (_, first) ) = first_caf.nth(0) {
                if first.left_exclusive | ( first.start != AlleleFreq(0.0) ) {
                    panic!("At DiscreteAlleleFreqs '{}', the lowest ContinuousAlleleFreqs starts at '{}' and left_exclusive is '{}'. Thus the Event space is not fully defined.", daf, first.start, first.left_exclusive);
                }
            } else {
                panic!("At DiscreteAlleleFreqs '{}', there are no ContinuousAlleleFreqs defined in Events, thus the Event space is not fully defined.", daf);
            }
            let mut peekable = cafs.into_iter().peekable();
            while let Some( (_, current) ) = peekable.next() {
                if let Some( (_, next) ) = peekable.peek() {
                    if current.end == next.start {
                        if current.right_exclusive ^ next.left_exclusive {
                            continue;
                        } else if current.right_exclusive & next.left_exclusive {
                            panic!("At DiscreteAlleleFreqs '{}', at {} the end of the previous ContinuousAlleleFreqs and the start of the ContinuousAlleleFreqs are exclusive. Thus this point in the Event space is not defined.", daf, current.end);
                        } else {
                            panic!("At DiscreteAlleleFreqs '{}', at {} the end of the previous ContinuousAlleleFreqs and the start of the ContinuousAlleleFreqs are inclusive. Thus this point in the Event space is multiply defined.", daf, current.end);
                        }
                    } else if current.end < next.start {
                        panic!("At DiscreteAlleleFreqs '{}', there is a gap between the end of one ContinuousAlleleFreqs at '{}' and the next starting at '{}'. Thus the Event space is not fully defined.", daf, current.end, next.start);
                    } else if current.end > next.start {
                        panic!("At DiscreteAlleleFreqs '{}', there is an overlap between the end of one ContinuousAlleleFreqs at '{}' and the next starting at '{}'. Thus, this part of the Event space is multiply defined.", daf, current.end, next.start);
                    }
                } else {
                    // last item on the iterator
                    if current.right_exclusive | ( current.end != AlleleFreq(1.0) ) {
                        panic!("At DiscreteAlleleFreqs '{}', the highest ContinuousAlleleFreqs ends at '{}' and right_exclusive is '{}'. Thus the Event space is not fully defined.", daf, current.end, current.right_exclusive);
                    }
                }
            }
        }
    }

    let prior_model = libprosic::priors::SingleCellBulkModel::new(
        ploidy,
        n_bulk_min,
        n_bulk_max
    );

    // init joint model
    let mut joint_model = libprosic::model::PairCaller::new(
        single_sample,
        bulk_sample,
        prior_model
    );

    // perform calling
    libprosic::call::pairwise::call::<
        _, _, _,
        libprosic::model::PairCaller<
            libprosic::model::DiscreteAlleleFreqs,
            libprosic::model::ContinuousAlleleFreqs,
            libprosic::model::priors::SingleCellBulkModel
        >, _, _, _, _>
    (
        candidates,
        output,
        &reference,
        &events,
        &mut joint_model,
        omit_snvs,
        omit_indels,
        Some(max_indel_len),
        observations.as_ref(),
        exclusive_end
    )
}
