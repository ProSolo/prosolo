use std::error::Error;

use clap;
use libprosic;
use libprosic::model::{AlleleFreq, ContinuousAlleleFreqs};
use rust_htslib::bam;
use bio::stats::Prob;


fn path_or_pipe(arg: Option<&str>) -> Option<&str> {
    arg.map_or(None, |f| if f == "-" { None } else { Some(f) })
}


pub fn single_cell_bulk(matches: &clap::ArgMatches) -> Result<(), Box<Error>> {
    // read command line parameters
    let single_mean_insert_size = value_t!(matches, "single-cell-insert-size-mean", f64).unwrap();
    let single_sd_insert_size = value_t!(matches, "single-cell-insert-size-sd", f64).unwrap();
    let bulk_mean_insert_size = value_t!(matches, "bulk-insert-size-mean", f64).unwrap_or(single_mean_insert_size);
    let bulk_sd_insert_size = value_t!(matches, "bulk-insert-size-sd", f64).unwrap_or(single_sd_insert_size);
//    let bulk_heterozygosity = try!(Prob::checked(value_t!(matches, "heterozygosity", f64).unwrap_or(1.25E-4)));
    let ploidy = value_t!(matches, "ploidy", u32).unwrap_or(2);
    let n_bulk_min = value_t!(matches, "bulk-min-n", usize).unwrap_or(8);
    let n_bulk_max = value_t!(matches, "bulk-max-n", usize).unwrap_or(100);
//    let bulk_effective_mutation_rate = value_t!(matches, "effective-mutation-rate", f64).unwrap_or(2000.0);
//    let deletion_factor = value_t!(matches, "deletion-factor", f64).unwrap_or(0.03);
//    let insertion_factor = value_t!(matches, "insertion-factor", f64).unwrap_or(0.01);
    let pileup_window = value_t!(matches, "pileup-window", u32).unwrap_or(2500);
    let no_fragment_evidence = matches.is_present("omit-fragment-evidence");
    let no_secondary = matches.is_present("omit-secondary-alignments");
    let no_mapq = matches.is_present("omit-mapq");
    let adjust_mapq = matches.is_present("adjust-mapq");
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
    let max_indel_overlap = value_t!(matches, "max-indel-overlap", u32).unwrap_or(25);

    let max_indel_len = value_t!(matches, "max-indel-len", u32).unwrap_or(1000);

    let single_bam = bam::IndexedReader::from_path(&single)?;
    let bulk_bam = bam::IndexedReader::from_path(&bulk)?;
/*
    let genome_size = (0..bulk_bam.header.target_count()).fold(0, |s, tid| {
        s + bulk_bam.header.target_len(tid).unwrap() as u64
    });
*/

    // dummy values until we implement indel calling for single cell data
    let prob_insertion_artifact = Prob(0.0);
    let prob_deletion_artifact = Prob(0.0);
    let prob_insertion_extend_artifact = Prob(0.0);
    let prob_deletion_extend_artifact = Prob(0.0);
    let indel_haplotype_window = 50;

    // init bulk sample
    let bulk_sample = libprosic::Sample::new(
        bulk_bam,
        pileup_window,
        !no_fragment_evidence,
        !no_secondary,
        !no_mapq,
        adjust_mapq,
        libprosic::InsertSize {
            mean: bulk_mean_insert_size,
            sd: bulk_sd_insert_size
        },
        libprosic::likelihood::LatentVariableModel::with_single_sample(),
        prob_insertion_artifact,
        prob_deletion_artifact,
        prob_insertion_extend_artifact,
        prob_deletion_extend_artifact,
        max_indel_overlap,
        indel_haplotype_window
    );

    // init single sample
    let single_sample = libprosic::Sample::new(
        single_bam,
        pileup_window,
        !no_fragment_evidence,
        !no_secondary,
        !no_mapq,
        adjust_mapq,
        libprosic::InsertSize {
            mean: single_mean_insert_size,
            sd: single_sd_insert_size
        },
        libprosic::likelihood::LatentVariableModel::with_single_sample(),
        prob_insertion_artifact,
        prob_deletion_artifact,
        prob_insertion_extend_artifact,
        prob_deletion_extend_artifact,
        max_indel_overlap,
        indel_haplotype_window
    );

    // setup events: case = single cell; control = bulk
    let events = [
        libprosic::call::pairwise::PairEvent {
            name: "hom_ref".to_owned(),
            af_case: vec![AlleleFreq(0.0)],
            af_control: ContinuousAlleleFreqs::right_exclusive( 0.0..0.5 )
        },
        libprosic::call::pairwise::PairEvent {
            name: "ADO_to_ref".to_owned(),
            af_case: vec![AlleleFreq(0.0)],
            af_control: ContinuousAlleleFreqs::inclusive( 0.5..1.0 )
        },
        libprosic::call::pairwise::PairEvent {
            name: "ADO_to_alt".to_owned(),
            af_case: vec![AlleleFreq(1.0)],
            af_control: ContinuousAlleleFreqs::inclusive( 0.0..0.5 )
        },
        libprosic::call::pairwise::PairEvent {
            name: "hom_alt".to_owned(),
            af_case: vec![AlleleFreq(1.0)],
            af_control: ContinuousAlleleFreqs::left_exclusive( 0.5..1.0 )
        },
        libprosic::call::pairwise::PairEvent {
            name: "err_alt".to_owned(),
            af_case: vec![AlleleFreq(0.5)],
            af_control: ContinuousAlleleFreqs::inclusive( 0.0..0.0 )
        },
        libprosic::call::pairwise::PairEvent {
            name: "het".to_owned(),
            af_case: vec![AlleleFreq(0.5)],
            af_control: ContinuousAlleleFreqs::exclusive( 0.0..1.0 )
        },
        libprosic::call::pairwise::PairEvent {
            name: "err_ref".to_owned(),
            af_case: vec![AlleleFreq(0.5)],
            af_control: ContinuousAlleleFreqs::inclusive( 1.0..1.0 )
        }
    ];

    let prior_model = libprosic::priors::SingleCellBulkModel::new(ploidy, n_bulk_min, n_bulk_max);

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
