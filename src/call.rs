use std::error::Error;

use clap;
use libprosic;
use libprosic::model::AlleleFreq;
use rust_htslib::bam;
use bio::stats::Prob;

pub fn single_cell_bulk(matches: &clap::ArgMatches) -> Result<(), Box<Error>> {
    // read command line parameters
    let single_mean_insert_size = value_t!(matches, "single-cell-insert-size-mean", f64).unwrap();
    let single_sd_insert_size = value_t!(matches, "single-cell-insert-size-sd", f64).unwrap();
    let bulk_mean_insert_size = value_t!(matches, "bulk-insert-size-mean", f64).unwrap_or(single_mean_insert_size);
    let bulk_sd_insert_size = value_t!(matches, "bulk-insert-size-sd", f64).unwrap_or(single_sd_insert_size);
    let bulk_heterozygosity = try!(Prob::checked(value_t!(matches, "heterozygosity", f64).unwrap_or(1.25E-4)));
    let ploidy = value_t!(matches, "ploidy", u32).unwrap_or(2);
    let bulk_effective_mutation_rate = value_t!(matches, "effective-mutation-rate", f64).unwrap_or(2000.0);
    let deletion_factor = value_t!(matches, "deletion-factor", f64).unwrap_or(0.03);
    let insertion_factor = value_t!(matches, "insertion-factor", f64).unwrap_or(0.01);
    let min_somatic_af = value_t!(matches, "min-somatic-af", f64).map(|af| AlleleFreq(af)).unwrap_or(AlleleFreq(0.05));
    let pileup_window = value_t!(matches, "pileup-window", u32).unwrap_or(2500);
    let no_fragment_evidence = matches.is_present("omit-fragment-evidence");
    let no_secondary = matches.is_present("omit-secondary-alignments");
    let no_mapq = matches.is_present("omit-mapq");
    let adjust_mapq = matches.is_present("adjust-mapq");
    let omit_snvs = matches.is_present("omit-snvs");
    let omit_indels = matches.is_present("omit-indels");
    let single = matches.value_of("single-cell").unwrap();
    let bulk = matches.value_of("bulk").unwrap();
    let candidates = matches.value_of("candidates").unwrap_or("-");
    let output = matches.value_of("output").unwrap_or("-");
    let reference = matches.value_of("reference").unwrap();
    let observations = matches.value_of("observations");
    let flat_priors = matches.is_present("flat-priors");
    let exclusive_end = matches.is_present("exclusive-end");

    let prob_spurious_isize = try!(Prob::checked(value_t!(matches, "prob-spurious-isize", f64).unwrap_or(0.0)));

    let max_indel_dist = value_t!(matches, "max-indel-dist", u32).unwrap_or(50);
    let max_indel_len_diff = value_t!(matches, "max-indel-len-diff", u32).unwrap_or(20);
    let max_indel_len = value_t!(matches, "max-indel-len", u32).unwrap_or(1000);

    let single_bam = try!(bam::IndexedReader::new(&single));
    let bulk_bam = try!(bam::IndexedReader::new(&bulk));
    let genome_size = (0..bulk_bam.header.target_count()).fold(0, |s, tid| {
        s + bulk_bam.header.target_len(tid).unwrap() as u64
    });

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
        prob_spurious_isize
    ).max_indel_dist(max_indel_dist).max_indel_len_diff(max_indel_len_diff);

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
        prob_spurious_isize
    ).max_indel_dist(max_indel_dist).max_indel_len_diff(max_indel_len_diff);

    // setup events: case = single cell; control = bulk
    let events = [
        libprosic::call::pairwise::PairEvent {
            name: "germline_hom_ref".to_owned(),
            af_case: vec![AlleleFreq(0.0)],
            af_control: AlleleFreq(0.0)..AlleleFreq(0.25000001)
        },
        libprosic::call::pairwise::PairEvent {
            name: "somatic_hom_ref".to_owned(),
            af_case: vec![AlleleFreq(0.0)],
            af_control: AlleleFreq(0.25000001)..AlleleFreq(0.5)
        },
        libprosic::call::pairwise::PairEvent {
            name: "alt_allele_dropout".to_owned(),
            af_case: vec![AlleleFreq(0.0)],
            af_control: AlleleFreq(0.5)..AlleleFreq(1.00000001)
        },
        libprosic::call::pairwise::PairEvent {
            name: "ref_allele_dropout".to_owned(),
            af_case: vec![AlleleFreq(1.0)],
            af_control: AlleleFreq(0.0)..AlleleFreq(0.50000001)
        },
        libprosic::call::pairwise::PairEvent {
            name: "somatic_hom_alt".to_owned(),
            af_case: vec![AlleleFreq(1.0)],
            af_control: AlleleFreq(0.50000001)..AlleleFreq(0.75000001)
        },
        libprosic::call::pairwise::PairEvent {
            name: "germline_hom_alt".to_owned(),
            af_case: vec![AlleleFreq(1.0)],
            af_control: AlleleFreq(0.75000001)..AlleleFreq(1.00000001)
        },
        libprosic::call::pairwise::PairEvent { //TODO: optimally, this would be a two-sided event?
            name: "amplification_error_to_alt".to_owned(),
            af_case: vec![AlleleFreq(0.5)],
            af_control: AlleleFreq(0.0)..AlleleFreq(0.00000001)
        },
        libprosic::call::pairwise::PairEvent { //TODO: optimally, this would be a two-sided event?
            name: "somatic_het_from_ref".to_owned(),
            af_case: vec![AlleleFreq(0.5)],
            af_control: AlleleFreq(0.00000001)..AlleleFreq(0.25000001)
        },
        libprosic::call::pairwise::PairEvent {
            name: "germline_het".to_owned(),
            af_case: vec![AlleleFreq(0.5)],
            af_control: AlleleFreq(0.25000001)..AlleleFreq(0.75000001)
        },
        libprosic::call::pairwise::PairEvent { //TODO: optimally, this would be a two-sided event?
            name: "somatic_het_from_alt".to_owned(),
            af_case: vec![AlleleFreq(0.5)],
            af_control: AlleleFreq(0.75000001)..AlleleFreq(1.0)
        },
        libprosic::call::pairwise::PairEvent { //TODO: optimally, this would be a two-sided event?
            name: "amplification_error_to_ref".to_owned(),
            af_case: vec![AlleleFreq(0.5)],
            af_control: AlleleFreq(1.0)..AlleleFreq(1.00000001)
        }
    ];
    // call absent variants as the complement of the other events
    let absent_event = libprosic::ComplementEvent { name: "absent".to_owned() };

    let prior_model = libprosic::priors::SingleCellBulkModel::new(ploidy);

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
        &candidates,
        &output,
        &reference,
        &events,
        Some(&absent_event),
        &mut joint_model,
        omit_snvs,
        omit_indels,
        Some(max_indel_len),
        observations.as_ref(),
        exclusive_end
    )
}
