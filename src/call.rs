use std::error::Error;

use clap;
use libprosic;
use libprosic::model::AlleleFreq;
use rust_htslib::bam;
use bio::stats::Prob;

pub fn single_cell_bulk(matches: &clap::ArgMatches) -> Result<(), Box<Error>> {
    // read command line parameters
    let single_mean_insert_size = value_t!(matches, "single-insert-size-mean", f64).unwrap();
    let single_sd_insert_size = value_t!(matches, "single-insert-size-sd", f64).unwrap();
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
    let omit_snvs = matches.is_present("omit-snvs");
    let omit_indels = matches.is_present("omit-indels");
    let single = matches.value_of("single").unwrap();
    let bulk = matches.value_of("bulk").unwrap();
    let candidates = matches.value_of("candidates").unwrap_or("-");
    let output = matches.value_of("output").unwrap_or("-");
    let observations = matches.value_of("observations");
    let flat_priors = matches.is_present("flat-priors");

    let prob_spurious_isize = try!(Prob::checked(value_t!(matches, "prob-spurious-isize", f64).unwrap_or(0.0)));
    let prob_missed_insertion_alignment = try!(Prob::checked(value_t!(matches, "prob-missed-insertion-alignment", f64).unwrap_or(0.0)));
    let prob_missed_deletion_alignment = try!(Prob::checked(value_t!(matches, "prob-missed-deletion-alignment", f64).unwrap_or(0.0)));
    let prob_spurious_indel_alignment = try!(Prob::checked(value_t!(matches, "prob-spurious-indel-alignment", f64).unwrap_or(0.0)));

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
        libprosic::InsertSize {
            mean: bulk_mean_insert_size,
            sd: bulk_sd_insert_size
        },
        libprosic::likelihood::LatentVariableModel::with_single_sample(),
        prob_spurious_isize,
        prob_missed_insertion_alignment,
        prob_missed_deletion_alignment,
        prob_spurious_indel_alignment
    ).max_indel_dist(max_indel_dist).max_indel_len_diff(max_indel_len_diff);

    // init single sample
    let single_sample = libprosic::Sample::new(
        single_bam,
        pileup_window,
        !no_fragment_evidence,
        !no_secondary,
        !no_mapq,
        libprosic::InsertSize {
            mean: single_mean_insert_size,
            sd: single_sd_insert_size
        },
        libprosic::likelihood::LatentVariableModel::with_single_sample(),
        prob_spurious_isize,
        prob_missed_insertion_alignment,
        prob_missed_deletion_alignment,
        prob_spurious_indel_alignment
    ).max_indel_dist(max_indel_dist).max_indel_len_diff(max_indel_len_diff);

    // setup events: case = single cell; control = bulk
    let events = [
        libprosic::call::pairwise::PairEvent {
            name: "germline_hom_ref".to_owned(),
            af_case: vec![AlleleFreq(0.0)],
            af_control: AlleleFreq(0.0)..AlleleFreq(0.5)
        },
        libprosic::call::pairwise::PairEvent {
            name: "somatic_hom_ref".to_owned(),
            af_case: vec![AlleleFreq(0.0)],
            af_control: AlleleFreq(0.5)..AlleleFreq(1.0)
        },
        libprosic::call::pairwise::PairEvent {
            name: "somatic_hom_alt".to_owned(),
            af_case: vec![AlleleFreq(1.0)],
            af_control: AlleleFreq(0.0)..AlleleFreq(0.5)
        },
        libprosic::call::pairwise::PairEvent {
            name: "germline_hom_alt".to_owned(),
            af_case: vec![AlleleFreq(1.0)],
            af_control: AlleleFreq(0.5)..AlleleFreq(1.0)
        },
        libprosic::call::pairwise::PairEvent { //TODO: optimally, this would be a two-sided event?
            name: "somatic_het_from_ref".to_owned(),
            af_case: vec![AlleleFreq(0.5)],
            af_control: AlleleFreq(0.0)..AlleleFreq(0.25)
        },
        libprosic::call::pairwise::PairEvent {
            name: "germline_het".to_owned(),
            af_case: vec![AlleleFreq(0.5)],
            af_control: AlleleFreq(0.25)..AlleleFreq(0.75)
        },
        libprosic::call::pairwise::PairEvent { //TODO: optimally, this would be a two-sided event?
            name: "somatic_het_from_alt".to_owned(),
            af_case: vec![AlleleFreq(0.5)],
            af_control: AlleleFreq(0.75)..AlleleFreq(1.0)
        }
    ];
    // call absent variants as the complement of the other events
    let absent_event = libprosic::ComplementEvent { name: "absent".to_owned() };

    let prior_model = libprosic::priors::SingleCellBulkModel::new(ploidy);

    // init joint model
    let mut joint_model = libprosic::model::PairModel::new(
        single_sample,
        bulk_sample,
        prior_model
    );

    // perform calling
    libprosic::call::pairwise::call::<
        _, _, _,
        libprosic::model::PairModel<
            libprosic::model::ContinuousAlleleFreqs,
            libprosic::model::DiscreteAlleleFreqs,
            libprosic::model::priors::SingleCellBulkModel
        >, _, _, _>
    (
        &candidates,
        &output,
        &events,
        Some(&absent_event),
        &mut joint_model,
        omit_snvs,
        omit_indels,
        Some(max_indel_len),
        observations.as_ref()
    )
}
