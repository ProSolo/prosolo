extern crate cli_test_dir;
extern crate regex;

use cli_test_dir::*;
use std::fs::File;
use std::io::prelude::*;

#[test]
fn test_simulated_omit_indels() {
    let omit_indels_run = TestDir::new("prosolo", "simulated-omit-indels");
    let _output = omit_indels_run
        .cmd()
        .arg("-v")
        .arg("single-cell-bulk")
        .arg("--omit-indels")
        .arg("--candidates")
        .arg(omit_indels_run.src_path("tests/candidates.bcf"))
        //        .arg("--obs").arg("observations_omit-indels.out")
        .arg("--output")
        .arg("test-out_omit-indels.bcf")
        .arg(omit_indels_run.src_path("tests/single-cell.bam"))
        .arg(omit_indels_run.src_path("tests/bulk.bam"))
        .arg(omit_indels_run.src_path("tests/ref.fa"))
        .output()
        .expect_success();
    //    omit_indels_run.expect_path("observations_omit-indels.out");
    omit_indels_run.expect_path("test-out_omit-indels.bcf");

    //test that the output BCF file is the same as before
    let bcf_path = omit_indels_run.src_path("tests/expected-out_omit-indels.bcf");
    let mut bcf_expected: Vec<u8> = Vec::new();
    if let Ok(_) = File::open(bcf_path).unwrap().read_to_end(&mut bcf_expected) {
        omit_indels_run.expect_file_contents("test-out_omit-indels.bcf", bcf_expected);
    } else {
        panic!("Couldn't read expected bcf file.");
    };
}
