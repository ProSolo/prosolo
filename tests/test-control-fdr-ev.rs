extern crate cli_test_dir;
extern crate regex;

use cli_test_dir::*;
use regex::Regex;
use std::fs::File;
use std::io::prelude::*;

#[test]
fn test_control_fdr_ev() {
    let control_fdr_ev = TestDir::new("prosolo", "control-fdr-ev");
    let output = control_fdr_ev.cmd()
        .arg("control-fdr")
        .arg(control_fdr_ev.src_path("tests/expected-out_omit-indels.bcf"))
        .arg("--events").arg("HOM_REF")
        .arg("--var").arg("SNV")
        .arg("--method").arg("ev")
        .output()
        .expect_success();

    //test that the output thresholds list is the same as before
    let expected = control_fdr_ev.src_path("tests/expected-out_fdr-ev.tsv");
    let mut str_expected = String::new();
    if let Ok(_) = File::open(expected).unwrap().read_to_string(&mut str_expected) {
        assert_eq!(output.stdout_str(), str_expected);
    } else {
        panic!("Couldn't read expected tsv file.");
    };
}
