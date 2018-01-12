extern crate cli_test_dir;
// TODO: consider changing to assert_cli crate, looks like a cleaner syntax and might allow removing the leading underscores for `output` in most examples to suppress unused variable warnings

use cli_test_dir::*;
use std::fs::File;
use std::io::prelude::*;

#[test]
fn test_control_fdr_ev() {
    let control_fdr_ev = TestDir::new("prosolo", "control-fdr-ev");
    let output = control_fdr_ev.cmd()
        .arg("control-fdr")
        .arg(control_fdr_ev.src_path("tests/expected-out_omit-indels.bcf"))
        .arg("--events").arg("ADO_TO_ALT,ADO_TO_REF")
        .arg("--var").arg("SNV")
        .arg("--method").arg("ev")
        .output()
        .expect_success();

    //test that the output thresholds list is the same as before
    let expected = control_fdr_ev.src_path("tests/expected-out_ADO-fdr-ev.tsv");
    let mut str_expected = String::new();
    if let Ok(_) = File::open(expected).unwrap().read_to_string(&mut str_expected) {
        assert_eq!(output.stdout_str(), str_expected);
    } else {
        panic!("Couldn't read expected TSV file.");
    };
}

#[test]
fn test_apply_fdr_ev() {
    let apply_fdr_ev = TestDir::new("prosolo", "apply-fdr");
    let _output = apply_fdr_ev.cmd()
        .arg("apply-fdr")
        .arg(apply_fdr_ev.src_path("tests/expected-out_omit-indels.bcf"))
        .arg("--threshold").arg("0.1819080263")
        .arg("--events").arg("ADO_TO_ALT,ADO_TO_REF")
        .arg("--var").arg("SNV")
        .arg("--output").arg("test-out_ADO_0-04.bcf")
        .output()
        .expect_success();

    //test that the output BCF file is the same as before
    let expected = apply_fdr_ev.src_path("tests/expected-out_ADO_0-04.bcf");
    let mut bcf_expected: Vec<u8> = Vec::new();
    if let Ok(_) = File::open(expected).unwrap().read_to_end(&mut bcf_expected) {
        apply_fdr_ev.expect_file_contents("test-out_ADO_0-04.bcf", bcf_expected );
    } else {
        panic!("Couldn't read expected BCF file.");
    };
}

#[test]
fn test_apply_fdr_ev_filter_1() {
    let apply_fdr_ev = TestDir::new("prosolo", "apply-fdr");
    let _output = apply_fdr_ev.cmd()
        .arg("apply-fdr")
        .arg(apply_fdr_ev.src_path("tests/control-fdr-tests/alt-prob-sum-above-one_test-sites.bcf"))
        .arg("--threshold").arg("0.000001")
        .arg("--events").arg("ADO_TO_REF,ADO_TO_ALT,HOM_ALT,HET,ERR_REF")
        .arg("--var").arg("SNV")
        .arg("--output").arg("test-out_ALT_filter-1.bcf")
        .output()
        .expect_success();

    //test that the output BCF file is the same as before
    let expected = apply_fdr_ev.src_path("tests/control-fdr-tests/expected-out_ALT_filter-1.bcf");
    let mut bcf_expected: Vec<u8> = Vec::new();
    if let Ok(_) = File::open(expected).unwrap().read_to_end(&mut bcf_expected) {
        apply_fdr_ev.expect_file_contents("test-out_ALT_filter-1.bcf", bcf_expected );
    } else {
        panic!("Couldn't read expected BCF file.");
    };
}

#[test]
fn test_apply_fdr_ev_filter_0() {
    let apply_fdr_ev = TestDir::new("prosolo", "apply-fdr");
    let _output = apply_fdr_ev.cmd()
        .arg("apply-fdr")
        .arg(apply_fdr_ev.src_path("tests/control-fdr-tests/alt-prob-sum-above-one_test-sites.bcf"))
        .arg("--threshold").arg("0.00001")
        .arg("--events").arg("ADO_TO_REF,ADO_TO_ALT,HOM_ALT,HET,ERR_REF")
        .arg("--var").arg("SNV")
        .arg("--output").arg("test-out_ALT_filter-0.bcf")
        .output()
        .expect_success();

    //test that the output BCF file is the same as before
    let expected = apply_fdr_ev.src_path("tests/control-fdr-tests/alt-prob-sum-above-one_test-sites.bcf");
    let mut bcf_expected: Vec<u8> = Vec::new();
    if let Ok(_) = File::open(expected).unwrap().read_to_end(&mut bcf_expected) {
        apply_fdr_ev.expect_file_contents("test-out_ALT_filter-0.bcf", bcf_expected );
    } else {
        panic!("Couldn't read expected BCF file.");
    };
}
