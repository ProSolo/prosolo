extern crate cli_test_dir;
// TODO: consider changing to assert_cli crate, looks like a cleaner syntax and might allow removing the leading underscores for `output` in most examples to suppress unused variable warnings

use cli_test_dir::*;
use std::fs::File;
use std::io::prelude::*;

#[test]
fn test_control_fdr_ado_0_04() {
    let control_fdr_ev = TestDir::new("prosolo", "control-fdr-ado-0-04");
    let _output = control_fdr_ev.cmd()
        .arg("control-fdr")
        .arg(control_fdr_ev.src_path("tests/expected-out_omit-indels.bcf"))
        .arg("--fdr").arg("0.04")
        .arg("--events").arg("ADO_TO_ALT,ADO_TO_REF")
        .arg("--var").arg("SNV")
        .arg("--output").arg("test-out_ADO_0-04.bcf")
        .output()
        .expect_success();

    //test that the output BCF file is the same as before
    let expected = control_fdr_ev.src_path("tests/expected-out_ADO_0-04.bcf");
    let mut bcf_expected: Vec<u8> = Vec::new();
    if let Ok(_) = File::open(expected).unwrap().read_to_end(&mut bcf_expected) {
        control_fdr_ev.expect_file_contents("test-out_ADO_0-04.bcf", bcf_expected );
    } else {
        panic!("Couldn't read expected BCF file.");
    };
}

#[test]
fn test_control_fdr_alt_filter_one() {
    let control_fdr_ev = TestDir::new("prosolo", "control-fdr-alt-filter-one");
    let _output = control_fdr_ev.cmd()
        .arg("control-fdr")
        .arg(control_fdr_ev.src_path("tests/control-fdr-tests/alt-prob-sum-above-one_test-sites.bcf"))
        .arg("--fdr").arg("0.000000001")
        .arg("--events").arg("ADO_TO_REF,ADO_TO_ALT,HOM_ALT,HET,ERR_REF")
        .arg("--var").arg("SNV")
        .arg("--output").arg("test-out_ALT_filter-one.bcf")
        .output()
        .expect_success();

    //test that the output BCF file is the same as before
    let expected = control_fdr_ev.src_path("tests/control-fdr-tests/expected-out_ALT_filter-one.bcf");
    let mut bcf_expected: Vec<u8> = Vec::new();
    if let Ok(_) = File::open(expected).unwrap().read_to_end(&mut bcf_expected) {
        control_fdr_ev.expect_file_contents("test-out_ALT_filter-one.bcf", bcf_expected );
    } else {
        panic!("Couldn't read expected BCF file.");
    };
}

#[test]
fn test_control_fdr_alt_filter_none() {
    let control_fdr_ev = TestDir::new("prosolo", "control-fdr-alt-filter_none");
    let _output = control_fdr_ev.cmd()
        .arg("control-fdr")
        .arg(control_fdr_ev.src_path("tests/control-fdr-tests/alt-prob-sum-above-one_test-sites.bcf"))
        .arg("--fdr").arg("1")
        .arg("--events").arg("ADO_TO_REF,ADO_TO_ALT,HOM_ALT,HET,ERR_REF")
        .arg("--var").arg("SNV")
        .arg("--output").arg("test-out_ALT_filter-none.bcf")
        .output()
        .expect_success();

    //test that the output BCF file is the same as before
    let expected = control_fdr_ev.src_path("tests/control-fdr-tests/expected-out_ALT_filter-none.bcf");
    let mut bcf_expected: Vec<u8> = Vec::new();
    if let Ok(_) = File::open(expected).unwrap().read_to_end(&mut bcf_expected) {
        control_fdr_ev.expect_file_contents("test-out_ALT_filter-none.bcf", bcf_expected );
    } else {
        panic!("Couldn't read expected BCF file.");
    };
}
