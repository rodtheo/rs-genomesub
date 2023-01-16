use assert_cmd::Command;
use predicates::prelude::*;
use std::path::Path;

type TestResult = Result<(), Box<dyn std::error::Error>>;

#[test]
fn tbltk_wot_argument() -> TestResult {
    let mut cmd = Command::cargo_bin("tbltk")?;
    cmd.arg("tofasta").arg("examples/input_multi.tbl");
    cmd.assert().failure().stderr(predicate::str::contains(
        "required arguments were not provided",
    ));

    Ok(())
}

#[test]
fn tbltk_test_stdout() -> TestResult {
    let predicate_file =
        predicate::path::eq_file(Path::new("examples/output_tbltk_test_stdout.fa"))
            .utf8()
            .unwrap();

    let mut cmd = Command::cargo_bin("tbltk")?;
    cmd.arg("tofasta")
        .arg("--genome")
        .arg("examples/input_multi_genome.fa")
        .arg("examples/input_multi.tbl");
    cmd.assert().success().stdout(predicate_file);

    Ok(())
}

#[test]
fn dies_no_args() -> TestResult {
    let mut cmd = Command::cargo_bin("tbl_rod")?;
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("USAGE"));

    Ok(())
}
