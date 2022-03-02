use assert_cmd::Command;
use predicates::prelude::*;

type TestResult = Result<(), Box<dyn std::error::Error>>;

#[test]
fn runs() -> TestResult {
    let mut cmd = Command::cargo_bin("tbl_rod")?;
    cmd.arg("-g")
        .arg("examples/input_multi.tbl")
        .assert()
        .success();

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
