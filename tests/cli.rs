use std::path::{Path, PathBuf};
use std::process::{Command, Output};

struct Fixture(PathBuf);
impl Fixture
{
    fn new() -> Self
    {
        let path = std::env::temp_dir().join(format!("fracsync_cli_{}", std::process::id()));
        std::fs::create_dir_all(&path).unwrap();
        Self(path)
    }
    fn path(&self, name: &str) -> PathBuf
    {
        self.0.join(name)
    }
}
impl Drop for Fixture
{
    fn drop(&mut self)
    {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}
fn cli(dir: &Path, args: &[&str]) -> Output
{
    Command::new(env!("CARGO_BIN_EXE_fracsync"))
        .current_dir(dir)
        .args(args)
        .output()
        .unwrap()
}
fn success(output: Output) -> String
{
    assert!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
    String::from_utf8(output.stdout).unwrap()
}

#[test]
fn sketch_contain_info_and_format_migration()
{
    let fixture = Fixture::new();
    let dna = "ACGTTGCAACGTACGTAACCGGTTACGTTTGGCCAATGCATGCAACGTACGATCGATCGGCTA";
    std::fs::write(fixture.path("reference.fa"), format!(">ref\n{dna}\n")).unwrap();
    let rna = dna.to_lowercase().replace('t', "u");
    std::fs::write(
        fixture.path("query.fq"),
        format!("@read\n{rna}\n+\n{}\n", "I".repeat(rna.len())),
    )
    .unwrap();
    for selector in [
        vec![],
        vec!["--offset", "5"],
        vec!["--selector", "minimizer", "--w", "1"],
    ]
    {
        let mut args = vec!["sketch", "reference.fa", "-o", "db.sig"];
        args.extend(selector);
        success(cli(&fixture.0, &args));
        let info = success(cli(&fixture.0, &["info", "db.sig"]));
        let count: usize = info
            .lines()
            .nth(1)
            .unwrap()
            .split('\t')
            .nth(1)
            .unwrap()
            .parse()
            .unwrap();
        assert!(count > 0);
        let output = success(cli(&fixture.0, &["contain", "db.sig", "query.fq"]));
        let fields: Vec<_> = output.lines().nth(1).unwrap().split('\t').collect();
        assert_eq!(fields[0], "reference");
        assert_eq!(fields[1], fields[2]);
        assert_eq!(&fields[3..], &["1.0000", "1.0000", "1.0000"]);
        assert_eq!(
            output,
            success(cli(&fixture.0, &["contain", "db.sig", "query.fq", "query.fq"]))
        );
    }
    let mut bytes = std::fs::read(fixture.path("db.sig")).unwrap();
    assert_eq!(&bytes[..8], b"FRACSYN3");
    bytes[..8].copy_from_slice(b"FRACSYN1");
    std::fs::write(fixture.path("legacy.sig"), bytes).unwrap();
    let rejected = cli(&fixture.0, &["contain", "legacy.sig", "query.fq"]);
    assert!(!rejected.status.success());
    assert!(String::from_utf8_lossy(&rejected.stderr).contains("re-sketch"));
    let bad = cli(&fixture.0, &["sketch", "reference.fa", "-o", "bad.sig", "--s", "0"]);
    assert!(!bad.status.success());
    assert!(!fixture.path("bad.sig").exists());
}
