use rayon::prelude::*;
use serde::*;
use std;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use std::time::{Duration, Instant};
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DataSet {
    pub records: Vec<Data>,
}

impl DataSet {
    pub fn from_json<P: AsRef<Path>>(path: P) -> serde_json::Result<DataSet> {
        let rdr = BufReader::new(File::open(path).unwrap());
        serde_json::de::from_reader(rdr)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Data {
    pub id: String,
    pub means: Vec<f32>,
}

pub fn setup_template_complement(path: &Path) -> std::result::Result<(String, String), ()> {
    let seq = bio::io::fasta::Reader::from_file(path).map_err(|_| ())?;
    let reference: Vec<u8> = seq
        .records()
        .filter_map(|e| e.ok())
        .fold(Vec::new(), |mut acc, e| {
            acc.extend_from_slice(&mut e.seq());
            acc
        });
    let reverse = bio::alphabets::dna::revcomp(&reference);
    let temp = String::from_utf8(reference).map_err(|_| ())?;
    let rev = String::from_utf8(reverse).map_err(|_| ())?;
    Ok((temp, rev))
}

pub fn speedcheck(
    queries: &Vec<(Vec<f32>, usize)>,
    temp: &Vec<f32>,
    querysize: usize,
    refsize: usize,
    power: f32,
    mode: dtw::Mode,
    is_hill: bool,
    threshold: &Option<f32>,
) -> f64 {
    let res: Vec<_> = queries
        .par_iter()
        .filter_map(|&(ref query, location)| {
            elapsed_time(
                &query, location, temp, querysize, refsize, power, mode, is_hill, threshold,
            )
        })
        .collect();
    let (duration, num): (_, u32) = res
        .into_iter()
        .fold((Duration::from_millis(0), 0), |(acc, num), x| {
            (x + acc, num + 1)
        });
    duration.as_secs() as f64 / num as f64
        + duration.subsec_nanos() as f64 / num as f64 / 10f64.powi(9)
}

pub fn elapsed_time(
    query: &Vec<f32>,
    location: usize,
    reference: &Vec<f32>,
    querysize: usize,
    refsize: usize,
    power: f32,
    mode: dtw::Mode,
    is_hill: bool,
    threshold: &Option<f32>,
) -> Option<Duration> {
    if query.len() < 50 + 2 * querysize
        || location + refsize > 200 + reference.len()
        || location < 200
    {
        return None;
    } else {
        let metric = if is_hill { "hill" } else { "normal" };
        let query = query[50..50 + 2 * querysize].to_vec();
        let query = dtw::normalize(&query, dtw::NormalizeType::Z);
        let query = squiggler::dedup(&query, power);
        let query: Vec<_> = query.into_iter().take(querysize).collect(); // algorithm modification
        let now = Instant::now();
        let res = dtw::utils::dtw_wrapper(
            &query,
            &reference[location - 200..location + refsize - 200],
            &mode,
            metric,
            &None,
            threshold,
        );
        let time = now.elapsed();
        eprint!("{:?}", res);
        Some(time)
    }
}
