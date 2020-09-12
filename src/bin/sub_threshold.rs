use histogram_minimizer::minimize;
use long_reference_speedcheck::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::path::Path;
const IS_HILL: bool = true;
const CHUNKING_POWER: f32 = 0.375;
fn main() {
    let args: Vec<_> = std::env::args().collect();
    let sam: HashMap<_, _> = BufReader::new(File::open(&args[1]).unwrap())
        .lines()
        .filter_map(|e| e.ok())
        .filter(|line| !line.starts_with('@'))
        .filter_map(|line| {
            let contents: Vec<_> = line.split('\t').take(4).collect();
            let id: String = contents[0].split('_').nth(0).unwrap().to_string();
            let flag: usize = contents[1].parse().unwrap();
            let pos: usize = contents[3].parse().unwrap();
            if flag == 0 {
                Some((id, pos))
            } else {
                None
            }
        })
        .collect();
    let queries = DataSet::from_json(&args[2]).unwrap();
    let queries: Vec<_> = queries
        .records
        .into_iter()
        .filter_map(|data| sam.get(&data.id).map(|&pos| (data.means, pos)))
        .filter(|(query, _)| query.len() > 1100)
        .take(700)
        .collect();
    let model = squiggler::Squiggler::new(&Path::new(&args[3])).unwrap();
    let (temp, rev) = setup_template_complement(&Path::new(&args[4])).unwrap();
    let temp: Vec<f32> = model
        .get_signal_from_fasta(&temp)
        .into_iter()
        .map(|e| e.2)
        .collect();
    let _temp = dtw::normalize(&temp, dtw::NormalizeType::Z);
    let rev: Vec<_> = model
        .get_signal_from_fasta(&rev)
        .into_iter()
        .map(|e| e.2)
        .collect();
    let rev = dtw::normalize(&rev, dtw::NormalizeType::Z);
    let refsize: usize = args[5].parse::<usize>().map(|e| e * 1000).expect("refsize");
    let querysize: usize = args[6].parse().expect("querysize");
    // add power when running fixed refsize
    println!("refsize,bandwidth,time,power");
    let cdata: Vec<_> = BufReader::new(File::open(&Path::new(&args[7])).expect("cdata"))
        .lines()
        .filter_map(|e| e.ok())
        .filter_map(|e| e.parse().ok())
        .collect();
    let udata: Vec<_> = BufReader::new(File::open(&Path::new(&args[8])).expect("udata"))
        .lines()
        .filter_map(|e| e.ok())
        .filter_map(|e| e.parse().ok())
        .collect();
    let (_, argmin) = minimize(udata, cdata);
    let power = CHUNKING_POWER;
    let result = speedcheck(
        &queries,
        &rev,
        querysize,
        refsize,
        power,
        dtw::Mode::QuickSub,
        IS_HILL,
        &Some(argmin as f32),
    );
    let bandwidth = 0;
    // add one culumn when runnning fixed refsize
    println!("{},{},{},{}", refsize, bandwidth, result, power);
}
