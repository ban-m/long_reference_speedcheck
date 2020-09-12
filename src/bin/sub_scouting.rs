const IS_HILL: bool = true;
use long_reference_speedcheck::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
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
    let queries: DataSet = DataSet::from_json(&args[2]).unwrap();
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
    let num_scouts: usize = args[7].parse().expect("# of scouts");
    let num_packs: usize = args[8].parse().expect("# of packs");
    let power: usize = args[9].parse().expect("power");
    // add power when running fixed refsize
    let result = speedcheck(
        &queries,
        &rev,
        querysize,
        refsize,
        power as f32 / 100.,
        dtw::Mode::Scouting(num_scouts, num_packs),
        IS_HILL,
        &None,
    );
    println!(
        "{},{},{},{},{}",
        refsize, num_scouts, num_packs, result, power
    );
}
