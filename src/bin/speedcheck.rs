const IS_HILL: bool = true;
use long_reference_speedcheck::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
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
    let temp = dtw::normalize(&temp, dtw::NormalizeType::Z);
    let rev: Vec<_> = model
        .get_signal_from_fasta(&rev)
        .into_iter()
        .map(|e| e.2)
        .collect();
    let _rev = dtw::normalize(&rev, dtw::NormalizeType::Z);
    let refsizes: Vec<_> = args[5]
        .split(',')
        .filter_map(|e| e.parse::<usize>().ok())
        .map(|e| e * 1000)
        .collect();
    let querysize: usize = args[6].parse().expect("querysize");
    let max_bandwidth: usize = args[7].parse().expect("bandwidth");
    let power: usize = args[8].parse().expect("power");
    // add power when running fixed refsize
    let mut parameters: Vec<_> = (1..max_bandwidth / 10)
        .map(|e| 10 * e + 1)
        .map(|e| dtw::Mode::SakoeChiba(e))
        .collect();
    if parameters.is_empty() {
        parameters.push(dtw::Mode::Sub);
    }
    for refsize in refsizes {
        for mode in parameters.iter() {
            let result = speedcheck(
                &queries,
                &temp,
                querysize,
                refsize,
                power as f32 / 100.,
                mode.clone(),
                IS_HILL,
                &None,
            );
            let bandwidth = match mode {
                &dtw::Mode::SakoeChiba(b) => b,
                &dtw::Mode::Sub => 0,
                _ => 0,
            };
            // add one culumn when runnning fixed refsize
            println!("{},{},{},{}", refsize, bandwidth, result, power);
        }
    }
}
