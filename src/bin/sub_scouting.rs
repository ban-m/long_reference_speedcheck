extern crate rayon;
extern crate bio;
extern crate fast5wrapper;
extern crate dtw;
extern crate squiggler;
extern crate rand;
extern crate histogram_minimizer;
use std::fs::File;
use std::io::{BufRead,BufReader};
use std::path::Path;
pub mod utils;
use utils::*;
const IS_HILL:bool = true;
fn main() {
    let args:Vec<_> = std::env::args().collect();
    let queries:Vec<_> = BufReader::new(File::open(&Path::new(&args[1])).unwrap())
        .lines().filter_map(|e|e.ok())
        .filter_map(|e|{
            let content:Vec<_> =  e.split(',').collect();
            let filename = content[0].to_string();
            let location = content[1].parse::<usize>().unwrap();
            Some((filename,location))
        })
        .filter_map(|(filename,location)|
                    match prepare_query(&filename){
                        Ok(res) => Some((res,location)),
                        Err(_) => None,
                    })
        .take(700)
        .collect();
    let model = squiggler::Squiggler::new(&Path::new(&args[2])).unwrap();
    let (temp,rev) = setup_template_complement(&Path::new(&args[3])).unwrap();
    let temp:Vec<f32> = model.get_signal_from_fasta(&temp).into_iter().map(|e|e.2).collect();
    let _temp = dtw::normalize(&temp,dtw::NormalizeType::Z);
    let rev :Vec<_> = model.get_signal_from_fasta(&rev)
        .into_iter().map(|e|e.2).collect();
    let rev = dtw::normalize(&rev,dtw::NormalizeType::Z);
    let refsize:usize = args[4].parse::<usize>().map(|e|e*1000).expect("refsize");
    let querysize:usize =  args[5].parse().expect("querysize");
    let num_scouts:usize = args[6].parse().expect("# of scouts");
    let num_packs:usize = args[7].parse().expect("# of packs");
    let power:usize = args[8].parse().expect("power");
    // add power when running fixed refsize
    let result = speedcheck(&queries,
                            &rev,
                            querysize,
                            refsize,
                            power as f32 / 100.,
                            dtw::Mode::Scouting(num_scouts,num_packs),
                            IS_HILL,
                            &None);
    println!("{},{},{},{},{}",
             refsize,num_scouts,num_packs,
             result,power);
}



