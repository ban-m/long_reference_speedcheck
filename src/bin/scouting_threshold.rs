extern crate rayon;
extern crate bio;
extern crate fast5wrapper;
extern crate dtw;
extern crate squiggler;
extern crate rand;
extern crate csv;
extern crate histogram_minimizer;
use histogram_minimizer::minimize;
use std::io::BufRead;
use std::io::BufReader;
use std::path::Path;
use std::fs::File;
const IS_HILL:bool = true;
mod utils;
use utils::*;
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
    let num_scouts:usize = args[6].parse().expect("scouts");
    let num_packs:usize = args[7].parse().expect("packs");
    let power:usize = args[8].parse().expect("power");
    let cdata:Vec<_> = BufReader::new(File::open(&Path::new(&args[9])).expect("cdata"))
        .lines()
        .filter_map(|e|e.ok())
        .filter_map(|e|e.parse().ok())
        .collect();
    let udata:Vec<_> = BufReader::new(File::open(&Path::new(&args[10])).expect("udata"))
        .lines()
        .filter_map(|e|e.ok())
        .filter_map(|e|e.parse().ok())
        .collect();
    let (_,argmin) = minimize(udata,cdata);
    let result = speedcheck(&queries,&rev,querysize,refsize,power as f32 / 100.,
                            dtw::Mode::Scouting(num_scouts,num_packs),
                            IS_HILL,
                            &Some(argmin as f32));
    // add one culumn when runnning fixed refsize
    println!("{},{},{},{},{}",
             refsize,
             num_scouts,
             num_packs,
             result,
             power)
}


