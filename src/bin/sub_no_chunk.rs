extern crate rayon;
extern crate bio;
extern crate fast5wrapper;
extern crate dtw;
extern crate squiggler;
extern crate rand;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::path::Path;
mod utils;
const IS_HILL:bool = false;
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
    let temp = dtw::normalize(&temp,dtw::NormalizeType::Z);
    let rev :Vec<_> = model.get_signal_from_fasta(&rev)
        .into_iter().map(|e|e.2).collect();
    let _rev = dtw::normalize(&rev,dtw::NormalizeType::Z);
    let refsizes :Vec<_>= args[4].split(',').filter_map(|e|e.parse::<usize>().ok())
        .map(|e|e*1000).collect();
    let querysize:usize =  args[5].parse().expect("querysize");
    // add power when running fixed refsize
    println!("refsize,time,power");
    let mode = dtw::Mode::Sub;
    for refsize in refsizes{
        let power = 0.0;
        let result = speedcheck(&queries,&temp,querysize,refsize,power,mode.clone(),
                                IS_HILL,&None);
        // add one culumn when runnning fixed refsize
        println!("{},{},{}",
                 refsize,
                 result,
                 power);
    }
}

