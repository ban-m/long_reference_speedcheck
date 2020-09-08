extern crate rayon;
extern crate bio;
extern crate fast5wrapper;
extern crate dtw;
extern crate squiggler;
extern crate rand;
extern crate csv;
use rayon::prelude::*;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::path::Path;
use csv::WriterBuilder;
const CHUNKING_POWER:f32 = 0.4;
const BINS:usize = 20;
fn prepare_query(path:&str)->fast5wrapper::result::Result<(String,Vec<f32>)>{
    let query:Vec<_> = fast5wrapper::get_event(path)?
    .into_iter().map(|e|e[2]).collect();
    let id:String = fast5wrapper::get_read_id(path)?;
    Ok((id,query))
}

fn main(){
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
                        Ok((id,events)) => Some((id,events,location)),
                        Err(_) => None,
                    })
        .take(10)
        .collect();
    let model = squiggler::Squiggler::new(&Path::new(&args[2])).unwrap();
    let (temp,rev) = setup_template_complement(&Path::new(&args[3])).unwrap();
    let temp:Vec<f32> = model.get_signal_from_fasta(&temp).into_iter().map(|e|e.2).collect();
    let temp = dtw::normalize(&temp,dtw::NormalizeType::Z);
    let rev :Vec<_> = model.get_signal_from_fasta(&rev)
        .into_iter().map(|e|e.2).collect();
    let _rev = dtw::normalize(&rev,dtw::NormalizeType::Z);
    let refsize:usize= args[4].parse::<usize>().expect("refsize") * 1_000;
    let querysize:usize =  args[5].parse().expect("querysize");
    let max_bandwidth:usize = args[6].parse().expect("bandwidth");
    let mut wtr = match WriterBuilder::new().from_path("./result/scorepath.csv"){
        Ok(res)=>res,
        Err(why)=>panic!("{}",why)
    };
    if let Err(why) = wtr.write_record(&vec!["position",
                                             "bandwidth",
                                             "id",
                                             "score",
                                             "kl_divergence"]){
        eprintln!("{}",why);
    }
    let mut parameters:Vec<_> = (1..max_bandwidth/2).map(|e|2*e+1)
        .map(|e|dtw::Mode::SakoeChiba(e)).collect();
    parameters.push(dtw::Mode::Sub);
    for (id,events,location) in queries{
        let result:Vec<_> = parameters.par_iter()
            .map(|mode|
                 scorepath(&events,location,&temp,querysize,refsize,mode))
            .reduce(||vec![],|mut acc,mut x|{acc.append(&mut x);acc});
        for (position,bandwidth,score,kl) in result{
            if let Err(w) =  wtr.write_record(&vec![format!("{}",position),
                                                    format!("{}",bandwidth),
                                                    format!("{}",id),
                                                    format!("{}",score),
                                                    format!("{}",kl)]){
                eprintln!("{:?}",w);
            }
        }
    }
}


fn setup_template_complement(path:&Path)->std::result::Result<(String,String),()>{
    let seq = bio::io::fasta::Reader::from_file(path).map_err(|_|())?;
    let reference:Vec<u8> = seq.records()
        .filter_map(|e|e.ok())
        .fold(Vec::new(),|mut acc,e|{acc.extend_from_slice(&mut e.seq());acc});
    let reverse = bio::alphabets::dna::revcomp(&reference);
    let temp = String::from_utf8(reference).map_err(|_|())?;
    let rev = String::from_utf8(reverse).map_err(|_|())?;
    Ok((temp,rev))
}



#[test]
fn kltest(){
    let p = vec![0.2;20];
    let q = vec![0.2;20];
    assert!(kl(&p,&q)<0.1)
}

#[test]
fn kltest2(){
    let p = noise(20);
    let q = noise(40);
    assert!(kl(&p,&q) > 0.,"{}",kl(&p,&q))
}

#[test]
fn kltest3(){
    let p = noise(20);
    let q = noise(40);
    assert!(!kl(&p,&q).is_nan(),"{}",kl(&p,&q))
}

#[test]
fn kltest5(){
    let p = noise(100);
    let q = noise(100);
    assert!(kl(&p,&q)>0.0,"{}",kl(&p,&q))
}

#[test]
fn kltest4(){
    let p = noise(20);
    let q = p.clone();
    assert!(kl(&p,&q) < 0.001,"{}",kl(&p,&q))
}    

#[inline]
fn histogram(p:&[f32])->Vec<f64>{
    let mut histogram = vec![0.;BINS];
    let len = p.len() as f64;
    for &x in  p{
        histogram[bin_index_of(x)] += 1.0;
    }
    histogram = histogram.into_iter()
        .map(|e|e/len).collect();
    debug_assert!((histogram.iter().sum::<f64>() - 1.0f64).abs()<0.01,"{}",histogram.iter().sum::<f64>());
    histogram
}

#[inline]
fn histogram_unif(p:&[f32])->Vec<f64>{
    // prior
    let mut histogram = vec![1. as f64;BINS];
    let len = (p.len()+BINS) as f64;
    for &x in p{
        histogram[bin_index_of(x)] += 1.0;
    }
    histogram = histogram.into_iter()
        .map(|e|e/len).collect();
    debug_assert!((histogram.iter().sum::<f64>() - 1.0f64).abs()<0.01,"{}",histogram.iter().sum::<f64>());
    histogram
}

#[inline]
fn kl(ps:&[f64],qs:&[f64]) -> f64{
    ps.iter().zip(qs.iter())
        .map(|(&p,&q)| if p == 0.{
            0.
        }else if q == 0. {
            1000000.
        }else{            
            p*(p/q).ln()
        })
        .sum()
}
#[inline]
fn bin_index_of(x:f32)->usize{
    let index = ((x+2.)/4. * BINS as f32).floor();
    if index < 0. {
        0
    }else if index >= BINS as f32 {
        BINS-1
    }else{
        index as usize
    }
}

fn kl_divergences(query:&[f32],reference:&[f32],bandwidth:usize) -> Vec<Vec<f64>>{
    // compute kl[i][j] = kl divergence between
    // query vs reference[i..querylength-3*bandwidth/2+j](j = 0..bandwidth)
    // i.e, sum(p_k log(pk\qk)) where p_k is the hidtogram of the query,
    // q_i is a histogram of the reference.
    // Note that we assume that query is true model, and
    // want to find the best matching interval in reference which
    // is simillar to query.
    // Currently, kl_divergence bandwidth is restricted in 2*bandsize
    // region only. To more information, loook at readu_otherpictures
    // and get intuition abount "corresponding region".
    // Also, binwidth is one of the most important part of kl divergence
    // computation. As far, I use [-2,2] region for compute histogram,
    // and 20 bins(width = 4/20 = 0.2).
    let query_histogram = histogram(&query);
    let reference_first_hitrogram=histogram_unif(&reference[0..query.len() -bandwidth*3/4]);
    //prior distribution
    // counting[i][k] = how many k in ref[0..i].
    let mut counting_table = vec![vec![0;BINS];reference.len()];
    for (index,&r) in reference.iter().enumerate(){
        let mark = bin_index_of(r);
        if index != 0{
            for k in 0..BINS{
                if k == mark{
                    counting_table[index][mark] = counting_table[index-1][mark]+1;
                }else{
                    counting_table[index][k] = counting_table[index-1][k];
                }
            }
        }else{
            for k in 0..BINS{
                if k == mark{
                    counting_table[0][mark] = 1; 
                }else{
                    counting_table[0][k] = 0;
                }
            }
        }
    }
    let mut res = vec![vec![0.;bandwidth];reference.len()];
    res[0][0] = kl(&query_histogram,&reference_first_hitrogram);
    debug_assert!(res[0][0]>0.);
//    eprintln!("{}",res[0][0]);
    for i in 0..reference.len(){
        for j in 0..bandwidth{
            if i == 0 && j == 0 {
                continue;
            }else if j + i + query.len() - 3*bandwidth/4 >= reference.len(){
                res[i][j] = 100000000.;//something big
            }else{
                // compute kl divergence between
                // query histogram vs reference[i..truej]
                let truej = j + i + query.len() - 3*bandwidth/4;
                let refhist = histogram_unif(&reference[i..truej]);
                res[i][j] = kl(&query_histogram,&refhist);
                // let index = bin_index_of(reference[truej]);
                // if j == 0{
                //     let update = ((truej-i+2+BINS) as f64/(truej-i+1+BINS) as f64).ln();
                //     res[i][j] = res[i-1][1] - update;
                // }else{
                //     let n = if i == 0{
                //         counting_table[truej-1][index]
                //     }else{
                //         counting_table[truej-1][index] - counting_table[i-1][index]
                //     };
                //     let update = ((truej-i+BINS) as f64/(truej-i+1+BINS) as f64).ln() +
                //         query_histogram[index]*(1.+1./(n+1) as f64).ln();
                //     res[i][j] = res[i][j-1] - update;
                // }
//                eprintln!("{},{}->{},",i,truej,res[i][j]);
            }
        }
    }
    for r in 0..res.len(){
        for q in 0..res[r].len(){
            if res[r][q] == 0.{
                eprintln!("{},{}->0.",r,q);
            }
        }
    }
    res
}

fn kl_divergence(query:&[f32],reference:&[f32],mode:&dtw::Mode)->Vec<f64>{
    match mode{
        &dtw::Mode::Sub => vec![0.],
        &dtw::Mode::SakoeChiba(b) => {
            let kl_divergences = kl_divergences(query,reference,b);
            (0..reference.len()-query.len()+1)
                .map(|i|kl_divergences[i].iter()
                     .fold(10000000.,|acc,&x|if acc < x { acc } else{ x }))
                .collect()
        }
        _ => vec![]
    }
}

fn scorepath(query:&Vec<f32>,location:usize,
             reference:&Vec<f32>,
             querysize:usize,refsize:usize,
             mode:&dtw::Mode)->Vec<(usize,usize,f32,f64)>{
    if query.len() < 50 + querysize || 
        location + refsize > 2000 + reference.len() || 
        location < 2000{
            return vec![]
        }else{
            let query = query[50..50+querysize].to_vec();
            let query = dtw::normalize(&query,dtw::NormalizeType::Z);
            let query = squiggler::dedup(&query,CHUNKING_POWER);
            let subref = &reference[location-2000..location+refsize-2000];
            let kl = kl_divergence(&query,
                                   subref,
                                   mode);
            match mode{
                &dtw::Mode::Sub => {
                    if let Ok(res) = dtw::dtw(&query,
                                              subref,
                                              mode.clone(),&hill){
                        vec![(res.2,0,res.0,0.)]
                    }else{
                        vec![]
                    }},
                &dtw::Mode::SakoeChiba(b) =>
                    subref
                    .windows(query.len())
                    .enumerate()
                    .filter_map(|(index,window)|
                                dtw::dtw(&query,
                                         &window,
                                         dtw::Mode::SakoeChiba(b),
                                         &hill)
                                .map(|e|(index,b,e.0))
                                .ok())
                    .zip(kl.iter())
                    .map(|((index,b,score),&kl)|(index,b,score,kl))
                    .collect(),
                _ => unreachable!(),
            }
        }
}


fn hill(x:&f32,y:&f32)->f32{
    let d = (x-y).powi(2);
    d/(1.0 + d)
}

