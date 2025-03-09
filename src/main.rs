use std::io::{self, Write};

use rand::rng;
use rand::Rng;
use rand_distr::{Distribution, Normal};

mod plots;
#[allow(unused_imports)]
use plots::resolution::{figure_1, figure_2};



#[allow(dead_code)]
fn input(prompt: &str) -> String {
    print!("{}", prompt);
    io::stdout().flush().expect("Failed to flush stdout.");
    let mut buffer = String::new();
    io::stdin().read_line(&mut buffer).expect("Failed to read line.");
    buffer.trim().to_string()
}


/// Generate a uniform partition on the time interval
/// [a, b] with 'k' values, 'k-1' increments.
fn linspace(a: f64, b: f64, k: usize) -> Vec<f64> {
    let step_size = (b - a) / (k as f64 - 1.0);
    let mut values = Vec::with_capacity(k);

    for i in 0..k {
        values.push(a + (i as f64) * step_size);
    }

    return values;
}


/// Fill Brownian motion values on the grid between indices
/// 'start' and 'end' recursively, using the independent
/// increments property, with variance scaled by sigma^2.
fn recursive_fill(time: &[f64], path: &mut [f64], start: usize,
                  end: usize, sigma: f64, rng: &mut impl Rng) {

    // base condition: no more intermediate points to fill
    if end - start <= 1 {
        return;
    }

    // determine the midpoint index.
    let mid = (start + end) / 2;

    let t_1 = time[start];
    let t = time[mid];
    let t_2 = time[end];

    let a = path[start];
    let b = path[end];

    // compute the conditional mean and variance for the midpoint
    // you find the eq. on https://en.wikipedia.org/wiki/Brownian_bridge
    let mean = a + ((t - t_1) / (t_2 - t_1)) * (b - a);
    let var = sigma * sigma * (t_2 - t) * (t - t_1) / (t_2 - t_1);

    // sample the midpoint value from a normal distribution
    let normal = Normal::new(mean, var.sqrt()).unwrap();
    let sample = normal.sample(rng);
    path[mid] = sample;

    // recursively fill the left and right intervals
    recursive_fill(time, path, start, mid, sigma, rng);
    recursive_fill(time, path, mid, end, sigma, rng);
}


/// Generate a 'k' higher resolution Brownian motion path in
/// time partition {a, b}.
fn brownian_bridge(k: usize, a: f64, b: f64, sigma: f64, start_val: f64,
                   end_val: f64, rng: &mut impl Rng) -> (Vec<f64>, Vec<f64>) {
    let num_points = k + 2;
    let time_partition = linspace(a, b, num_points);
    let mut path = vec![0.0; num_points];

    path[0] = start_val;
    path[num_points - 1] = end_val;

    recursive_fill(&time_partition, &mut path, 0, num_points - 1, sigma, rng);
    return (time_partition, path);
}


/// Given a time partition with more than 1 interval: {a, b,
/// c, ..., n}, generate a 'k' higher resolution Brownian
/// motion path in-between all n-1 intervals.
fn refine_all_increments(time: &[f64], path: &[f64], k: usize, sigma: f64,
                         rng: &mut impl Rng) -> (Vec<f64>, Vec<f64>) {

    let mut new_time = Vec::new();
    let mut new_path = Vec::new();
    let num_increments = time.len() - 1;

    // in each interval, build a Brownian bridge
    for i in 0..num_increments {
        let (seg_time, seg_path) = brownian_bridge(k, time[i], time[i + 1],
                                                   sigma, path[i], path[i + 1], rng);

        // for the first segment include all points: [t_0, t_i],
        // for next segments skip the first point i.e: (t_i, t_j]
        // to avoid duplicating the common endpoint 't_i'
        if i == 0 {
            new_time.extend(seg_time);
            new_path.extend(seg_path);
        }
        else {
            new_time.extend(seg_time.into_iter().skip(1));
            new_path.extend(seg_path.into_iter().skip(1));
        }
    }

    return (new_time, new_path);
}


/// Refines the Brownian motion path using a vector of
/// incremental resolution improvement values 'k'.
fn incremental_refinment(resolution_increases: &[usize], a: f64, b: f64,
                         sigma: f64, start_val: f64, end_val: f64,
                         rng: &mut impl Rng)
                         -> Vec<(usize, Vec<f64>, Vec<f64>)> {

    // store resoluts here
    let mut refinements = Vec::new();

    // start with the original time partition {a, b}
    let mut current_time = vec![a, b];
    let mut current_path = vec![start_val, end_val];

    // refine path incrementally
    for &k in resolution_increases {
        let (new_time, new_path) = refine_all_increments(&current_time,
                                   &current_path, k, sigma, rng);

        refinements.push((k, new_time.clone(), new_path.clone()));

        // update current path for next refinement
        current_time = new_time;
        current_path = new_path;
    }

    return refinements;
}


// Calculate sample mean.
pub fn arithmetic_mean(sample: &[f64]) -> f64 {
    let mut sum = 0.0;
    let n = sample.len();

    for i in 0..n {
        sum += sample[i];
    }

    return sum / (n as f64);
}


/// Calculate sample standard deviation.
pub fn standard_deviation(sample: &[f64]) -> f64 {
    let mut sum_sq = 0.0;
    let n = sample.len();
    let mean = arithmetic_mean(sample);

    for i in 0..n {
        let diff = sample[i] - mean;
        sum_sq += diff * diff;
    }

    return (sum_sq / ((n - 1) as f64)).sqrt();
}


/// Compute the standard deviation of the Brownian motion
/// increments, to show that, in theory, σ_measured converges
/// to the real σ.
fn increment_stddev(time: &[f64], path: &[f64]) -> f64 {
    // calculate the increments between successive points
    let mut increments = Vec::with_capacity(path.len() - 1);

    for i in 0..(path.len() - 1) {
        increments.push(path[i + 1] - path[i]);
    }

    // assumes uniform partition (which it is when calling
    // 'refine_all_increments'), thus dt is constant
    let dt = time[1] - time[0];
    let incremental_stddev = standard_deviation(&increments);

    // scale by 1/sqrt(dt) because increments ~ N(0, σ²dt)
    incremental_stddev / dt.sqrt()
}


fn main() {
    let mut rng = rng();

    let a: f64 = 0.0;
    let b: f64 = 2.0;
    let sigma = 14.0;


    //let res_increase = vec![1, 1, 1, 1, 1, 10, 10];
    //let res_increase = vec![1, 2, 4, 8];
    let res_increase = vec![3, 3, 3];

    // sample starting and ending values
    let normal = Normal::new(0.0, sigma * (b - a).sqrt()).unwrap();
    let start_val = normal.sample(&mut rng);
    let end_val = normal.sample(&mut rng);

    // call recursive_bridge to get refined partitions
    let results = incremental_refinment(&res_increase, a, b, sigma,
                                        start_val, end_val, &mut rng);

    figure_2(&results);

    for result in &results {
        let (_, time, path) = result;
        let n = path.len() - 1;

        let stddev = increment_stddev(&time, &path);

        println!("num increments = {}, \t σ_measured = {}, \t σ_real = {}",
                 n, stddev, sigma);

    }

    let _ = input("[ENTER] to exit");
}
