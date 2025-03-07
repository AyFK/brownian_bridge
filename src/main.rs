use std::io::{self, Write};

use rand::rng;
use rand::Rng;
use rand_distr::{Distribution, Normal};

mod plots;
use plots::resolution::figure;

#[allow(dead_code)]
fn input(prompt: &str) -> String {
    print!("{}", prompt);
    io::stdout().flush().expect("Failed to flush stdout.");
    let mut buffer = String::new();
    io::stdin().read_line(&mut buffer).expect("Failed to read line.");
    buffer.trim().to_string()
}

/// Generate a uniform partition on the time interval
/// [a, b] with 'n' values, 'n-1' increments.
fn linspace(a: f64, b: f64, n: usize) -> Vec<f64> {
    let step_size = (b - a) / (n as f64 - 1.0);
    let mut values = Vec::with_capacity(n);

    for i in 0..n {
        values.push(a + (i as f64) * step_size);
    }

    return values;
}


/// Fill Brownian motion values on the grid between indices
/// 'start' and 'end' recursively, using the independent
/// increments property, with variance scaled by sigma^2.
fn recursive_fill(time_partition: &[f64], brownian: &mut [f64], start: usize,
                  end: usize, sigma: f64, rng: &mut impl Rng) {

    // base condition: no more intermediate points to fill
    if end - start <= 1 {
        return;
    }

    // determine the midpoint index.
    let mid = (start + end) / 2;

    let t_1 = time_partition[start];
    let t = time_partition[mid];
    let t_2 = time_partition[end];

    let a = brownian[start];
    let b = brownian[end];

    // compute the conditional mean and variance for the midpoint
    // you find the eq. on https://en.wikipedia.org/wiki/Brownian_bridge
    let mean = a + ((t - t_1) / (t_2 - t_1)) * (b - a);
    let var = sigma * sigma * (t_2 - t) * (t - t_1) / (t_2 - t_1);

    // sample the midpoint value from a normal distribution
    let normal = Normal::new(mean, var.sqrt()).unwrap();
    let sample = normal.sample(rng);
    brownian[mid] = sample;

    // recursively fill the left and right intervals
    recursive_fill(time_partition, brownian, start, mid, sigma, rng);
    recursive_fill(time_partition, brownian, mid, end, sigma, rng);

}


/// Generate a 'n' higher resolution Brownian motion path in
/// time interval [a, b].
fn brownian_bridge(n: usize, a: f64, b: f64, sigma: f64, start_val: f64,
                   end_val: f64, rng: &mut impl Rng) -> (Vec<f64>, Vec<f64>) {
    let num_points = n + 2;
    let time_partition = linspace(a, b, num_points);
    let mut brownian = vec![0.0; num_points];

    brownian[0] = start_val;
    brownian[num_points - 1] = end_val;

    recursive_fill(&time_partition, &mut brownian, 0, num_points - 1, sigma, rng);
    return (time_partition, brownian);
}


fn main() {

    let mut rng = rng();

    let a = 0.0;
    let b = 1.0;


    let sigma = 1.0;
    let extra_res = 5;

    let start_val = 0.0;
    let normal = Normal::new(0.0, sigma).unwrap();
    let end_val = normal.sample(&mut rng);


    let (t, w) = brownian_bridge(extra_res, a, b, sigma, start_val, end_val,
                                 &mut rng);

    figure(&[a, b], &[start_val, end_val], &t, &w);
    let _ = input("[ENTER] to exit");
}
