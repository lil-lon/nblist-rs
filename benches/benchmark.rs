use std::fs;
use std::io::Write;
use std::path::PathBuf;
use std::time::Instant;

use nblist::types::{Neighbor, System, UnitCell, Vector3};

const LATTICE_CONST: f64 = 5.43;
const CUTOFF: f64 = 6.0;
const SIZES: &[usize] = &[2, 3, 5, 8, 10, 15, 20];
const N_RUNS: usize = 5;

fn mean_std(values: &[f64]) -> (f64, f64) {
    let n = values.len() as f64;
    let mean = values.iter().sum::<f64>() / n;
    let var = values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / n;
    (mean, var.sqrt())
}

fn simple_cubic(n: usize, a: f64) -> System {
    let mut pos = Vec::with_capacity(n * n * n);
    for iz in 0..n {
        for iy in 0..n {
            for ix in 0..n {
                pos.push(Vector3::new(ix as f64 * a, iy as f64 * a, iz as f64 * a));
            }
        }
    }
    let l = n as f64 * a;
    System {
        pos,
        cell: UnitCell::new(
            Vector3::new(l, 0.0, 0.0),
            Vector3::new(0.0, l, 0.0),
            Vector3::new(0.0, 0.0, l),
        )
        .unwrap(),
        pbc: [true, true, true],
    }
}

fn bench(
    system: &System,
    cutoff: f64,
    method: fn(&System, f64, bool) -> Vec<Neighbor>,
) -> (f64, f64) {
    // Warmup
    let _ = method(system, cutoff, false);

    let mut times = Vec::with_capacity(N_RUNS);
    for _ in 0..N_RUNS {
        let start = Instant::now();
        let _ = method(system, cutoff, false);
        times.push(start.elapsed().as_secs_f64());
    }
    mean_std(&times)
}

fn main() {
    let mut csv = String::from("n_atoms,method,mean_s,std_s\n");

    for &n in SIZES {
        let n_atoms = n * n * n;
        let system = simple_cubic(n, LATTICE_CONST);
        print!("N={n:>2} ({n_atoms:>5} atoms) ... ");
        std::io::stdout().flush().unwrap();

        let (mean, std) = bench(&system, CUTOFF, nblist::naive::build_neighbor_list);
        print!("naive {mean:.4e}s (±{std:.4e})  ");
        csv.push_str(&format!("{n_atoms},naive,{mean:.6e},{std:.6e}\n"));

        let (mean, std) = bench(&system, CUTOFF, nblist::cell_list::build_neighbor_list);
        println!("cell_list {mean:.4e}s (±{std:.4e})");
        csv.push_str(&format!("{n_atoms},cell_list,{mean:.6e},{std:.6e}\n"));
    }

    let out_path =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("benchmarks/results/bench_rust_results.csv");
    fs::write(&out_path, &csv).expect("failed to write CSV");
    println!("\nResults written to {}", out_path.display());
}
