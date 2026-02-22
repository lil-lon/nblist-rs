use nblist::types::{System, UnitCell, Vector3};

fn main() {
    let sys = System {
        pos: vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(2.0, 0.0, 0.0),
            Vector3::new(5.0, 5.0, 5.0),
            Vector3::new(8.0, 0.0, 0.0),
        ],
        cell: UnitCell {
            a: 10.0,
            b: 10.0,
            c: 10.0,
        },
    };
    let cutoff = 3.0;

    let naive_result = nblist::naive::build_neighbor_list(&sys, cutoff);
    println!("Naive: {} neighbors", naive_result.len());

    let cell_list_result = nblist::cell_list::build_neighbor_list(&sys, cutoff);
    println!("Cell list: {} neighbors", cell_list_result.len());
}
