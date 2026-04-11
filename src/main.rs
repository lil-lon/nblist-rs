use nblist::types::{System, UnitCell, Vector3};

fn main() {
    let sys = System {
        pos: vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(2.0, 0.0, 0.0),
            Vector3::new(5.0, 5.0, 5.0),
            Vector3::new(8.0, 0.0, 0.0),
        ],
        cell: UnitCell::new(
            Vector3::new(6.0, 1.0, 0.5),
            Vector3::new(0.5, 8.0, 1.5),
            Vector3::new(1.0, 0.5, 5.0),
        )
        .unwrap(),
        pbc: [true, true, true],
    };
    let cutoff = 3.0;

    let naive_result = nblist::naive::build_neighbor_list(&sys, cutoff, false);
    println!("Naive: {} neighbors", naive_result.len());

    let cell_list_result = nblist::cell_list::build_neighbor_list(&sys, cutoff, false);
    println!("Cell list: {} neighbors", cell_list_result.len());

    let cell_list_half = nblist::cell_list::build_neighbor_list(&sys, cutoff, true);
    println!("Cell list (half): {} neighbors", cell_list_half.len());
}
