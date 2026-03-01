use crate::types::{Neighbor, System, Vector3};

pub fn build_neighbor_list(system: &System, cutoff: f64) -> Vec<Neighbor> {
    let mut result = Vec::new();
    let heights = system.cell.get_heights();
    let a_replicas = (cutoff / heights[0]).ceil() as i32;
    let b_replicas = (cutoff / heights[1]).ceil() as i32;
    let c_replicas = (cutoff / heights[2]).ceil() as i32;
    // Wrap positions into [0, 1) so that ceil(cutoff / L) replicas suffice.
    let wrapped: Vec<Vector3> = system
        .pos
        .iter()
        .map(|p| {
            let f = system.cell.get_fractional(p);
            let f_wrapped = Vector3::new(
                f.x.rem_euclid(1.0),
                f.y.rem_euclid(1.0),
                f.z.rem_euclid(1.0),
            );
            system.cell.get_cartesian(&f_wrapped)
        })
        .collect();
    let n = wrapped.len();

    for i in 0..n {
        for j in 0..n {
            for a_ind in -a_replicas..=a_replicas {
                for b_ind in -b_replicas..=b_replicas {
                    for c_ind in -c_replicas..=c_replicas {
                        if i == j && (a_ind, b_ind, c_ind) == (0, 0, 0) {
                            continue;
                        }
                        let pos_j = wrapped[j]
                            + system.cell.a * a_ind as f64
                            + system.cell.b * b_ind as f64
                            + system.cell.c * c_ind as f64;
                        let distance = (pos_j - wrapped[i]).norm();
                        if distance <= cutoff {
                            result.push(Neighbor {
                                i,
                                j,
                                offset: [a_ind, b_ind, c_ind],
                                distance,
                            })
                        }
                    }
                }
            }
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use crate::test_neighbors;

    #[test]
    fn basic() {
        test_neighbors::test_basic(super::build_neighbor_list);
    }

    #[test]
    fn no_neighbors() {
        test_neighbors::test_no_neighbors(super::build_neighbor_list);
    }

    #[test]
    fn pbc_corner() {
        test_neighbors::test_pbc_corner(super::build_neighbor_list);
    }

    #[test]
    fn self_image() {
        test_neighbors::test_self_image(super::build_neighbor_list);
    }

    #[test]
    fn multiple_images() {
        test_neighbors::test_multiple_images(super::build_neighbor_list);
    }

    #[test]
    fn non_cubic() {
        test_neighbors::test_non_cubic(super::build_neighbor_list);
    }

    #[test]
    fn unwrapped_coordinates() {
        test_neighbors::test_unwrapped_coordinates(super::build_neighbor_list);
    }
}
