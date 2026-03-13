use crate::types::{Neighbor, System, Vector3};

pub fn build_neighbor_list(system: &System, cutoff: f64) -> Vec<Neighbor> {
    let mut result = Vec::new();
    let heights = system.cell.get_heights();
    let a_replicas = if system.pbc[0] {
        (cutoff / heights[0]).ceil() as i32
    } else {
        0
    };
    let b_replicas = if system.pbc[1] {
        (cutoff / heights[1]).ceil() as i32
    } else {
        0
    };
    let c_replicas = if system.pbc[2] {
        (cutoff / heights[2]).ceil() as i32
    } else {
        0
    };
    // Wrap positions into [0, 1) so that ceil(cutoff / L) replicas suffice.
    let wrapped: Vec<Vector3> = system
        .pos
        .iter()
        .map(|p| {
            let f = system.cell.get_fractional(p);
            let f_wrapped = Vector3::new(
                if system.pbc[0] {
                    f.x.rem_euclid(1.0)
                } else {
                    f.x
                },
                if system.pbc[1] {
                    f.y.rem_euclid(1.0)
                } else {
                    f.y
                },
                if system.pbc[2] {
                    f.z.rem_euclid(1.0)
                } else {
                    f.z
                },
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
    fn orthorhombic() {
        test_neighbors::test_orthorhombic(super::build_neighbor_list);
    }

    #[test]
    fn unwrapped_coordinates() {
        test_neighbors::test_unwrapped_coordinates(super::build_neighbor_list);
    }

    #[test]
    fn monoclinic_basic() {
        test_neighbors::test_monoclinic_basic(super::build_neighbor_list);
    }

    #[test]
    fn monoclinic_pbc_corner() {
        test_neighbors::test_monoclinic_pbc_corner(super::build_neighbor_list);
    }

    #[test]
    fn monoclinic_unwrapped() {
        test_neighbors::test_monoclinic_unwrapped(super::build_neighbor_list);
    }

    #[test]
    fn monoclinic_self_image() {
        test_neighbors::test_monoclinic_self_image(super::build_neighbor_list);
    }

    #[test]
    fn monoclinic_skewed() {
        test_neighbors::test_monoclinic_skewed(super::build_neighbor_list);
    }

    #[test]
    fn triclinic_full_3d() {
        test_neighbors::test_triclinic_full_3d(super::build_neighbor_list);
    }

    #[test]
    fn triclinic_full_3d_self_image() {
        test_neighbors::test_triclinic_full_3d_self_image(super::build_neighbor_list);
    }

    #[test]
    fn monoclinic_boundary_distance() {
        test_neighbors::test_monoclinic_boundary_distance(super::build_neighbor_list);
    }

    #[test]
    fn monoclinic_left_handed() {
        test_neighbors::test_monoclinic_left_handed(super::build_neighbor_list);
    }

    #[test]
    fn triclinic_multi_height_replicas() {
        test_neighbors::test_triclinic_multi_height_replicas(super::build_neighbor_list);
    }

    #[test]
    fn triclinic_multi_atom() {
        test_neighbors::test_triclinic_multi_atom(super::build_neighbor_list);
    }

    #[test]
    fn triclinic_near_degenerate() {
        test_neighbors::test_triclinic_near_degenerate(super::build_neighbor_list);
    }

    #[test]
    fn slab_self_image() {
        test_neighbors::test_slab_self_image(super::build_neighbor_list);
    }

    #[test]
    fn wire_self_image() {
        test_neighbors::test_wire_self_image(super::build_neighbor_list);
    }

    #[test]
    fn isolated_no_images() {
        test_neighbors::test_isolated_no_images(super::build_neighbor_list);
    }

    #[test]
    fn slab_two_atoms() {
        test_neighbors::test_slab_two_atoms(super::build_neighbor_list);
    }

    #[test]
    fn isolated_two_atoms() {
        test_neighbors::test_isolated_two_atoms(super::build_neighbor_list);
    }

    #[test]
    fn slab_no_wrap_non_periodic() {
        test_neighbors::test_slab_no_wrap_non_periodic(super::build_neighbor_list);
    }
}
