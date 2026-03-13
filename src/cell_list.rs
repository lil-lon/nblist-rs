use crate::types::{Neighbor, System, UnitCell, Vector3};

pub fn build_neighbor_list(system: &System, cutoff: f64) -> Vec<Neighbor> {
    if cutoff == 0.0 {
        return Vec::new();
    }
    // Wrap positions into [0, 1) using fractional coordinates
    let wrapped: Vec<Vector3> = system
        .pos
        .iter()
        .map(|p| {
            let f = system.cell.get_fractional(p);
            Vector3::new(
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
            )
        })
        .collect();

    let (bins, n_bins, bin_sizes) = assign_to_bins(&wrapped, &system.cell, cutoff);
    let heights = system.cell.get_heights();
    let bin_heights: [f64; 3] = std::array::from_fn(|i| heights[i] * bin_sizes[i]);
    // When bin_size is smaller than cutoff, need to search more than 1 (when cutoff is greater than lattice constant)
    let search_ranges: [i32; 3] = std::array::from_fn(|i| {
        if system.pbc[i] {
            (cutoff / bin_heights[i]).ceil() as i32
        } else {
            0
        }
    });

    let mut result = Vec::new();
    for bz in 0..n_bins[2] as i32 {
        for by in 0..n_bins[1] as i32 {
            for bx in 0..n_bins[0] as i32 {
                let cur_idx = bx + n_bins[0] as i32 * (by + n_bins[1] as i32 * bz);
                let i_idxs: &Vec<usize> = &bins[cur_idx as usize];

                for dz in -search_ranges[2]..=search_ranges[2] {
                    for dy in -search_ranges[1]..=search_ranges[1] {
                        for dx in -search_ranges[0]..=search_ranges[0] {
                            for &i_idx in i_idxs {
                                let offset = [
                                    (bx + dx).div_euclid(n_bins[0] as i32),
                                    (by + dy).div_euclid(n_bins[1] as i32),
                                    (bz + dz).div_euclid(n_bins[2] as i32),
                                ];

                                let nb_bins = [
                                    (bx + dx).rem_euclid(n_bins[0] as i32),
                                    (by + dy).rem_euclid(n_bins[1] as i32),
                                    (bz + dz).rem_euclid(n_bins[2] as i32),
                                ];
                                let j_idx = nb_bins[0]
                                    + n_bins[0] as i32
                                        * (nb_bins[1] + n_bins[1] as i32 * nb_bins[2]);
                                let j_idxs: &Vec<usize> = &bins[j_idx as usize];
                                for &j_idx in j_idxs {
                                    if i_idx == j_idx && offset == [0, 0, 0] {
                                        continue;
                                    }
                                    let diff = wrapped[j_idx]
                                        + Vector3::new(
                                            offset[0] as f64,
                                            offset[1] as f64,
                                            offset[2] as f64,
                                        )
                                        - wrapped[i_idx];
                                    let distance = system.cell.get_cartesian(&diff).norm();
                                    if distance <= cutoff {
                                        result.push(Neighbor {
                                            i: i_idx,
                                            j: j_idx,
                                            offset,
                                            distance,
                                        })
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    result
}

/// Expects positions already wrapped into [0, 1).
/// Wrapping is the caller's responsibility.
fn assign_to_bins(
    wrapped: &[Vector3],
    cell: &UnitCell,
    cutoff: f64,
) -> (Vec<Vec<usize>>, [usize; 3], [f64; 3]) {
    fn get_bin_idx(pos: &Vector3, bin_sizes: [f64; 3], n_bins: [usize; 3]) -> usize {
        let x_bin_idx = (pos.x / bin_sizes[0]).floor().min((n_bins[0] - 1) as f64) as usize;
        let y_bin_idx = (pos.y / bin_sizes[1]).floor().min((n_bins[1] - 1) as f64) as usize;
        let z_bin_idx = (pos.z / bin_sizes[2]).floor().min((n_bins[2] - 1) as f64) as usize;
        x_bin_idx + n_bins[0] * (y_bin_idx + n_bins[1] * z_bin_idx)
    }

    let heights = cell.get_heights();
    let n_bins = heights.map(|h| (h / cutoff).floor().max(1.0) as usize);
    // fractional bin sizes
    let bin_sizes = n_bins.map(|bin| 1.0 / bin as f64);

    let atom_bins: Vec<usize> = wrapped
        .iter()
        .map(|p| get_bin_idx(p, bin_sizes, n_bins))
        .collect();

    let total_bins = n_bins[0] * n_bins[1] * n_bins[2];
    let mut bins: Vec<Vec<usize>> = vec![vec![]; total_bins];
    for (atom_idx, &bin) in atom_bins.iter().enumerate() {
        bins[bin].push(atom_idx);
    }
    (bins, n_bins, bin_sizes)
}

#[cfg(test)]
mod tests {
    use super::*;
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

    // --- assign_to_bins tests (cell list specific) ---

    #[test]
    fn test_assign_to_bins_basic() {
        // Non-cubic cell: a=6, b=10, c=15, cutoff=3
        // heights = [6, 10, 15] (orthogonal)
        // n_bins = [2, 3, 5], bin_sizes(frac) = [1/2, 1/3, 1/5]
        // linear index = bx + 2*(by + 3*bz)
        // Atom0 frac(0,0,0)              → bin(0,0,0) → 0
        // Atom1 frac(2/3, 0.5, 13/15)    → bin(1,1,4) → 1+2*(1+3*4) = 27
        // Atom2 frac(1/6, 0.8, 7/15)     → bin(0,2,2) → 0+2*(2+3*2) = 16
        // Atom3 frac(11/12, 0.05, 14/15) → bin(1,0,4) → 1+2*(0+3*4) = 25
        // Atom4 frac(1/12, 0.05, 1/30)   → bin(0,0,0) → 0 (same as Atom0)
        let wrapped = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(2.0 / 3.0, 0.5, 13.0 / 15.0),
            Vector3::new(1.0 / 6.0, 0.8, 7.0 / 15.0),
            Vector3::new(11.0 / 12.0, 0.05, 14.0 / 15.0),
            Vector3::new(1.0 / 12.0, 0.05, 1.0 / 30.0),
        ];
        let cell = UnitCell::new(
            Vector3::new(6.0, 0.0, 0.0),
            Vector3::new(0.0, 10.0, 0.0),
            Vector3::new(0.0, 0.0, 15.0),
        )
        .unwrap();
        let (bins, n_bins, bin_sizes) = assign_to_bins(&wrapped, &cell, 3.0);
        assert_eq!(n_bins, [2, 3, 5]);
        assert!((bin_sizes[0] - 0.5).abs() < 1e-10);
        assert!((bin_sizes[1] - 1.0 / 3.0).abs() < 1e-10);
        assert!((bin_sizes[2] - 0.2).abs() < 1e-10);
        assert_eq!(bins[0], vec![0, 4]);
        assert_eq!(bins[27], vec![1]);
        assert_eq!(bins[16], vec![2]);
        assert_eq!(bins[25], vec![3]);
    }

    #[test]
    fn test_assign_to_bins_same_bin() {
        // Two atoms close together → both in bin(0,0,0)
        // L=10, cutoff=3 → n_bins=[3,3,3], bin_sizes=[1/3,1/3,1/3]
        let wrapped = vec![Vector3::new(0.1, 0.1, 0.1), Vector3::new(0.2, 0.2, 0.2)];
        let cell = UnitCell::new(
            Vector3::new(10.0, 0.0, 0.0),
            Vector3::new(0.0, 10.0, 0.0),
            Vector3::new(0.0, 0.0, 10.0),
        )
        .unwrap();
        let (bins, _, _) = assign_to_bins(&wrapped, &cell, 3.0);
        assert_eq!(bins[0], vec![0, 1]);
    }

    #[test]
    fn test_assign_to_bins_cutoff_larger_than_cell() {
        // L=3, cutoff=5 → n_bins=[1,1,1], single bin
        let wrapped = vec![
            Vector3::new(1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0),
            Vector3::new(2.0 / 3.0, 1.0 / 3.0, 5.0 / 6.0),
        ];
        let cell = UnitCell::new(
            Vector3::new(3.0, 0.0, 0.0),
            Vector3::new(0.0, 3.0, 0.0),
            Vector3::new(0.0, 0.0, 3.0),
        )
        .unwrap();
        let (bins, n_bins, _) = assign_to_bins(&wrapped, &cell, 5.0);
        assert_eq!(n_bins, [1, 1, 1]);
        assert_eq!(bins.len(), 1);
        assert_eq!(bins[0], vec![0, 1]);
    }

    #[test]
    fn test_assign_to_bins_triclinic() {
        // Triclinic cell: a=(5,0,0), b=(4,3,0), c=(0,0,5), cutoff=2
        // heights: d_a=3.0, d_b=3.0, d_c=5.0  (differ from norms |a|=5, |b|=5, |c|=5)
        // n_bins = [floor(3/2)=1, floor(3/2)=1, floor(5/2)=2]
        //   (norm-based would give [2, 2, 2] — wrong)
        // bin_sizes(frac) = [1.0, 1.0, 0.5]
        // linear index = bx + 1*(by + 1*bz) = bz
        // Atom0 frac(0.1, 0.2, 0.1) → bin(0,0,0) → 0
        // Atom1 frac(0.5, 0.8, 0.7) → bin(0,0,1) → 1
        // Atom2 frac(0.9, 0.9, 0.3) → bin(0,0,0) → 0
        let wrapped = vec![
            Vector3::new(0.1, 0.2, 0.1),
            Vector3::new(0.5, 0.8, 0.7),
            Vector3::new(0.9, 0.9, 0.3),
        ];
        let cell = UnitCell::new(
            Vector3::new(5.0, 0.0, 0.0),
            Vector3::new(4.0, 3.0, 0.0),
            Vector3::new(0.0, 0.0, 5.0),
        )
        .unwrap();
        let (bins, n_bins, bin_sizes) = assign_to_bins(&wrapped, &cell, 2.0);
        assert_eq!(n_bins, [1, 1, 2]);
        assert!((bin_sizes[0] - 1.0).abs() < 1e-10);
        assert!((bin_sizes[1] - 1.0).abs() < 1e-10);
        assert!((bin_sizes[2] - 0.5).abs() < 1e-10);
        assert_eq!(bins[0], vec![0, 2]);
        assert_eq!(bins[1], vec![1]);
    }

    #[test]
    fn test_assign_to_bins_boundary() {
        // L=10, cutoff=5 → n_bins=[2,2,2], bin_size(frac)=0.5
        // frac 0.4999 → bin 0, frac 0.5 → bin 1, frac 0.9999 → bin 1 (clamped)
        let wrapped = vec![
            Vector3::new(0.4999, 0.0, 0.0),
            Vector3::new(0.5, 0.0, 0.0),
            Vector3::new(0.9999, 0.0, 0.0),
        ];
        let cell = UnitCell::new(
            Vector3::new(10.0, 0.0, 0.0),
            Vector3::new(0.0, 10.0, 0.0),
            Vector3::new(0.0, 0.0, 10.0),
        )
        .unwrap();
        let (bins, _, _) = assign_to_bins(&wrapped, &cell, 5.0);
        assert_eq!(bins[0], vec![0]);
        assert_eq!(bins[1], vec![1, 2]);
    }
}
