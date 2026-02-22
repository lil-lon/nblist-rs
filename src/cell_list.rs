use crate::types::{Neighbor, System, UnitCell, Vector3};

pub fn build_neighbor_list(system: &System, cutoff: f64) -> Vec<Neighbor> {
    if cutoff == 0.0 {
        return Vec::new();
    }
    // Wrap positions into [0, L) so that ceil(cutoff / L) replicas suffice.
    let wrapped: Vec<Vector3> = system
        .pos
        .iter()
        .map(|p| {
            Vector3::new(
                p.x.rem_euclid(system.cell.a),
                p.y.rem_euclid(system.cell.b),
                p.z.rem_euclid(system.cell.c),
            )
        })
        .collect();

    let (bins, n_bins, bin_sizes) = assign_to_bins(&wrapped, &system.cell, cutoff);
    // When bin_size is smaller than cutoff, need to search more than 1 (when cutoff is greater than lattice constant)
    let search_ranges: [i32; 3] = bin_sizes.map(|b| (cutoff / b).ceil() as i32);

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
                                let shift_vec = Vector3::new(
                                    offset[0] as f64 * system.cell.a,
                                    offset[1] as f64 * system.cell.b,
                                    offset[2] as f64 * system.cell.c,
                                );

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
                                    let distance =
                                        (wrapped[j_idx] + shift_vec - wrapped[i_idx]).norm();
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

/// Expects positions already wrapped into [0, L).
/// Wrapping is the caller's responsibility.
pub(crate) fn assign_to_bins(
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

    let lattice: [f64; 3] = [cell.a, cell.b, cell.c];
    let n_bins: [usize; 3] = lattice.map(|l| (l / cutoff).floor().max(1.0) as usize);
    let bin_sizes: [f64; 3] = std::array::from_fn(|i| lattice[i] / n_bins[i] as f64);

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
    fn non_cubic() {
        test_neighbors::test_non_cubic(super::build_neighbor_list);
    }

    #[test]
    fn unwrapped_coordinates() {
        test_neighbors::test_unwrapped_coordinates(super::build_neighbor_list);
    }

    // --- assign_to_bins tests (cell list specific) ---

    #[test]
    fn test_assign_to_bins_basic() {
        // Non-cubic cell: a=6, b=10, c=15, cutoff=3
        // n_bins = [2, 3, 5], bin_sizes = [3.0, 10/3, 3.0]
        // linear index = bx + 2*(by + 3*bz)
        // Atom0(0,0,0)     → bin(0,0,0) → 0
        // Atom1(4,5,13)    → bin(1,1,4) → 1+2*(1+3*4) = 27
        // Atom2(1,8,7)     → bin(0,2,2) → 0+2*(2+3*2) = 16
        // Atom3(5.5,0.5,14)→ bin(1,0,4) → 1+2*(0+3*4) = 25
        // Atom4(0.5,0.5,0.5)→ bin(0,0,0) → 0 (same as Atom0)
        let wrapped = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(4.0, 5.0, 13.0),
            Vector3::new(1.0, 8.0, 7.0),
            Vector3::new(5.5, 0.5, 14.0),
            Vector3::new(0.5, 0.5, 0.5),
        ];
        let cell = UnitCell {
            a: 6.0,
            b: 10.0,
            c: 15.0,
        };
        let (bins, n_bins, bin_sizes) = assign_to_bins(&wrapped, &cell, 3.0);
        assert_eq!(n_bins, [2, 3, 5]);
        assert!((bin_sizes[0] - 3.0).abs() < 1e-10);
        assert!((bin_sizes[1] - 10.0 / 3.0).abs() < 1e-10);
        assert!((bin_sizes[2] - 3.0).abs() < 1e-10);
        assert_eq!(bins[0], vec![0, 4]);
        assert_eq!(bins[27], vec![1]);
        assert_eq!(bins[16], vec![2]);
        assert_eq!(bins[25], vec![3]);
    }

    #[test]
    fn test_assign_to_bins_same_bin() {
        // Two atoms close together → both in bin(0,0,0)
        let wrapped = vec![Vector3::new(1.0, 1.0, 1.0), Vector3::new(2.0, 2.0, 2.0)];
        let cell = UnitCell {
            a: 10.0,
            b: 10.0,
            c: 10.0,
        };
        let (bins, _, _) = assign_to_bins(&wrapped, &cell, 3.0);
        assert_eq!(bins[0], vec![0, 1]);
    }

    #[test]
    fn test_assign_to_bins_cutoff_larger_than_cell() {
        // L=3, cutoff=5 → n_bins=[1,1,1], single bin
        let wrapped = vec![Vector3::new(0.5, 0.5, 0.5), Vector3::new(2.0, 1.0, 2.5)];
        let cell = UnitCell {
            a: 3.0,
            b: 3.0,
            c: 3.0,
        };
        let (bins, n_bins, _) = assign_to_bins(&wrapped, &cell, 5.0);
        assert_eq!(n_bins, [1, 1, 1]);
        assert_eq!(bins.len(), 1);
        assert_eq!(bins[0], vec![0, 1]);
    }

    #[test]
    fn test_assign_to_bins_boundary() {
        // L=10, cutoff=5 → n_bins=[2,2,2], bin_size=5.0
        // Atom at 4.999 → bin 0, atom at 5.0 → bin 1, atom at 9.999 → bin 1 (clamped)
        let wrapped = vec![
            Vector3::new(4.999, 0.0, 0.0),
            Vector3::new(5.0, 0.0, 0.0),
            Vector3::new(9.999, 0.0, 0.0),
        ];
        let cell = UnitCell {
            a: 10.0,
            b: 10.0,
            c: 10.0,
        };
        let (bins, _, _) = assign_to_bins(&wrapped, &cell, 5.0);
        assert_eq!(bins[0], vec![0]);
        assert_eq!(bins[1], vec![1, 2]);
    }
}
