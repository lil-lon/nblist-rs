use crate::types::{Neighbor, System, UnitCell, Vector3};

pub type BuildFn = fn(&System, f64) -> Vec<Neighbor>;

pub fn test_basic(build: BuildFn) {
    // 4 atoms in a 10x10x10 box, cutoff = 3.0
    // A(0,0,0) B(2,0,0) C(5,5,5) D(8,0,0)
    // A-B: 2.0 ✓, B-A: 2.0 ✓
    // A-D via PBC: 2.0 ✓, D-A via PBC: 2.0 ✓
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
    let result = build(&sys, 3.0);
    // Both directions: (0,1), (1,0), (0,3), (3,0)
    let pairs: Vec<(usize, usize)> = result.iter().map(|n| (n.i, n.j)).collect();
    assert!(pairs.contains(&(0, 1)));
    assert!(pairs.contains(&(1, 0)));
    assert!(pairs.contains(&(0, 3)));
    assert!(pairs.contains(&(3, 0)));
    assert_eq!(result.len(), 4);
}

pub fn test_no_neighbors(build: BuildFn) {
    let sys = System {
        pos: vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(5.0, 5.0, 5.0)],
        cell: UnitCell {
            a: 20.0,
            b: 20.0,
            c: 20.0,
        },
    };
    let result = build(&sys, 1.0);
    assert!(result.is_empty());
}

pub fn test_pbc_corner(build: BuildFn) {
    // A(0.5, 0.5, 0.5) B(9.5, 9.5, 9.5) in 10x10x10 box
    // dist via PBC = sqrt(3) ≈ 1.73
    let sys = System {
        pos: vec![Vector3::new(0.5, 0.5, 0.5), Vector3::new(9.5, 9.5, 9.5)],
        cell: UnitCell {
            a: 10.0,
            b: 10.0,
            c: 10.0,
        },
    };
    let result = build(&sys, 2.0);
    // Both directions
    let pairs: Vec<(usize, usize)> = result.iter().map(|n| (n.i, n.j)).collect();
    assert_eq!(pairs.len(), 2);
    assert!(pairs.contains(&(0, 1)));
    assert!(pairs.contains(&(1, 0)));
}

pub fn test_self_image(build: BuildFn) {
    // Single atom at origin, L=3, cutoff=3.0 (== L, boundary case)
    // Self-images at distance 3.0 along each axis (±a, ±b, ±c) → 6 neighbors
    // Diagonal images e.g. (3,3,0) at dist=sqrt(18)≈4.24 → outside cutoff
    let sys = System {
        pos: vec![Vector3::new(0.0, 0.0, 0.0)],
        cell: UnitCell {
            a: 3.0,
            b: 3.0,
            c: 3.0,
        },
    };
    let result = build(&sys, 3.0);
    // All neighbors are self-images (i==0, j==0)
    assert!(result.iter().all(|n| n.i == 0 && n.j == 0));
    assert_eq!(result.len(), 6);
    // Each should have distance 3.0
    assert!(result.iter().all(|n| (n.distance - 3.0).abs() < 1e-10));
    // Offsets should be the 6 face-adjacent cells
    let offsets: Vec<[i32; 3]> = result.iter().map(|n| n.offset).collect();
    assert!(offsets.contains(&[1, 0, 0]));
    assert!(offsets.contains(&[-1, 0, 0]));
    assert!(offsets.contains(&[0, 1, 0]));
    assert!(offsets.contains(&[0, -1, 0]));
    assert!(offsets.contains(&[0, 0, 1]));
    assert!(offsets.contains(&[0, 0, -1]));
}

pub fn test_multiple_images(build: BuildFn) {
    // 2 atoms in L=4 box, cutoff=5
    // A(0,0,0) B(1,0,0)
    let sys = System {
        pos: vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 0.0, 0.0)],
        cell: UnitCell {
            a: 4.0,
            b: 4.0,
            c: 4.0,
        },
    };
    let result = build(&sys, 5.0);

    // Check A→B has multiple offsets
    let a_to_b: Vec<&Neighbor> = result.iter().filter(|n| n.i == 0 && n.j == 1).collect();
    assert_eq!(a_to_b.len(), 11);
    assert!(
        a_to_b
            .iter()
            .any(|n| n.offset == [0, 0, 0] && (n.distance - 1.0).abs() < 1e-10)
    );
    assert!(
        a_to_b
            .iter()
            .any(|n| n.offset == [-1, 0, 0] && (n.distance - 3.0).abs() < 1e-10)
    );

    // Check self-images exist (A→A with offset != [0,0,0])
    let a_self: Vec<&Neighbor> = result.iter().filter(|n| n.i == 0 && n.j == 0).collect();
    assert_eq!(a_self.len(), 6);
    assert!(a_self.iter().all(|n| (n.distance - 4.0).abs() < 1e-10));
}

pub fn test_non_cubic(build: BuildFn) {
    // Non-cubic cell: a=3, b=10, c=10, cutoff=3.5
    // A(0,0,0) B(1,0,0)
    let sys = System {
        pos: vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 0.0, 0.0)],
        cell: UnitCell {
            a: 3.0,
            b: 10.0,
            c: 10.0,
        },
    };
    let result = build(&sys, 3.5);

    // A→B: offset [0,0,0] (dist=1) and [-1,0,0] (dist=2)
    let a_to_b: Vec<&Neighbor> = result.iter().filter(|n| n.i == 0 && n.j == 1).collect();
    assert_eq!(a_to_b.len(), 2);
    assert!(
        a_to_b
            .iter()
            .any(|n| n.offset == [0, 0, 0] && (n.distance - 1.0).abs() < 1e-10)
    );
    assert!(
        a_to_b
            .iter()
            .any(|n| n.offset == [-1, 0, 0] && (n.distance - 2.0).abs() < 1e-10)
    );

    // A self-images: only along a-axis (dist=3), not b or c (dist=10)
    let a_self: Vec<&Neighbor> = result.iter().filter(|n| n.i == 0 && n.j == 0).collect();
    assert_eq!(a_self.len(), 2);
    assert!(a_self.iter().any(|n| n.offset == [1, 0, 0]));
    assert!(a_self.iter().any(|n| n.offset == [-1, 0, 0]));
}

pub fn test_unwrapped_coordinates(build: BuildFn) {
    // Atom j is far outside the primary cell (diff > L).
    // L=5, cutoff=3, i=(1,1,1), j=(14,-8,22)
    // j wrapped = (4, 2, 2), minimum-image dist = sqrt(9+1+1) = sqrt(11) ≈ 3.32 → outside
    // But i=(1,1,1) j_wrapped=(4,2,2): diff=(3,1,1)
    //   offset [0,0,0]: dist=sqrt(11)≈3.32 → outside cutoff
    //   offset [-1,0,0]: j at (-1,2,2), diff=(-2,1,1), dist=sqrt(6)≈2.45 ✓
    let sys = System {
        pos: vec![Vector3::new(1.0, 1.0, 1.0), Vector3::new(14.0, -8.0, 22.0)],
        cell: UnitCell {
            a: 5.0,
            b: 5.0,
            c: 5.0,
        },
    };
    let result = build(&sys, 3.0);
    assert_eq!(result.len(), 2); // i→j and j→i
    let i_to_j: Vec<&Neighbor> = result.iter().filter(|n| n.i == 0 && n.j == 1).collect();
    assert_eq!(i_to_j.len(), 1);
    assert_eq!(i_to_j[0].offset, [-1, 0, 0]);
    assert!((i_to_j[0].distance - 6.0_f64.sqrt()).abs() < 1e-10);
}
