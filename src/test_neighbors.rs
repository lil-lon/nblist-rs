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
        cell: UnitCell::new(
            Vector3::new(10.0, 0.0, 0.0),
            Vector3::new(0.0, 10.0, 0.0),
            Vector3::new(0.0, 0.0, 10.0),
        )
        .unwrap(),
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
        cell: UnitCell::new(
            Vector3::new(20.0, 0.0, 0.0),
            Vector3::new(0.0, 20.0, 0.0),
            Vector3::new(0.0, 0.0, 20.0),
        )
        .unwrap(),
    };
    let result = build(&sys, 1.0);
    assert!(result.is_empty());
}

pub fn test_pbc_corner(build: BuildFn) {
    // A(0.5, 0.5, 0.5) B(9.5, 9.5, 9.5) in 10x10x10 box
    // dist via PBC = sqrt(3) ≈ 1.73
    let sys = System {
        pos: vec![Vector3::new(0.5, 0.5, 0.5), Vector3::new(9.5, 9.5, 9.5)],
        cell: UnitCell::new(
            Vector3::new(10.0, 0.0, 0.0),
            Vector3::new(0.0, 10.0, 0.0),
            Vector3::new(0.0, 0.0, 10.0),
        )
        .unwrap(),
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
        cell: UnitCell::new(
            Vector3::new(3.0, 0.0, 0.0),
            Vector3::new(0.0, 3.0, 0.0),
            Vector3::new(0.0, 0.0, 3.0),
        )
        .unwrap(),
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
        cell: UnitCell::new(
            Vector3::new(4.0, 0.0, 0.0),
            Vector3::new(0.0, 4.0, 0.0),
            Vector3::new(0.0, 0.0, 4.0),
        )
        .unwrap(),
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
        cell: UnitCell::new(
            Vector3::new(3.0, 0.0, 0.0),
            Vector3::new(0.0, 10.0, 0.0),
            Vector3::new(0.0, 0.0, 10.0),
        )
        .unwrap(),
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
        cell: UnitCell::new(
            Vector3::new(5.0, 0.0, 0.0),
            Vector3::new(0.0, 5.0, 0.0),
            Vector3::new(0.0, 0.0, 5.0),
        )
        .unwrap(),
    };
    let result = build(&sys, 3.0);
    assert_eq!(result.len(), 2); // i→j and j→i
    let i_to_j: Vec<&Neighbor> = result.iter().filter(|n| n.i == 0 && n.j == 1).collect();
    assert_eq!(i_to_j.len(), 1);
    assert_eq!(i_to_j[0].offset, [-1, 0, 0]);
    assert!((i_to_j[0].distance - 6.0_f64.sqrt()).abs() < 1e-10);
}

pub fn test_triclinic_basic(build: BuildFn) {
    // Triclinic cell: a=(5,0,0), b=(2,5,0), c=(0,0,5)
    // A=(1,1,0), B=(6,4,0), cutoff=2.9
    // Both already inside cell in fractional coords.
    // Nearest image via offset [-1,-1,0]:
    //   B + (-1)*a + (-1)*b - A = (6-5-2-1, 4-0-5-1, 0) = (-2,-2,0)
    //   dist = sqrt(8) ≈ 2.83 < 2.9 ✓
    // All other offsets give dist > 2.9
    let sys = System {
        pos: vec![Vector3::new(1.0, 1.0, 0.0), Vector3::new(6.0, 4.0, 0.0)],
        cell: UnitCell::new(
            Vector3::new(5.0, 0.0, 0.0),
            Vector3::new(2.0, 5.0, 0.0),
            Vector3::new(0.0, 0.0, 5.0),
        )
        .unwrap(),
    };
    let result = build(&sys, 2.9);
    assert_eq!(result.len(), 2);
    let a_to_b: Vec<&Neighbor> = result.iter().filter(|n| n.i == 0 && n.j == 1).collect();
    assert_eq!(a_to_b.len(), 1);
    assert_eq!(a_to_b[0].offset, [-1, -1, 0]);
    assert!((a_to_b[0].distance - 8.0_f64.sqrt()).abs() < 1e-10);
    let b_to_a: Vec<&Neighbor> = result.iter().filter(|n| n.i == 1 && n.j == 0).collect();
    assert_eq!(b_to_a.len(), 1);
    assert_eq!(b_to_a[0].offset, [1, 1, 0]);
    assert!((b_to_a[0].distance - 8.0_f64.sqrt()).abs() < 1e-10);
}

pub fn test_triclinic_pbc_corner(build: BuildFn) {
    // Triclinic cell: a=(5,0,0), b=(2,5,0), c=(0,0,5)
    // A at frac=(0.05,0.05,0.05), B at frac=(0.95,0.95,0.95)
    // A_cart = 0.05*(5,0,0)+0.05*(2,5,0)+0.05*(0,0,5) = (0.35, 0.25, 0.25)
    // B_cart = 0.95*(5,0,0)+0.95*(2,5,0)+0.95*(0,0,5) = (6.65, 4.75, 4.75)
    // Offset [-1,-1,-1]: B-a-b-c-A = (-0.7,-0.5,-0.5), dist=sqrt(0.99)
    let cell = UnitCell::new(
        Vector3::new(5.0, 0.0, 0.0),
        Vector3::new(2.0, 5.0, 0.0),
        Vector3::new(0.0, 0.0, 5.0),
    )
    .unwrap();
    let a_cart = cell.get_cartesian(&Vector3::new(0.05, 0.05, 0.05));
    let b_cart = cell.get_cartesian(&Vector3::new(0.95, 0.95, 0.95));
    let sys = System {
        pos: vec![a_cart, b_cart],
        cell,
    };
    let result = build(&sys, 1.5);
    assert_eq!(result.len(), 2);
    let a_to_b: Vec<&Neighbor> = result.iter().filter(|n| n.i == 0 && n.j == 1).collect();
    assert_eq!(a_to_b.len(), 1);
    assert_eq!(a_to_b[0].offset, [-1, -1, -1]);
    assert!((a_to_b[0].distance - 0.99_f64.sqrt()).abs() < 1e-10);
    let b_to_a: Vec<&Neighbor> = result.iter().filter(|n| n.i == 1 && n.j == 0).collect();
    assert_eq!(b_to_a.len(), 1);
    assert_eq!(b_to_a[0].offset, [1, 1, 1]);
    assert!((b_to_a[0].distance - 0.99_f64.sqrt()).abs() < 1e-10);
}

pub fn test_triclinic_unwrapped(build: BuildFn) {
    // Triclinic cell: a=(5,0,0), b=(2,5,0), c=(0,0,5)
    // A=(1,1,1), B=(14,6,1) — B is far outside the cell.
    // B_frac = (2.32, 1.2, 0.2) → wrapped to (0.32, 0.2, 0.2)
    // B_wrapped_cart = (2.0, 1.0, 1.0)
    // A stays at (1,1,1) after wrapping.
    // offset [0,0,0]: diff=(1,0,0), dist=1.0 < 1.5 ✓
    // Naive Cartesian wrap would give B=(4,1,1), dist=3.0 → MISS
    let sys = System {
        pos: vec![Vector3::new(1.0, 1.0, 1.0), Vector3::new(14.0, 6.0, 1.0)],
        cell: UnitCell::new(
            Vector3::new(5.0, 0.0, 0.0),
            Vector3::new(2.0, 5.0, 0.0),
            Vector3::new(0.0, 0.0, 5.0),
        )
        .unwrap(),
    };
    let result = build(&sys, 1.5);
    assert_eq!(result.len(), 2);
    let a_to_b: Vec<&Neighbor> = result.iter().filter(|n| n.i == 0 && n.j == 1).collect();
    assert_eq!(a_to_b.len(), 1);
    assert_eq!(a_to_b[0].offset, [0, 0, 0]);
    assert!((a_to_b[0].distance - 1.0).abs() < 1e-10);
}

pub fn test_triclinic_self_image(build: BuildFn) {
    // Triclinic cell: a=(5,0,0), b=(2,5,0), c=(0,0,5)
    // Single atom at origin, cutoff=5.5
    // Self-images:
    //   [±1,0,0]: dist=|a|=5.0
    //   [0,±1,0]: dist=|b|=sqrt(29)≈5.39
    //   [0,0,±1]: dist=|c|=5.0
    // All 6 within cutoff. Diagonals like [±1,∓1,0] have dist=sqrt(34)≈5.83 > 5.5.
    let sys = System {
        pos: vec![Vector3::new(0.0, 0.0, 0.0)],
        cell: UnitCell::new(
            Vector3::new(5.0, 0.0, 0.0),
            Vector3::new(2.0, 5.0, 0.0),
            Vector3::new(0.0, 0.0, 5.0),
        )
        .unwrap(),
    };
    let result = build(&sys, 5.5);
    assert!(result.iter().all(|n| n.i == 0 && n.j == 0));
    assert_eq!(result.len(), 6);
    let offsets: Vec<[i32; 3]> = result.iter().map(|n| n.offset).collect();
    assert!(offsets.contains(&[1, 0, 0]));
    assert!(offsets.contains(&[-1, 0, 0]));
    assert!(offsets.contains(&[0, 1, 0]));
    assert!(offsets.contains(&[0, -1, 0]));
    assert!(offsets.contains(&[0, 0, 1]));
    assert!(offsets.contains(&[0, 0, -1]));
    // ±a and ±c have dist 5.0, ±b have dist sqrt(29)
    for n in &result {
        match n.offset {
            [1, 0, 0] | [-1, 0, 0] | [0, 0, 1] | [0, 0, -1] => {
                assert!((n.distance - 5.0).abs() < 1e-10);
            }
            [0, 1, 0] | [0, -1, 0] => {
                assert!((n.distance - 29.0_f64.sqrt()).abs() < 1e-10);
            }
            _ => panic!("unexpected offset {:?}", n.offset),
        }
    }
}

pub fn test_triclinic_skewed(build: BuildFn) {
    // Highly skewed cell: a=(10,0,0), b=(9.5,1,0), c=(0,0,5)
    // d_b = V/|a×c| = 50/50 = 1.0 (very thin slab in b-direction)
    // |b| = sqrt(91.25) ≈ 9.55
    // Height-based replicas: ceil(2.5/1.0) = 3 for b
    // Norm-based would give: ceil(2.5/9.55) = 1 → misses [±2,∓2,0]
    //
    // Single atom at origin, cutoff=2.5
    // Self-images:
    //   [1,-1,0]: |a-b| = |(0.5,-1,0)| = sqrt(1.25) ✓
    //   [-1,1,0]: |b-a| = sqrt(1.25) ✓
    //   [2,-2,0]: |2a-2b| = |(1,-2,0)| = sqrt(5) ✓
    //   [-2,2,0]: |2b-2a| = sqrt(5) ✓
    // All other offsets give dist > 2.5
    let sys = System {
        pos: vec![Vector3::new(0.0, 0.0, 0.0)],
        cell: UnitCell::new(
            Vector3::new(10.0, 0.0, 0.0),
            Vector3::new(9.5, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 5.0),
        )
        .unwrap(),
    };
    let result = build(&sys, 2.5);
    assert!(result.iter().all(|n| n.i == 0 && n.j == 0));
    assert_eq!(result.len(), 4);
    let offsets: Vec<[i32; 3]> = result.iter().map(|n| n.offset).collect();
    assert!(offsets.contains(&[1, -1, 0]));
    assert!(offsets.contains(&[-1, 1, 0]));
    assert!(offsets.contains(&[2, -2, 0]));
    assert!(offsets.contains(&[-2, 2, 0]));
    for n in &result {
        match n.offset {
            [1, -1, 0] | [-1, 1, 0] => {
                assert!((n.distance - 1.25_f64.sqrt()).abs() < 1e-10);
            }
            [2, -2, 0] | [-2, 2, 0] => {
                assert!((n.distance - 5.0_f64.sqrt()).abs() < 1e-10);
            }
            _ => panic!("unexpected offset {:?}", n.offset),
        }
    }
}
