use crate::types::{Neighbor, System, UnitCell, Vector3};

pub type BuildFn = fn(&System, f64, bool) -> Vec<Neighbor>;

/// Helper to sort neighbors for deterministic comparison.
fn sorted_neighbors(neighbors: &[Neighbor]) -> Vec<(usize, usize, [i32; 3])> {
    let mut v: Vec<(usize, usize, [i32; 3])> =
        neighbors.iter().map(|n| (n.i, n.j, n.offset)).collect();
    v.sort();
    v
}

pub fn test_basic(build: BuildFn) {
    // [Cubic] 4 atoms in a 10x10x10 box, cutoff = 3.0
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
        pbc: [true, true, true],
    };
    let result = build(&sys, 3.0, false);
    // Both directions: (0,1), (1,0), (0,3), (3,0)
    let pairs: Vec<(usize, usize)> = result.iter().map(|n| (n.i, n.j)).collect();
    assert!(pairs.contains(&(0, 1)));
    assert!(pairs.contains(&(1, 0)));
    assert!(pairs.contains(&(0, 3)));
    assert!(pairs.contains(&(3, 0)));
    assert_eq!(result.len(), 4);
}

pub fn test_no_neighbors(build: BuildFn) {
    // [Cubic] 20x20x20 box, cutoff=1.0 — no pairs within range
    let sys = System {
        pos: vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(5.0, 5.0, 5.0)],
        cell: UnitCell::new(
            Vector3::new(20.0, 0.0, 0.0),
            Vector3::new(0.0, 20.0, 0.0),
            Vector3::new(0.0, 0.0, 20.0),
        )
        .unwrap(),
        pbc: [true, true, true],
    };
    let result = build(&sys, 1.0, false);
    assert!(result.is_empty());
}

pub fn test_pbc_corner(build: BuildFn) {
    // [Cubic] A(0.5, 0.5, 0.5) B(9.5, 9.5, 9.5) in 10x10x10 box
    // dist via PBC = sqrt(3) ≈ 1.73
    let sys = System {
        pos: vec![Vector3::new(0.5, 0.5, 0.5), Vector3::new(9.5, 9.5, 9.5)],
        cell: UnitCell::new(
            Vector3::new(10.0, 0.0, 0.0),
            Vector3::new(0.0, 10.0, 0.0),
            Vector3::new(0.0, 0.0, 10.0),
        )
        .unwrap(),
        pbc: [true, true, true],
    };
    let result = build(&sys, 2.0, false);
    // Both directions
    let pairs: Vec<(usize, usize)> = result.iter().map(|n| (n.i, n.j)).collect();
    assert_eq!(pairs.len(), 2);
    assert!(pairs.contains(&(0, 1)));
    assert!(pairs.contains(&(1, 0)));
}

pub fn test_self_image(build: BuildFn) {
    // [Cubic] Single atom at origin, L=3, cutoff=3.0 (== L, boundary case)
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
        pbc: [true, true, true],
    };
    let result = build(&sys, 3.0, false);
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
    // [Cubic] 2 atoms in L=4 box, cutoff=5
    // A(0,0,0) B(1,0,0)
    let sys = System {
        pos: vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 0.0, 0.0)],
        cell: UnitCell::new(
            Vector3::new(4.0, 0.0, 0.0),
            Vector3::new(0.0, 4.0, 0.0),
            Vector3::new(0.0, 0.0, 4.0),
        )
        .unwrap(),
        pbc: [true, true, true],
    };
    let result = build(&sys, 5.0, false);

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

pub fn test_orthorhombic(build: BuildFn) {
    // [Orthorhombic] a=3, b=10, c=10, cutoff=3.5
    // A(0,0,0) B(1,0,0)
    let sys = System {
        pos: vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 0.0, 0.0)],
        cell: UnitCell::new(
            Vector3::new(3.0, 0.0, 0.0),
            Vector3::new(0.0, 10.0, 0.0),
            Vector3::new(0.0, 0.0, 10.0),
        )
        .unwrap(),
        pbc: [true, true, true],
    };
    let result = build(&sys, 3.5, false);

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
    // [Cubic] Atom j is far outside the primary cell (diff > L).
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
        pbc: [true, true, true],
    };
    let result = build(&sys, 3.0, false);
    assert_eq!(result.len(), 2); // i→j and j→i
    let i_to_j: Vec<&Neighbor> = result.iter().filter(|n| n.i == 0 && n.j == 1).collect();
    assert_eq!(i_to_j.len(), 1);
    assert_eq!(i_to_j[0].offset, [-1, 0, 0]);
    assert!((i_to_j[0].distance - 6.0_f64.sqrt()).abs() < 1e-10);
}

pub fn test_monoclinic_basic(build: BuildFn) {
    // [Monoclinic] a=(5,0,0), b=(2,5,0), c=(0,0,5) — α=β=90°, γ≠90°
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
        pbc: [true, true, true],
    };
    let result = build(&sys, 2.9, false);
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

pub fn test_monoclinic_pbc_corner(build: BuildFn) {
    // Monoclinic: a=(5,0,0), b=(2,5,0), c=(0,0,5) — α=β=90°, γ≠90°
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
        pbc: [true, true, true],
    };
    let result = build(&sys, 1.5, false);
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

pub fn test_monoclinic_unwrapped(build: BuildFn) {
    // Monoclinic: a=(5,0,0), b=(2,5,0), c=(0,0,5) — α=β=90°, γ≠90°
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
        pbc: [true, true, true],
    };
    let result = build(&sys, 1.5, false);
    assert_eq!(result.len(), 2);
    let a_to_b: Vec<&Neighbor> = result.iter().filter(|n| n.i == 0 && n.j == 1).collect();
    assert_eq!(a_to_b.len(), 1);
    assert_eq!(a_to_b[0].offset, [0, 0, 0]);
    assert!((a_to_b[0].distance - 1.0).abs() < 1e-10);
}

pub fn test_monoclinic_self_image(build: BuildFn) {
    // Monoclinic: a=(5,0,0), b=(2,5,0), c=(0,0,5) — α=β=90°, γ≠90°
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
        pbc: [true, true, true],
    };
    let result = build(&sys, 5.5, false);
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

pub fn test_monoclinic_skewed(build: BuildFn) {
    // Monoclinic (highly skewed): a=(10,0,0), b=(9.5,1,0), c=(0,0,5) — α=β=90°, γ≠90°
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
        pbc: [true, true, true],
    };
    let result = build(&sys, 2.5, false);
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

pub fn test_triclinic_full_3d(build: BuildFn) {
    // Full 3D triclinic: a=(6,0,0), b=(1,5,0), c=(0.5,1,4)
    // α≈75°, β≈83°, γ≈79° — all angles ≠ 90°
    // A=frac(0.1,0.1,0.1), B=frac(0.9,0.9,0.9), cutoff=3.5
    // A_cart = (0.75, 0.6, 0.4), B_cart = (6.75, 5.4, 3.6)
    //
    // 0→1 offset [-1,-1,-1]: diff=(-1.5,-1.2,-0.8), dist=sqrt(4.33)≈2.08 ✓
    // 0→1 offset [-1,-1, 0]: diff=(-1.0,-0.2, 3.2), dist=sqrt(11.28)≈3.36 ✓
    // 1→0: reverse offsets [1,1,1] and [1,1,0], same distances
    let cell = UnitCell::new(
        Vector3::new(6.0, 0.0, 0.0),
        Vector3::new(1.0, 5.0, 0.0),
        Vector3::new(0.5, 1.0, 4.0),
    )
    .unwrap();
    let a_cart = cell.get_cartesian(&Vector3::new(0.1, 0.1, 0.1));
    let b_cart = cell.get_cartesian(&Vector3::new(0.9, 0.9, 0.9));
    let sys = System {
        pos: vec![a_cart, b_cart],
        cell,
        pbc: [true, true, true],
    };
    let result = build(&sys, 3.5, false);
    assert_eq!(result.len(), 4);

    let a_to_b: Vec<&Neighbor> = result.iter().filter(|n| n.i == 0 && n.j == 1).collect();
    assert_eq!(a_to_b.len(), 2);
    let offsets_ab: Vec<[i32; 3]> = a_to_b.iter().map(|n| n.offset).collect();
    assert!(offsets_ab.contains(&[-1, -1, -1]));
    assert!(offsets_ab.contains(&[-1, -1, 0]));
    for n in &a_to_b {
        match n.offset {
            [-1, -1, -1] => assert!((n.distance - (4.33_f64).sqrt()).abs() < 1e-10),
            [-1, -1, 0] => assert!((n.distance - (11.28_f64).sqrt()).abs() < 1e-10),
            _ => panic!("unexpected offset {:?}", n.offset),
        }
    }

    let b_to_a: Vec<&Neighbor> = result.iter().filter(|n| n.i == 1 && n.j == 0).collect();
    assert_eq!(b_to_a.len(), 2);
    let offsets_ba: Vec<[i32; 3]> = b_to_a.iter().map(|n| n.offset).collect();
    assert!(offsets_ba.contains(&[1, 1, 1]));
    assert!(offsets_ba.contains(&[1, 1, 0]));
}

pub fn test_triclinic_full_3d_self_image(build: BuildFn) {
    // Full 3D triclinic: a=(4,0,0), b=(1,4,0), c=(1,1,4)
    // α≈73°, β≈76°, γ≈76°
    // Single atom at origin, cutoff=4.3
    // Self-images:
    //   [±1,0,0]: |a|=4.0 ✓
    //   [0,±1,0]: |b|=sqrt(17)≈4.12 ✓
    //   [0,0,±1]: |c|=sqrt(18)≈4.24 ✓
    // Diagonals: [1,-1,0]=|(3,-4,0)|=5.0 > 4.3, etc.
    let sys = System {
        pos: vec![Vector3::new(0.0, 0.0, 0.0)],
        cell: UnitCell::new(
            Vector3::new(4.0, 0.0, 0.0),
            Vector3::new(1.0, 4.0, 0.0),
            Vector3::new(1.0, 1.0, 4.0),
        )
        .unwrap(),
        pbc: [true, true, true],
    };
    let result = build(&sys, 4.3, false);
    assert!(result.iter().all(|n| n.i == 0 && n.j == 0));
    assert_eq!(result.len(), 6);
    let offsets: Vec<[i32; 3]> = result.iter().map(|n| n.offset).collect();
    assert!(offsets.contains(&[1, 0, 0]));
    assert!(offsets.contains(&[-1, 0, 0]));
    assert!(offsets.contains(&[0, 1, 0]));
    assert!(offsets.contains(&[0, -1, 0]));
    assert!(offsets.contains(&[0, 0, 1]));
    assert!(offsets.contains(&[0, 0, -1]));
    for n in &result {
        match n.offset {
            [1, 0, 0] | [-1, 0, 0] => assert!((n.distance - 4.0).abs() < 1e-10),
            [0, 1, 0] | [0, -1, 0] => assert!((n.distance - 17.0_f64.sqrt()).abs() < 1e-10),
            [0, 0, 1] | [0, 0, -1] => assert!((n.distance - 18.0_f64.sqrt()).abs() < 1e-10),
            _ => panic!("unexpected offset {:?}", n.offset),
        }
    }
}

pub fn test_monoclinic_boundary_distance(build: BuildFn) {
    // Monoclinic: a=(5,0,0), b=(2,5,0), c=(0,0,5) — α=β=90°, γ≠90°
    // A=(0,0,0), B=(2,0,0), cutoff=3.0
    //
    // 0→1: [0,0,0] dist=2.0 ✓, [-1,0,0] dist=3.0 == cutoff ✓
    // 1→0: [0,0,0] dist=2.0 ✓, [1,0,0] dist=3.0 == cutoff ✓
    let sys = System {
        pos: vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(2.0, 0.0, 0.0)],
        cell: UnitCell::new(
            Vector3::new(5.0, 0.0, 0.0),
            Vector3::new(2.0, 5.0, 0.0),
            Vector3::new(0.0, 0.0, 5.0),
        )
        .unwrap(),
        pbc: [true, true, true],
    };
    let result = build(&sys, 3.0, false);
    assert_eq!(result.len(), 4);

    let a_to_b: Vec<&Neighbor> = result.iter().filter(|n| n.i == 0 && n.j == 1).collect();
    assert_eq!(a_to_b.len(), 2);
    assert!(
        a_to_b
            .iter()
            .any(|n| n.offset == [0, 0, 0] && (n.distance - 2.0).abs() < 1e-10)
    );
    assert!(
        a_to_b
            .iter()
            .any(|n| n.offset == [-1, 0, 0] && (n.distance - 3.0).abs() < 1e-10)
    );

    let b_to_a: Vec<&Neighbor> = result.iter().filter(|n| n.i == 1 && n.j == 0).collect();
    assert_eq!(b_to_a.len(), 2);
    assert!(
        b_to_a
            .iter()
            .any(|n| n.offset == [0, 0, 0] && (n.distance - 2.0).abs() < 1e-10)
    );
    assert!(
        b_to_a
            .iter()
            .any(|n| n.offset == [1, 0, 0] && (n.distance - 3.0).abs() < 1e-10)
    );
}

pub fn test_monoclinic_left_handed(build: BuildFn) {
    // Monoclinic (left-handed, det < 0): a=(5,0,0), b=(0,0,5), c=(2,5,0) — α=γ=90°, β≠90°
    // det = -125
    // Atoms: frac(0.05,0.05,0.05), frac(0.95,0.95,0.95), cutoff=1.5
    // A_cart = (0.35, 0.25, 0.25), B_cart = (6.65, 4.75, 4.75)
    // 0→1 offset [-1,-1,-1]: diff=(-0.7,-0.5,-0.5), dist=sqrt(0.99)≈0.995 ✓
    let cell = UnitCell::new(
        Vector3::new(5.0, 0.0, 0.0),
        Vector3::new(0.0, 0.0, 5.0),
        Vector3::new(2.0, 5.0, 0.0),
    )
    .unwrap();
    let a_cart = cell.get_cartesian(&Vector3::new(0.05, 0.05, 0.05));
    let b_cart = cell.get_cartesian(&Vector3::new(0.95, 0.95, 0.95));
    let sys = System {
        pos: vec![a_cart, b_cart],
        cell,
        pbc: [true, true, true],
    };
    let result = build(&sys, 1.5, false);
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

pub fn test_triclinic_multi_height_replicas(build: BuildFn) {
    // Highly skewed 3D cell: a=(10,0,0), b=(9.5,1,0), c=(9.5,0,1)
    // V=10, d_a≈0.74, d_b=1.0, d_c=1.0
    // cutoff=1.5, replicas=(3,2,2) — all 3 directions need >1 replica
    //
    // Single atom at origin, self-images within 1.5:
    //   [1,-1,0]/[-1,1,0]: |(0.5,-1,0)| = sqrt(1.25) ✓
    //   [1,0,-1]/[-1,0,1]: |(0.5,0,-1)| = sqrt(1.25) ✓
    //   [0,1,-1]/[0,-1,1]: |(0,1,-1)|   = sqrt(2.0)   ✓
    let sys = System {
        pos: vec![Vector3::new(0.0, 0.0, 0.0)],
        cell: UnitCell::new(
            Vector3::new(10.0, 0.0, 0.0),
            Vector3::new(9.5, 1.0, 0.0),
            Vector3::new(9.5, 0.0, 1.0),
        )
        .unwrap(),
        pbc: [true, true, true],
    };
    let result = build(&sys, 1.5, false);
    assert!(result.iter().all(|n| n.i == 0 && n.j == 0));
    assert_eq!(result.len(), 6);
    let offsets: Vec<[i32; 3]> = result.iter().map(|n| n.offset).collect();
    assert!(offsets.contains(&[1, -1, 0]));
    assert!(offsets.contains(&[-1, 1, 0]));
    assert!(offsets.contains(&[1, 0, -1]));
    assert!(offsets.contains(&[-1, 0, 1]));
    assert!(offsets.contains(&[0, 1, -1]));
    assert!(offsets.contains(&[0, -1, 1]));
    for n in &result {
        match n.offset {
            [1, -1, 0] | [-1, 1, 0] | [1, 0, -1] | [-1, 0, 1] => {
                assert!((n.distance - 1.25_f64.sqrt()).abs() < 1e-10);
            }
            [0, 1, -1] | [0, -1, 1] => {
                assert!((n.distance - 2.0_f64.sqrt()).abs() < 1e-10);
            }
            _ => panic!("unexpected offset {:?}", n.offset),
        }
    }
}

pub fn test_triclinic_multi_atom(build: BuildFn) {
    // Full 3D triclinic with 3 atoms + boundary distance
    // Cell: a=(4,0,0), b=(1,4,0), c=(1,1,4)
    // P0=(0,0,0)=frac(0,0,0), P1=(2,0,0)=frac(0.5,0,0), P2=(0.5,2,0)=frac(0,0.5,0)
    // cutoff=2.5
    //
    // 0→1: [0,0,0] dist=2.0, [-1,0,0] dist=2.0 → 2 neighbors
    // 0→2: [0,0,0] dist=sqrt(4.25), [0,-1,0] dist=sqrt(4.25) → 2 neighbors
    // 1→2: [0,0,0] dist=2.5, [1,-1,0] dist=2.5 (both == cutoff) → 2 neighbors
    // Reverse pairs: 1→0: 2, 2→0: 2, 2→1: 2
    // No self-images (|a|=4, |b|=sqrt(17), |c|=sqrt(18) all > 2.5)
    // Total: 12 neighbors
    let sys = System {
        pos: vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(2.0, 0.0, 0.0),
            Vector3::new(0.5, 2.0, 0.0),
        ],
        cell: UnitCell::new(
            Vector3::new(4.0, 0.0, 0.0),
            Vector3::new(1.0, 4.0, 0.0),
            Vector3::new(1.0, 1.0, 4.0),
        )
        .unwrap(),
        pbc: [true, true, true],
    };
    let result = build(&sys, 2.5, false);
    assert_eq!(result.len(), 12);

    // 0→1
    let p01: Vec<&Neighbor> = result.iter().filter(|n| n.i == 0 && n.j == 1).collect();
    assert_eq!(p01.len(), 2);
    assert!(
        p01.iter()
            .any(|n| n.offset == [0, 0, 0] && (n.distance - 2.0).abs() < 1e-10)
    );
    assert!(
        p01.iter()
            .any(|n| n.offset == [-1, 0, 0] && (n.distance - 2.0).abs() < 1e-10)
    );

    // 0→2
    let p02: Vec<&Neighbor> = result.iter().filter(|n| n.i == 0 && n.j == 2).collect();
    assert_eq!(p02.len(), 2);
    assert!(
        p02.iter()
            .any(|n| n.offset == [0, 0, 0] && (n.distance - 4.25_f64.sqrt()).abs() < 1e-10)
    );
    assert!(
        p02.iter()
            .any(|n| n.offset == [0, -1, 0] && (n.distance - 4.25_f64.sqrt()).abs() < 1e-10)
    );

    // 1→2: both at boundary distance == cutoff
    let p12: Vec<&Neighbor> = result.iter().filter(|n| n.i == 1 && n.j == 2).collect();
    assert_eq!(p12.len(), 2);
    assert!(
        p12.iter()
            .any(|n| n.offset == [0, 0, 0] && (n.distance - 2.5).abs() < 1e-10)
    );
    assert!(
        p12.iter()
            .any(|n| n.offset == [1, -1, 0] && (n.distance - 2.5).abs() < 1e-10)
    );
}

pub fn test_triclinic_near_degenerate(build: BuildFn) {
    // Near-degenerate cell: a=(5,0,0), b=(4.9,0.5,0), c=(4.9,0,0.5)
    // V=1.25, d_a≈0.36, d_b=0.5, d_c=0.5
    // cutoff=1.1, replicas=(4,3,3)
    // Norm-based replicas(1,1,1) would miss offsets with |na|=2
    //
    // Single atom at origin, 12 self-images:
    //   [±1,∓1,0], [±1,0,∓1]: dist=sqrt(0.26) ≈ 0.51
    //   [0,±1,∓1]:             dist=sqrt(0.50) ≈ 0.71
    //   [±2,∓1,∓1]:            dist=sqrt(0.54) ≈ 0.73
    //   [±2,∓2,0], [±2,0,∓2]: dist=sqrt(1.04) ≈ 1.02
    let sys = System {
        pos: vec![Vector3::new(0.0, 0.0, 0.0)],
        cell: UnitCell::new(
            Vector3::new(5.0, 0.0, 0.0),
            Vector3::new(4.9, 0.5, 0.0),
            Vector3::new(4.9, 0.0, 0.5),
        )
        .unwrap(),
        pbc: [true, true, true],
    };
    let result = build(&sys, 1.1, false);
    assert!(result.iter().all(|n| n.i == 0 && n.j == 0));
    assert_eq!(result.len(), 12);
    let offsets: Vec<[i32; 3]> = result.iter().map(|n| n.offset).collect();
    // sqrt(0.26) group
    assert!(offsets.contains(&[1, -1, 0]));
    assert!(offsets.contains(&[-1, 1, 0]));
    assert!(offsets.contains(&[1, 0, -1]));
    assert!(offsets.contains(&[-1, 0, 1]));
    // sqrt(0.50) group
    assert!(offsets.contains(&[0, 1, -1]));
    assert!(offsets.contains(&[0, -1, 1]));
    // sqrt(0.54) group
    assert!(offsets.contains(&[2, -1, -1]));
    assert!(offsets.contains(&[-2, 1, 1]));
    // sqrt(1.04) group
    assert!(offsets.contains(&[2, -2, 0]));
    assert!(offsets.contains(&[-2, 2, 0]));
    assert!(offsets.contains(&[2, 0, -2]));
    assert!(offsets.contains(&[-2, 0, 2]));

    for n in &result {
        match n.offset {
            [1, -1, 0] | [-1, 1, 0] | [1, 0, -1] | [-1, 0, 1] => {
                assert!((n.distance - 0.26_f64.sqrt()).abs() < 1e-10);
            }
            [0, 1, -1] | [0, -1, 1] => {
                assert!((n.distance - 0.50_f64.sqrt()).abs() < 1e-10);
            }
            [2, -1, -1] | [-2, 1, 1] => {
                assert!((n.distance - 0.54_f64.sqrt()).abs() < 1e-10);
            }
            [2, -2, 0] | [-2, 2, 0] | [2, 0, -2] | [-2, 0, 2] => {
                assert!((n.distance - 1.04_f64.sqrt()).abs() < 1e-10);
            }
            _ => panic!("unexpected offset {:?}", n.offset),
        }
    }
}

pub fn test_slab_self_image(build: BuildFn) {
    // [Cubic] pbc=[true,true,false] — slab: periodic in a,b only
    // Single atom at origin, L=3, cutoff=3.0
    // Self-images along a (±1,0,0) and b (0,±1,0) → 4 neighbors
    // No c-axis images because pbc[2]=false
    let sys = System {
        pos: vec![Vector3::new(0.0, 0.0, 0.0)],
        cell: UnitCell::new(
            Vector3::new(3.0, 0.0, 0.0),
            Vector3::new(0.0, 3.0, 0.0),
            Vector3::new(0.0, 0.0, 3.0),
        )
        .unwrap(),
        pbc: [true, true, false],
    };
    let result = build(&sys, 3.0, false);
    assert!(result.iter().all(|n| n.i == 0 && n.j == 0));
    assert_eq!(result.len(), 4);
    assert!(result.iter().all(|n| (n.distance - 3.0).abs() < 1e-10));
    let offsets: Vec<[i32; 3]> = result.iter().map(|n| n.offset).collect();
    assert!(offsets.contains(&[1, 0, 0]));
    assert!(offsets.contains(&[-1, 0, 0]));
    assert!(offsets.contains(&[0, 1, 0]));
    assert!(offsets.contains(&[0, -1, 0]));
    // c-axis images must NOT appear
    assert!(!offsets.contains(&[0, 0, 1]));
    assert!(!offsets.contains(&[0, 0, -1]));
}

pub fn test_wire_self_image(build: BuildFn) {
    // [Cubic] pbc=[true,false,false] — wire: periodic in a only
    // Single atom at origin, L=3, cutoff=3.0
    // Self-images along a only → 2 neighbors
    let sys = System {
        pos: vec![Vector3::new(0.0, 0.0, 0.0)],
        cell: UnitCell::new(
            Vector3::new(3.0, 0.0, 0.0),
            Vector3::new(0.0, 3.0, 0.0),
            Vector3::new(0.0, 0.0, 3.0),
        )
        .unwrap(),
        pbc: [true, false, false],
    };
    let result = build(&sys, 3.0, false);
    assert!(result.iter().all(|n| n.i == 0 && n.j == 0));
    assert_eq!(result.len(), 2);
    assert!(result.iter().all(|n| (n.distance - 3.0).abs() < 1e-10));
    let offsets: Vec<[i32; 3]> = result.iter().map(|n| n.offset).collect();
    assert!(offsets.contains(&[1, 0, 0]));
    assert!(offsets.contains(&[-1, 0, 0]));
}

pub fn test_isolated_no_images(build: BuildFn) {
    // [Cubic] pbc=[false,false,false] — isolated: no periodicity
    // Single atom at origin, L=3, cutoff=3.0
    // No self-images at all → 0 neighbors
    let sys = System {
        pos: vec![Vector3::new(0.0, 0.0, 0.0)],
        cell: UnitCell::new(
            Vector3::new(3.0, 0.0, 0.0),
            Vector3::new(0.0, 3.0, 0.0),
            Vector3::new(0.0, 0.0, 3.0),
        )
        .unwrap(),
        pbc: [false, false, false],
    };
    let result = build(&sys, 3.0, false);
    assert!(result.is_empty());
}

pub fn test_slab_two_atoms(build: BuildFn) {
    // [Cubic] pbc=[true,true,false] — slab, 10x10x10, cutoff=2.0
    // A(0.5,0.5,0.5), B(9.5,0.5,0.5)
    // a-axis wraps: diff_x = -1.0, dist = 1.0 ✓
    // C(0.5,0.5,9.5) — c-axis does NOT wrap: diff_z = 9.0, dist = 9.0 → NOT neighbor
    let sys = System {
        pos: vec![
            Vector3::new(0.5, 0.5, 0.5),
            Vector3::new(9.5, 0.5, 0.5),
            Vector3::new(0.5, 0.5, 9.5),
        ],
        cell: UnitCell::new(
            Vector3::new(10.0, 0.0, 0.0),
            Vector3::new(0.0, 10.0, 0.0),
            Vector3::new(0.0, 0.0, 10.0),
        )
        .unwrap(),
        pbc: [true, true, false],
    };
    let result = build(&sys, 2.0, false);
    // A↔B via a-axis PBC: 2 pairs
    assert_eq!(result.len(), 2);
    let pairs: Vec<(usize, usize)> = result.iter().map(|n| (n.i, n.j)).collect();
    assert!(pairs.contains(&(0, 1)));
    assert!(pairs.contains(&(1, 0)));
    // A-C should NOT be neighbors (c-axis non-periodic, dist=9.0)
    assert!(!pairs.contains(&(0, 2)));
    assert!(!pairs.contains(&(2, 0)));
}

pub fn test_isolated_two_atoms(build: BuildFn) {
    // [Cubic] pbc=[false,false,false] — isolated, L=10, cutoff=3.0
    // A(1,1,1), B(3,1,1): direct dist=2.0 ✓ (no wrapping needed)
    // C(9,1,1): direct dist=8.0 → NOT neighbor (no PBC wrapping)
    //   (with full PBC, dist via wrapping would be 2.0)
    let sys = System {
        pos: vec![
            Vector3::new(1.0, 1.0, 1.0),
            Vector3::new(3.0, 1.0, 1.0),
            Vector3::new(9.0, 1.0, 1.0),
        ],
        cell: UnitCell::new(
            Vector3::new(10.0, 0.0, 0.0),
            Vector3::new(0.0, 10.0, 0.0),
            Vector3::new(0.0, 0.0, 10.0),
        )
        .unwrap(),
        pbc: [false, false, false],
    };
    let result = build(&sys, 3.0, false);
    // Only A↔B (direct distance 2.0)
    assert_eq!(result.len(), 2);
    let a_to_b: Vec<&Neighbor> = result.iter().filter(|n| n.i == 0 && n.j == 1).collect();
    assert_eq!(a_to_b.len(), 1);
    assert_eq!(a_to_b[0].offset, [0, 0, 0]);
    assert!((a_to_b[0].distance - 2.0).abs() < 1e-10);
    // C should NOT appear as neighbor of A
    assert!(result.iter().all(|n| n.i != 2 || n.j != 0));
    assert!(result.iter().all(|n| n.i != 0 || n.j != 2));
}

pub fn test_slab_no_wrap_non_periodic(build: BuildFn) {
    // [Cubic] pbc=[true,true,false], L=5, cutoff=2.0
    // Atom at frac=(0.1, 0.1, 1.5) — outside cell along c-axis
    // With pbc[2]=false, c-coordinate should NOT be wrapped
    // So Cartesian z = 7.5 (stays at 7.5, not wrapped to 2.5)
    // Another atom at (0.5, 0.5, 7.0)
    // Direct dist = sqrt(0.16+0.16+0.25) = sqrt(0.57) ≈ 0.75 ✓
    // If c were wrapped: atom0 would be at z=2.5, dist to atom1 = sqrt(0.16+0.16+20.25) → miss
    let cell = UnitCell::new(
        Vector3::new(5.0, 0.0, 0.0),
        Vector3::new(0.0, 5.0, 0.0),
        Vector3::new(0.0, 0.0, 5.0),
    )
    .unwrap();
    let sys = System {
        pos: vec![
            cell.get_cartesian(&Vector3::new(0.1, 0.1, 1.5)),
            Vector3::new(0.5, 0.5, 7.0),
        ],
        cell,
        pbc: [true, true, false],
    };
    let result = build(&sys, 2.0, false);
    assert_eq!(result.len(), 2);
    let a_to_b: Vec<&Neighbor> = result.iter().filter(|n| n.i == 0 && n.j == 1).collect();
    assert_eq!(a_to_b.len(), 1);
    assert_eq!(a_to_b[0].offset, [0, 0, 0]);
}

pub fn test_slab_cross_bin_boundary(build: BuildFn) {
    // [Cubic] pbc=[true,true,false], L=10, cutoff=3.0
    // n_bins along z = floor(10/3) = 3, bin boundaries at frac 0, 1/3, 2/3, 1
    //   → Cartesian z = 0, 3.333, 6.667, 10
    // A(5, 5, 3.3) → frac_z=0.33 → bin 0
    // B(5, 5, 3.4) → frac_z=0.34 → bin 1
    // distance = 0.1, well within cutoff
    // Bug: if search_ranges[2]=0, bin 0 won't look at bin 1 → pair missed
    let sys = System {
        pos: vec![Vector3::new(5.0, 5.0, 3.3), Vector3::new(5.0, 5.0, 3.4)],
        cell: UnitCell::new(
            Vector3::new(10.0, 0.0, 0.0),
            Vector3::new(0.0, 10.0, 0.0),
            Vector3::new(0.0, 0.0, 10.0),
        )
        .unwrap(),
        pbc: [true, true, false],
    };
    let result = build(&sys, 3.0, false);
    assert_eq!(result.len(), 2);
    let a_to_b: Vec<&Neighbor> = result.iter().filter(|n| n.i == 0 && n.j == 1).collect();
    assert_eq!(a_to_b.len(), 1);
    assert_eq!(a_to_b[0].offset, [0, 0, 0]);
    assert!((a_to_b[0].distance - 0.1).abs() < 1e-10);
}

pub fn test_isolated_atoms_outside_cell(build: BuildFn) {
    // [Cubic] pbc=[false,false,false], L=10, cutoff=2.0
    // A at frac(-0.5, 0.5, 0.5) → Cartesian(-5, 5, 5) — outside cell
    // B at frac(-0.3, 0.5, 0.5) → Cartesian(-3, 5, 5) — outside cell
    // Direct distance = 2.0, within cutoff
    // Both atoms are outside the cell on the negative side.
    // Without bin index clamping, this would panic or produce wrong results.
    let sys = System {
        pos: vec![Vector3::new(-5.0, 5.0, 5.0), Vector3::new(-3.0, 5.0, 5.0)],
        cell: UnitCell::new(
            Vector3::new(10.0, 0.0, 0.0),
            Vector3::new(0.0, 10.0, 0.0),
            Vector3::new(0.0, 0.0, 10.0),
        )
        .unwrap(),
        pbc: [false, false, false],
    };
    let result = build(&sys, 2.0, false);
    assert_eq!(result.len(), 2);
    let a_to_b: Vec<&Neighbor> = result.iter().filter(|n| n.i == 0 && n.j == 1).collect();
    assert_eq!(a_to_b.len(), 1);
    assert_eq!(a_to_b[0].offset, [0, 0, 0]);
    assert!((a_to_b[0].distance - 2.0).abs() < 1e-10);
}

/// Shared check: given a full and half list for the same system, verify that
/// the half list is exactly the canonical half of the full list.
fn verify_half_invariants(full: &[Neighbor], half: &[Neighbor]) {
    // Every half entry must be canonical.
    for n in half {
        assert!(
            n.is_half_canonical(),
            "non-canonical entry in half list: (i={}, j={}, offset={:?})",
            n.i,
            n.j,
            n.offset,
        );
    }

    // In a full list, every physical pair appears twice, so the length must be even
    // and the half list must be exactly half of it.
    assert_eq!(full.len() % 2, 0, "full list length ({}) is odd", full.len());
    assert_eq!(
        half.len() * 2,
        full.len(),
        "half list length ({}) is not half of full ({})",
        half.len(),
        full.len(),
    );

    // The canonical entries extracted from the full list must match the half list exactly.
    let mut expected: Vec<(usize, usize, [i32; 3])> = full
        .iter()
        .filter(|n| n.is_half_canonical())
        .map(|n| (n.i, n.j, n.offset))
        .collect();
    expected.sort();
    assert_eq!(
        sorted_neighbors(half),
        expected,
        "half list does not match canonical entries of full list",
    );

    // Distances must match between full and half for each canonical pair.
    for hn in half {
        let fne = full
            .iter()
            .find(|fne| fne.i == hn.i && fne.j == hn.j && fne.offset == hn.offset)
            .expect("half pair not found in full list");
        assert!(
            (fne.distance - hn.distance).abs() < 1e-12,
            "distance mismatch for (i={}, j={}, offset={:?}): full={}, half={}",
            hn.i,
            hn.j,
            hn.offset,
            fne.distance,
            hn.distance,
        );
    }
}

/// Comprehensive half-neighbor-list test suite. Runs `build` in both full and half
/// mode over several representative systems and verifies the half-list invariants.
pub fn test_half_neighbor_list(build: BuildFn) {
    // --- Case 1: small cubic system with cross-image pairs ---
    {
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
            pbc: [true, true, true],
        };
        let full = build(&sys, 3.0, false);
        let half = build(&sys, 3.0, true);
        verify_half_invariants(&full, &half);
        // Full has (0,1),(1,0),(0,3),(3,0) = 4; half keeps (0,1) and (0,3).
        assert_eq!(full.len(), 4);
        assert_eq!(half.len(), 2);
        let half_pairs: Vec<(usize, usize)> = half.iter().map(|n| (n.i, n.j)).collect();
        assert!(half_pairs.contains(&(0, 1)));
        assert!(half_pairs.contains(&(0, 3)));
    }

    // --- Case 2: single atom with face-adjacent self-images (cubic) ---
    {
        let sys = System {
            pos: vec![Vector3::new(0.0, 0.0, 0.0)],
            cell: UnitCell::new(
                Vector3::new(3.0, 0.0, 0.0),
                Vector3::new(0.0, 3.0, 0.0),
                Vector3::new(0.0, 0.0, 3.0),
            )
            .unwrap(),
            pbc: [true, true, true],
        };
        let full = build(&sys, 3.0, false);
        let half = build(&sys, 3.0, true);
        verify_half_invariants(&full, &half);
        // Full has 6 self-images (±a, ±b, ±c); half keeps +a, +b, +c.
        assert_eq!(full.len(), 6);
        assert_eq!(half.len(), 3);
        let half_offsets: Vec<[i32; 3]> = half.iter().map(|n| n.offset).collect();
        assert!(half_offsets.contains(&[1, 0, 0]));
        assert!(half_offsets.contains(&[0, 1, 0]));
        assert!(half_offsets.contains(&[0, 0, 1]));
    }

    // --- Case 3: full triclinic multi-atom ---
    {
        let cell = UnitCell::new(
            Vector3::new(6.0, 0.0, 0.0),
            Vector3::new(1.0, 5.0, 0.0),
            Vector3::new(0.5, 1.0, 4.0),
        )
        .unwrap();
        let sys = System {
            pos: vec![
                cell.get_cartesian(&Vector3::new(0.1, 0.1, 0.1)),
                cell.get_cartesian(&Vector3::new(0.9, 0.9, 0.9)),
                cell.get_cartesian(&Vector3::new(0.5, 0.2, 0.8)),
                cell.get_cartesian(&Vector3::new(0.3, 0.7, 0.4)),
            ],
            cell,
            pbc: [true, true, true],
        };
        let full = build(&sys, 4.0, false);
        let half = build(&sys, 4.0, true);
        verify_half_invariants(&full, &half);
        assert!(!half.is_empty());
    }

    // --- Case 4: slab geometry (mixed PBC) ---
    {
        let sys = System {
            pos: vec![
                Vector3::new(0.5, 0.5, 5.0),
                Vector3::new(3.0, 0.5, 5.0),
                Vector3::new(0.5, 3.0, 5.0),
            ],
            cell: UnitCell::new(
                Vector3::new(4.0, 0.0, 0.0),
                Vector3::new(0.0, 4.0, 0.0),
                Vector3::new(0.0, 0.0, 20.0),
            )
            .unwrap(),
            pbc: [true, true, false],
        };
        let full = build(&sys, 3.0, false);
        let half = build(&sys, 3.0, true);
        verify_half_invariants(&full, &half);
        assert!(!half.is_empty());
    }

    // --- Case 5: large cutoff forcing multi-image replicas (small cubic) ---
    {
        let sys = System {
            pos: vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 0.0, 0.0)],
            cell: UnitCell::new(
                Vector3::new(4.0, 0.0, 0.0),
                Vector3::new(0.0, 4.0, 0.0),
                Vector3::new(0.0, 0.0, 4.0),
            )
            .unwrap(),
            pbc: [true, true, true],
        };
        let full = build(&sys, 5.0, false);
        let half = build(&sys, 5.0, true);
        verify_half_invariants(&full, &half);
    }
}

#[cfg(test)]
mod tests {
    use crate::cell_list;
    use crate::naive;
    use crate::types::{System, UnitCell, Vector3};

    use super::sorted_neighbors;

    #[test]
    fn test_naive_cell_list_consistency() {
        // Full 3D triclinic, 4 atoms, cutoff=4.0
        // Verify naive and cell_list produce identical neighbor lists
        let cell = UnitCell::new(
            Vector3::new(6.0, 0.0, 0.0),
            Vector3::new(1.0, 5.0, 0.0),
            Vector3::new(0.5, 1.0, 4.0),
        )
        .unwrap();
        let sys = System {
            pos: vec![
                cell.get_cartesian(&Vector3::new(0.1, 0.1, 0.1)),
                cell.get_cartesian(&Vector3::new(0.9, 0.9, 0.9)),
                cell.get_cartesian(&Vector3::new(0.5, 0.2, 0.8)),
                cell.get_cartesian(&Vector3::new(0.3, 0.7, 0.4)),
            ],
            cell,
            pbc: [true, true, true],
        };
        let cutoff = 4.0;
        let naive_result = naive::build_neighbor_list(&sys, cutoff, false);
        let cell_list_result = cell_list::build_neighbor_list(&sys, cutoff, false);

        let naive_sorted = sorted_neighbors(&naive_result);
        let cell_list_sorted = sorted_neighbors(&cell_list_result);
        assert_eq!(
            naive_sorted, cell_list_sorted,
            "naive and cell_list produced different neighbor lists"
        );

        // Also verify distances match
        let mut naive_dists: Vec<(usize, usize, [i32; 3], i64)> = naive_result
            .iter()
            .map(|n| (n.i, n.j, n.offset, (n.distance * 1e10).round() as i64))
            .collect();
        naive_dists.sort();
        let mut cell_list_dists: Vec<(usize, usize, [i32; 3], i64)> = cell_list_result
            .iter()
            .map(|n| (n.i, n.j, n.offset, (n.distance * 1e10).round() as i64))
            .collect();
        cell_list_dists.sort();
        assert_eq!(
            naive_dists, cell_list_dists,
            "naive and cell_list produced different distances"
        );
    }

    #[test]
    fn test_naive_cell_list_half_consistency() {
        // Same triclinic system as the full-list consistency test: verify that
        // naive and cell_list produce identical half lists when invoked with
        // `half = true`.
        let cell = UnitCell::new(
            Vector3::new(6.0, 0.0, 0.0),
            Vector3::new(1.0, 5.0, 0.0),
            Vector3::new(0.5, 1.0, 4.0),
        )
        .unwrap();
        let sys = System {
            pos: vec![
                cell.get_cartesian(&Vector3::new(0.1, 0.1, 0.1)),
                cell.get_cartesian(&Vector3::new(0.9, 0.9, 0.9)),
                cell.get_cartesian(&Vector3::new(0.5, 0.2, 0.8)),
                cell.get_cartesian(&Vector3::new(0.3, 0.7, 0.4)),
            ],
            cell,
            pbc: [true, true, true],
        };
        let cutoff = 4.0;
        let naive_half = naive::build_neighbor_list(&sys, cutoff, true);
        let cell_list_half = cell_list::build_neighbor_list(&sys, cutoff, true);

        assert_eq!(
            sorted_neighbors(&naive_half),
            sorted_neighbors(&cell_list_half),
            "naive and cell_list produced different half neighbor lists"
        );

        let mut naive_dists: Vec<(usize, usize, [i32; 3], i64)> = naive_half
            .iter()
            .map(|n| (n.i, n.j, n.offset, (n.distance * 1e10).round() as i64))
            .collect();
        naive_dists.sort();
        let mut cell_list_dists: Vec<(usize, usize, [i32; 3], i64)> = cell_list_half
            .iter()
            .map(|n| (n.i, n.j, n.offset, (n.distance * 1e10).round() as i64))
            .collect();
        cell_list_dists.sort();
        assert_eq!(
            naive_dists, cell_list_dists,
            "naive and cell_list produced different half-list distances"
        );
    }
}
