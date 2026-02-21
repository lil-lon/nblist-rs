use std::ops::{Add, Mul, Sub};

fn main() {
    println!("Hello, world!");
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct Vector3 {
    x: f64,
    y: f64,
    z: f64,
}

impl Add for Vector3 {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl Sub for Vector3 {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl Mul<f64> for Vector3 {
    type Output = Self;
    fn mul(self, other: f64) -> Self {
        Self {
            x: self.x * other,
            y: self.y * other,
            z: self.z * other,
        }
    }
}

impl Vector3 {
    fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    fn norm(&self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }
}

#[derive(Debug, Clone, Copy)]
struct UnitCell {
    a: f64,
    b: f64,
    c: f64,
}

// Positions are stored in Cartesian coordinates.
#[derive(Debug, Clone)]
struct System {
    pos: Vec<Vector3>,
    cell: UnitCell,
}

impl System {
    fn build_neighbor_list(&self, cutoff: f64) -> Vec<Neighbor> {
        let mut result = Vec::new();
        let a_replicas = (cutoff / self.cell.a).ceil() as i32;
        let b_replicas = (cutoff / self.cell.b).ceil() as i32;
        let c_replicas = (cutoff / self.cell.c).ceil() as i32;
        let n = self.pos.len();
        for i in 0..n {
            for j in 0..n {
                for a_ind in -a_replicas..=a_replicas {
                    for b_ind in -b_replicas..=b_replicas {
                        for c_ind in -c_replicas..=c_replicas {
                            if i == j && (a_ind, b_ind, c_ind) == (0, 0, 0) {
                                continue;
                            }
                            let pos_j = self.pos[j]
                                + Vector3::new(
                                    a_ind as f64 * self.cell.a,
                                    b_ind as f64 * self.cell.b,
                                    c_ind as f64 * self.cell.c,
                                );
                            let distance = (pos_j - self.pos[i]).norm();
                            if distance < cutoff {
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
}

#[derive(Debug)]
struct Neighbor {
    i: usize,
    j: usize,
    offset: [i32; 3],
    distance: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_vector3_norm() {
        let vec = Vector3::new(3.0, 4.0, 5.0);
        let result = vec.norm();
        let expected = 50.0_f64.sqrt();
        assert!((result - expected).abs() < 1e-10);
    }

    #[test]
    fn test_vector3_add() {
        let a = Vector3::new(1.0, 2.0, 3.0);
        let b = Vector3::new(4.0, 5.0, 6.0);
        assert_eq!(a + b, Vector3::new(5.0, 7.0, 9.0));
    }

    #[test]
    fn test_vector3_sub() {
        let a = Vector3::new(1.0, 2.0, 3.0);
        let b = Vector3::new(4.0, 8.0, 16.0);
        assert_eq!(b - a, Vector3::new(3.0, 6.0, 13.0));
    }

    #[test]
    fn test_vector3_mul() {
        let a = Vector3::new(1.0, 2.0, 3.0);
        assert_eq!(a * 5.5, Vector3::new(5.5, 11.0, 16.5));
    }

    #[test]
    fn test_neighbor_list_basic() {
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
        let result = sys.build_neighbor_list(3.0);
        // Both directions: (0,1), (1,0), (0,3), (3,0)
        let pairs: Vec<(usize, usize)> = result.iter().map(|n| (n.i, n.j)).collect();
        assert!(pairs.contains(&(0, 1)));
        assert!(pairs.contains(&(1, 0)));
        assert!(pairs.contains(&(0, 3)));
        assert!(pairs.contains(&(3, 0)));
        assert_eq!(result.len(), 4);
    }

    #[test]
    fn test_neighbor_list_no_neighbors() {
        let sys = System {
            pos: vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(5.0, 5.0, 5.0)],
            cell: UnitCell {
                a: 20.0,
                b: 20.0,
                c: 20.0,
            },
        };
        let result = sys.build_neighbor_list(1.0);
        assert!(result.is_empty());
    }

    #[test]
    fn test_neighbor_list_pbc_corner() {
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
        let result = sys.build_neighbor_list(2.0);
        // Both directions
        let pairs: Vec<(usize, usize)> = result.iter().map(|n| (n.i, n.j)).collect();
        assert_eq!(pairs.len(), 2);
        assert!(pairs.contains(&(0, 1)));
        assert!(pairs.contains(&(1, 0)));
    }

    #[test]
    fn test_neighbor_list_self_image() {
        // Single atom at origin, L=3, cutoff=3.5
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
        let result = sys.build_neighbor_list(3.5);
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

    #[test]
    fn test_neighbor_list_multiple_images() {
        // 2 atoms in L=4 box, cutoff=5
        // A(0,0,0) B(1,0,0)
        // cutoff > L/2, so multiple images of B are neighbors of A:
        //   offset [0,0,0]: dist=1  ✓
        //   offset [-1,0,0]: B at (-3,0,0), dist=3  ✓
        //   offset [1,0,0]: B at (5,0,0), dist=5 → not < cutoff ✗
        //   offset [0,1,0]: B at (1,4,0), dist=sqrt(17)≈4.12  ✓
        //   offset [0,-1,0]: B at (1,-4,0), dist=sqrt(17)≈4.12  ✓
        //   offset [0,0,1]: B at (1,0,4), dist=sqrt(17)≈4.12  ✓
        //   offset [0,0,-1]: B at (1,0,-4), dist=sqrt(17)≈4.12  ✓
        // Also A self-images at distance 4.0 along each axis:
        //   (±4,0,0), (0,±4,0), (0,0,±4) → dist=4  ✓ (6 total)
        let sys = System {
            pos: vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 0.0, 0.0)],
            cell: UnitCell {
                a: 4.0,
                b: 4.0,
                c: 4.0,
            },
        };
        let result = sys.build_neighbor_list(5.0);

        // Check A→B has multiple offsets
        let a_to_b: Vec<&Neighbor> = result.iter().filter(|n| n.i == 0 && n.j == 1).collect();
        // [0,0,0] dist=1, [-1,0,0] dist=3,
        // [0,1,0], [0,-1,0], [0,0,1], [0,0,-1] each dist=√17≈4.12 → 6 total
        assert_eq!(a_to_b.len(), 6);
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

    #[test]
    fn test_neighbor_list_non_cubic() {
        // Non-cubic cell: a=3, b=10, c=10, cutoff=3.5
        // A(0,0,0) B(1,0,0)
        // Along a-axis (L=3): replicas needed. Along b,c (L=10): no replicas.
        // A→B: offset [0,0,0] dist=1 ✓, offset [-1,0,0] dist=|1-3|=2 ✓, offset [1,0,0] dist=|1+3|=4 ✗
        // A→A self-image: offset [±1,0,0] dist=3 ✓, offset [0,±1,0] dist=10 ✗
        let sys = System {
            pos: vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 0.0, 0.0)],
            cell: UnitCell {
                a: 3.0,
                b: 10.0,
                c: 10.0,
            },
        };
        let result = sys.build_neighbor_list(3.5);

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
}
