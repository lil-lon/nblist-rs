# Step 1: Data Structures & Minimum Image

## Data Structures
- [x] `Vector3`: struct with x, y, z (f64)
  - [x] `new(x, y, z) -> Self`
  - [x] `norm(&self) -> f64`
  - [x] `impl Add<Vector3>` for `+`
  - [x] `impl Sub<Vector3>` for `-`
  - [x] `impl Mul<f64>` for scalar `*`
- [x] `UnitCell`: struct with a, b, c (f64, orthogonal only)
- [x] `System`: struct with `pos: Vec<Vector3>` and `cell: UnitCell`

## Minimum Image Convention
- [x] `UnitCell::minimum_image(&self, diff: Vector3) -> Vector3`
  - Apply PBC: each component wraps to [-L/2, L/2) (handles diff > L)
- [x] `System::distance(&self, i: usize, j: usize) -> f64`
  - Compute `pos[j] - pos[i]`, apply minimum image, return norm

## Validation
- [x] `#[test]` for Vector3 norm
- [x] `#[test]` for minimum image
- [x] `#[test]` for System::distance with PBC
- [x] `#[test]` for System::distance with large diff (> L)

# Step 2: Naive O(N^2) Neighbor List

- [x] `System::build_neighbor_list(&self, cutoff: f64) -> Vec<(usize, usize)>`
  - Double loop over all pairs (i < j)
  - Use System::distance, collect pairs within cutoff
- [x] `#[test]` for basic, all neighbors, no neighbors, PBC corner cases

# Step 3: Cell List (O(N))

- [ ] Divide box into cells of size >= cutoff
- [ ] Assign atoms to cells
- [ ] Search only neighboring cells (27 cells in 3D)
- [ ] `#[test]` verify same result as naive
- [ ] Benchmark and compare with naive

# Step 4: Triclinic Cell Support & Fractional Coordinates

- [ ] `UnitCell` with 3 lattice vectors (Vector3 x3)
- [ ] Fractional <-> Cartesian coordinate conversion
- [ ] Minimum image in fractional coordinates
- [ ] Update cell list for triclinic

# Step 5: Further Optimization

- [ ] Verlet list (skin distance, rebuild check)
- [ ] Data layout optimization (AoS vs SoA)
- [ ] Parallelization with rayon