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

## Validation
- [x] `#[test]` for Vector3 norm, add, sub, mul

# Step 2: Naive O(N^2) Neighbor List

- [x] `System::build_neighbor_list(&self, cutoff: f64) -> Vec<(usize, usize)>`
  - Double loop over all pairs (i < j)
  - Use System::distance, collect pairs within cutoff
- [x] `#[test]` for basic, no neighbors, PBC corner cases

# Step 3: Periodic Image Replication (cutoff > L/2)

- [x] Replicate periodic images: replica count per axis = ceil(cutoff / L)
- [x] Update build_neighbor_list to iterate over all periodic offsets
- [x] Return `Vec<Neighbor>` with indices, offset, and distance (both directions)
- [x] Remove minimum_image and distance (no longer needed)
- [x] `#[test]` self-image, multiple images, non-cubic cell

# Step 4: Cell List (O(N))

- [x] Divide box into cells of size >= cutoff
- [x] Assign atoms to cells
- [x] Search only neighboring cells (27 cells in 3D)
- [x] `#[test]` verify same result as naive
- [x] Benchmark and compare with naive

# Step 5: Triclinic Cell Support & Fractional Coordinates

- [ ] `UnitCell` with 3 lattice vectors (Vector3 x3)
- [ ] Fractional <-> Cartesian coordinate conversion
- [ ] Minimum image in fractional coordinates
- [ ] Update cell list for triclinic

# Step 6: Further Optimization

- [ ] Verlet list (skin distance, rebuild check)
- [ ] Data layout optimization (AoS vs SoA)
- [ ] Parallelization with rayon