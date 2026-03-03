# Nick Test Directory Organization

This directory contains tests for the SpaceHub AGN disk model development.

## Directory Structure

### `disk_model/` - Active Disk Model Tests (Jan-Feb 2026)
**Status:** Active development - disk model implementation is ongoing

Tests for the AGN disk model implementation (`src/interaction/disk-model.hpp` and `disk-model_damping.hpp`):
- `AutomatedDiskModelTest.cpp` - Component tests for disk model (vertical density, velocities, Mach regimes)
- `KeplerSimple-DiskModel.cpp` - Basic integration test for disk model
- `KeplerSimple-DiskModelSelfReg.cpp` - Self-regulating zone test
- `KeplerSimple-DiskModeRetrograde.cpp` - Retrograde orbit test
- `TestDiskVelocity.cpp` - Disk velocity calculation tests

**Output subdirectories:**
- `AnalyticalTests/` - Analytical test data (circular, retrograde, self-reg configurations)
- `disk_model_output/` - Component test CSVs (force components, Mach regimes, density profiles)

### `fairbairn_tests/` - Fairbairn Data Cube Tests (Feb 2026)
**Status:** Recently implemented - 4D interpolation for Fairbairn et al. drag forces

- `TestFairbairnDataCube.cpp` - Standalone test for 4D interpolation engine
  - Tests against Python reference implementation in `disktab/accessDataCube.ipynb`
  - Validates interpolation, extrapolation, and edge cases
- `fairbairn_output/` - Test output CSVs for comparison with Python

### `basic_tests/` - Basic Kepler Orbit Tests (Jan 2026)
**Status:** Foundational tests - StaticGasField force testing

- `KeplerSimple.cpp` - Basic Kepler orbit with static gas field
- `KeplerSimple-hc.cpp` - High-cadence variant for detailed output
- `KeplerSimple.dat` - Sample output data

### `archived/` - Parameter Sweep Tests (Dec 2025 - Jan 2026)
**Status:** Superseded by more recent disk model implementation

Early parameter sweep tests, likely replaced by newer disk model tests:
- `ParamSweepTest.cpp` - Old parameter sweep framework
- `KeplerSimpleSweep.cpp` - Sweep over binary parameters
- `sweep_results.txt`, `run_q*.dat` - Old sweep output data

## Building Tests

All executables have been removed - recompile as needed:

```bash
# From repository root:
g++ -std=c++17 -O3 test/nick_test/disk_model/KeplerSimple-DiskModel.cpp -o test/nick_test/KeplerSimple-DiskModel

# Or from test/nick_test:
g++ -std=c++17 -O3 fairbairn_tests/TestFairbairnDataCube.cpp -o TestFairbairnDataCube
```

## Notes on Implementation Status

- **Disk Model**: Active development - currently implementing analytical gradients in Self-Reg zone (see plan: `~/.claude/plans/concurrent-discovering-salamander.md`)
- **Fairbairn DataCube**: Recently implemented and tested against Python reference
- **Basic Tests**: Foundational - StaticGasField is implemented and working
- **Parameter Sweeps**: Archived - superseded by component-based disk model testing

---
*Last organized: 2026-02-16*
