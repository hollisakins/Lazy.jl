# Lazy.jl Implementation Plan & TODO List

## Phase 1: Critical Performance & Thread Safety Fixes 🚀
*Priority: URGENT - Performance bottlenecks and thread safety issues*

### 1.1 Thread Safety Fixes
- [x] **CRITICAL: Fix nested threading in template building** (`src/main.jl:416`) ✅ **COMPLETED**
  - Removed inner `@threads` loop to prevent thread pool exhaustion
  - This can cause thread pool exhaustion and performance degradation
  - **Impact**: Major performance improvement for large template sets

### 1.2 Memory Optimization
- [x] **Optimize transpose operations** (lines 501, 504, 548) ✅ **COMPLETED**
  - Replaced `transpose(templgrid_i)` with pre-allocated buffers and manual loops
  - **Impact**: Reduce memory allocations during fitting
- [x] **Add memory usage estimation/reporting** ✅ **COMPLETED**
  - Added comprehensive memory estimation with warnings
  - **Impact**: Prevent out-of-memory crashes
- [x] **Implement streaming for very large template sets** ✅ **SKIPPED**
  - Not needed - template files are small and current approach is efficient
  - **Impact**: N/A

### 1.3 Numerical Stability
- [x] **Add type annotations to performance-critical functions** ✅ **COMPLETED**
  - Added type annotations to `load_data()`, `set_sys_err()`, `get_filter()`, etc.
  - **Impact**: Better type inference and performance
- [x] **Add numerical stability checks** ✅ **COMPLETED**
  - Added checks for non-finite coefficients and invalid chi2 values
  - **Impact**: More robust fitting for edge cases

## Phase 2: Code Architecture & Maintainability 🏗️
*Priority: HIGH - Code quality and maintainability*

### 2.1 Modular Refactoring
- [ ] **Break down monolithic `main.jl`** (1000+ lines)
  - Extract `CLI.jl` - command line interface
  - Extract `Fitting.jl` - core fitting algorithms
  - Extract `Templates.jl` - template management
  - Extract `Filters.jl` - filter management
  - Extract `IO.jl` - file I/O operations
  - Extract `Utils.jl` - utility functions
  - **Impact**: Easier maintenance, testing, and contribution

### 2.2 CLI Architecture Improvement
- [ ] **Implement ArgParse.jl for robust argument parsing**
  - Replace manual CLI parsing with ArgParse.jl
  - Add automatic help generation and validation
  - Support for subcommands, optional arguments, and type checking
  - Better error messages for invalid arguments
  - **Impact**: More maintainable and extensible CLI, foundation for advanced features

### 2.3 Error Handling Consistency
- [x] **Replace `error()` calls with `LazyError`** (lines 239, 362, 432, 436)
  - Consistent error handling throughout codebase
  - Better error messages for users
- [ ] **Add file handle safety with try-finally blocks**
  - Ensure FITS files are properly closed on errors
  - **Impact**: Prevent resource leaks

### 2.4 Documentation Standardization
- [ ] **Convert docstrings to Julia format**
  - Replace Python-style `"""` with proper Julia docstrings
  - Add parameter and return type documentation
  - **Impact**: Better IDE support and API documentation

## Phase 3: Enhanced User Experience 👥
*Priority: MEDIUM - Usability for non-Julia users*

### 3.1 CLI Enhancements
- [ ] **Add configuration validation command**
  - `lazy validate-config <param_file>`
  - **Impact**: Help users catch configuration errors early
- [ ] **Add memory estimation command**
  - `lazy estimate-memory <param_file>`
  - **Impact**: Help users plan resource requirements
- [x] **Add version command** ✅ **COMPLETED**
  - `lazy --version` 
- [ ] **Add citation command**
  - `lazy --cite`
  - **Impact**: Better reproducibility and attribution

### 3.2 Better Error Messages
- [ ] **Improve error message quality**
  - Show available options when items not found
  - Suggest corrections for common mistakes
  - **Impact**: Better user experience for non-Julia users

### 3.3 Progress Reporting
- [ ] **Enhanced progress bars with time estimates**
  - Use ProgressMeter with ETA instead of basic ProgressBar
  - **Impact**: Better user feedback during long runs

## Phase 4: Testing & Validation 🧪
*Priority: MEDIUM - Reliability and correctness*

### 4.1 Comprehensive Test Suite
- [ ] **Add scientific computing tests**
  - Template interpolation accuracy
  - IGM attenuation correctness
  - Fitting algorithm validation
  - **Impact**: Catch regressions and validate scientific accuracy
- [ ] **Add integration tests with real data**
  - End-to-end fitting pipeline tests
  - **Impact**: Ensure real-world functionality

### 4.2 Performance Benchmarking
- [ ] **Create benchmark suite**
  - Track performance across versions
  - Identify performance regressions
  - **Impact**: Maintain performance standards

### 4.3 Configuration Validation
- [ ] **Comprehensive parameter validation**
  - Check for required sections and parameters
  - Validate file paths and data formats
  - **Impact**: Prevent runtime errors from bad configurations

## Phase 5: Advanced Features & Optimization 🔬
*Priority: LOW - Nice-to-have improvements*

### 5.1 Algorithm Optimization
- [ ] **SIMD optimization for hot loops**
  - Add `@simd` to template interpolation loops
  - **Impact**: Minor performance gains
- [ ] **Memory mapping for large template files**
  - Use memory-mapped files instead of loading everything
  - **Impact**: Support even larger datasets

### 5.2 Advanced CLI Features
- [ ] **Format conversion utilities**
  - `lazy convert-catalog <format>`
  - **Impact**: Easier migration from other tools
- [ ] **Configuration templates**
  - Generate config files for common use cases
  - **Impact**: Lower barrier to entry

## Implementation Notes

### Performance Impact Achieved/Estimates
- **Phase 1**: 🚀 ✅ **COMPLETED** - Major performance improvements achieved:
  - **Thread safety**: Fixed nested threading preventing thread pool exhaustion
  - **Memory optimization**: Eliminated expensive transpose operations in hot loops
  - **Memory management**: Added proactive memory estimation and warnings
  - **Numerical stability**: Added robust error detection and recovery
  - **Type inference**: Improved with comprehensive type annotations
- **Phase 2**: 🏗️ No performance impact, better maintainability
- **Phase 3**: 👥 Minor performance impact, major usability improvement
- **Phase 4**: 🧪 No performance impact, better reliability
- **Phase 5**: 🔬 Minor performance improvements

### Resource Requirements
- **Phase 1**: ~1-2 weeks, requires deep Julia threading knowledge
- **Phase 2**: ~2-3 weeks, major refactoring effort
- **Phase 3**: ~1 week, straightforward CLI improvements
- **Phase 4**: ~1-2 weeks, requires test data and validation
- **Phase 5**: ~1 week, optional enhancements

### Compatibility Notes
- All phases maintain backward compatibility for CLI interface
- Parameter file format remains unchanged
- Output format remains unchanged
- Only internal code organization changes

### Quick Wins (Can be done immediately)
- [x] Fix nested threading issue (Phase 1.1) ✅ **COMPLETED**
- [x] Add memory usage warnings (Phase 1.2) ✅ **COMPLETED** 
- [x] Fix error() -> LazyError calls (Phase 2.3) ✅ **COMPLETED**
- [x] Add --version flag (Phase 3.1) ✅ **COMPLETED**