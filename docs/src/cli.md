# CLI Reference

## Commands

### `lazy fit`

Run photometric redshift fitting.

```bash
lazy fit -p <param_file> [options]
```

| Flag | Description |
|------|-------------|
| `-p`, `--param` | Path to the TOML parameter file (required) |
| `-t`, `--threads` | Number of threads: an integer or `auto` for all available cores. Default: 1 |
| `-y`, `--yes` | Non-interactive mode. Automatically resumes interrupted jobs without prompting |
| `-h`, `--help` | Show help for the fit command |

**Examples:**

```bash
lazy fit -p params.toml              # Single-threaded
lazy fit -p params.toml -t 8         # 8 threads
lazy fit -p params.toml -t auto      # All available threads
lazy fit -p params.toml -y           # Auto-resume without prompts
```

### `lazy list-templates`

Display all available template sets and their template counts.

```bash
lazy list-templates
```

### `lazy list-filters`

Display all available filter transmission curves with descriptions.

```bash
lazy list-filters
```

### `lazy params`

Print an annotated example parameter file to stdout. Redirect to create a starting configuration:

```bash
lazy params > my_params.toml
```

### `lazy cache-clear`

Remove all cached template grids. Useful if you've modified template files or want to reclaim disk space.

```bash
lazy cache-clear
```

### Global Options

| Flag | Description |
|------|-------------|
| `-v`, `--version` | Show version information |
| `-h`, `--help` | Show help message |

## Threading

Thread count is controlled by the `-t` flag on the `lazy fit` command. This is passed to Julia's `-t` option internally.

- `-t 1` (default): Single-threaded execution
- `-t N`: Use N threads
- `-t auto`: Use all available CPU threads

Threading applies to both template grid construction and the per-object fitting loop. See [Advanced Usage](@ref) for performance scaling benchmarks.
