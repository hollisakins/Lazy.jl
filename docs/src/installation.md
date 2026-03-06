# Installation

Lazy.jl requires Julia 1.12+ and is installed from source.

## Prerequisites

Install Julia 1.12+ using [juliaup](https://github.com/JuliaLang/juliaup):

```bash
curl -fsSL https://install.julialang.org | sh
juliaup add 1.12
juliaup default 1.12
```

## Install Lazy.jl

Clone the repository and run the install script:

```bash
git clone https://github.com/hollisakins/Lazy.jl.git
cd Lazy.jl
bash install.sh
```

The install script:
1. Verifies Julia 1.12+ is available
2. Installs the Lazy.jl package into a dedicated Julia environment (`@lazy`)
3. Copies the `lazy` CLI executable to `~/.julia/bin/`

## Add to PATH

Add the following to your shell configuration (`~/.bashrc` or `~/.zshrc`):

```bash
export PATH="$PATH:$HOME/.julia/bin"
```

Then reload your shell:

```bash
source ~/.bashrc  # or source ~/.zshrc
```

## Verify Installation

```bash
lazy --version
lazy list-templates
```

The first run will precompile packages (1-2 minutes). Subsequent runs start instantly.

## Updating

To update to the latest version:

```bash
cd Lazy.jl
git pull
bash install.sh
```

## Troubleshooting

**First run is slow**: Julia precompiles all dependencies on first use. This is a one-time cost.

**`lazy: command not found`**: Make sure `~/.julia/bin` is in your PATH. Check with `echo $PATH`.

**Julia version mismatch**: Run `julia --version` to verify you have 1.12+. Use `juliaup default 1.12` to set the default.
