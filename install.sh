#!/bin/bash

# Script to check if Julia is installed on the system

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Define the executable to copy
EXECUTABLE="bin/lazy"
TARGET_DIR=~/.julia/bin


echo "Running environment checks..."

# Check if current directory is named Lazy.jl
current_dir=$(basename "$(pwd)")
echo "Current directory: $current_dir"

if [ "$current_dir" != "Lazy.jl" ]; then
    echo "✗ ERROR: Current directory is not Lazy.jl" >&2
    echo "Please run this script from the Lazy.jl directory" >&2
    exit 1
else
    echo "✓ Current directory is Lazy.jl"
fi

echo "Checking for Julia installation..."

# First check if julia is installed 
if command_exists julia; then
    julia_version=$(julia --version)
    # Check if version contains "1.12"
    echo "✓ Julia is installed: $julia_version"
    if [[ $julia_version == *"1.12"* ]]; then
        echo "✓ Required version 1.12 is installed"
    else
        echo "✗ ERROR: Julia version 1.12 is required" >&2
        echo "Current version doesn't match the required version" >&2
        exit 1
    fi
else
    echo "✗ ERROR: Julia is not installed or not in PATH" >&2
    echo "Please install Julia from https://julialang.org/install/" >&2
    exit 1
fi

echo "Installing Lazy.jl..."
julia --project=@lazy -e 'using Pkg; try Pkg.rm("Lazy") catch; "" end; Pkg.develop(path=".")'
if [ $? -eq 0 ]; then
    echo "✓ Successfully installed Lazy.jl"
else
    echo "✗ ERROR: Failed to install Lazy.jl" >&2
    exit 1
fi


# Copy executable to ~/.julia/bin
echo "Attempting to copy $EXECUTABLE to $TARGET_DIR..."

# Check if target directory exists
if [ ! -d "$TARGET_DIR" ]; then
    echo "Directory $TARGET_DIR does not exist"
    read -p "Would you like to create it? (y/n): " choice
    case "$choice" in
        y|Y|yes|Yes|YES)
            mkdir -p "$TARGET_DIR"
            echo "✓ Created directory $TARGET_DIR"
            ;;
        *)
            echo "✗ Installation aborted" >&2
            exit 1
            ;;
    esac
fi

# Check if executable exists
if [ ! -f "$EXECUTABLE" ]; then
    echo "✗ ERROR: Executable '$EXECUTABLE' not found" >&2
    exit 1
fi

# Copy the executables
TARGET_EXECUTABLE="$TARGET_DIR/$(basename "$EXECUTABLE")"
cp "$EXECUTABLE" "$TARGET_DIR/"
if [ $? -eq 0 ]; then
    echo "✓ Successfully copied $EXECUTABLE to $TARGET_DIR"
    
    # Make it executable
    chmod +x "$TARGET_EXECUTABLE"
    echo "✓ Made $TARGET_EXECUTABLE executable"
    
    # Check if the directory is in PATH
    if [[ ":$PATH:" != *":$TARGET_DIR:"* ]]; then
        echo "NOTE: $TARGET_DIR is not in your PATH"
        echo "You may want to add it by adding this line to your shell profile:"
        echo "export PATH=\"\$PATH:$TARGET_DIR\""
    fi
    
    exit 0
else
    echo "✗ ERROR: Failed to copy $EXECUTABLE to $TARGET_DIR" >&2
    exit 1
fi
