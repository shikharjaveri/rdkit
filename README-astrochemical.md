# RDKit Astrochemical Extension

This extension to RDKit adds support for molecules with non-standard valences that are commonly found in astrochemistry research. These molecules would typically fail sanitization in standard RDKit, but this extension implements a relaxed sanitization approach that allows for their processing.

## Features

- Modified sanitization routine that skips strict valence checking (`sanitizeAstroMol`)
- Convenience functions for loading SMILES as astrochemical molecules (`MolFromSmilesAstro`)
- Support for unusual electron configurations found in space chemistry
- Compatibility with standard RDKit functionality (fingerprints, 3D conformer generation, etc.)

## Examples of Molecules Supported

- Polyyne anions (e.g., `C#CC#CC#CC#CC#[C-]`)
- Cyanide radicals (e.g., `[C]#N`)
- Carbon monoxide with unusual valence (e.g., `[C]=O`)
- Hypervalent carbon species (e.g., `[C](=O)=O`)
- Various astrochemical ions (e.g., `[CH+]=O`, `[NH+]#N`)

## Installation

### Building From Source

To use this extension, you need to build RDKit from source with the astrochemical modifications:

#### Prerequisites

- A C++ compiler (GCC, Clang, or MSVC)
- CMake (version 3.13.4 or higher)
- Boost libraries
- Python (and its development headers)
- NumPy
- Additional dependencies as per standard RDKit requirements

#### Build Steps (Linux/macOS)

1. Clone the repository:
   ```bash
   git clone https://github.com/shikharjaveri/rdkit.git
   cd rdkit
   git checkout explore-sanitization
   ```

2. Create a build directory:
   ```bash
   mkdir build
   cd build
   ```

3. Configure the build with CMake:
   ```bash
   cmake .. \
     -DCMAKE_INSTALL_PREFIX=/path/to/install \
     -DPYTHON_EXECUTABLE=/path/to/python \
     -DRDK_BUILD_PYTHON_WRAPPERS=ON \
     -DCMAKE_BUILD_TYPE=Release
   ```

4. Build and install:
   ```bash
   make -j4  # Adjust the number based on your CPU cores
   make install
   ```

#### Build Steps (Windows)

1. Clone the repository:
   ```powershell
   git clone https://github.com/shikharjaveri/rdkit.git
   cd rdkit
   git checkout explore-sanitization
   ```

2. Create a build directory:
   ```powershell
   mkdir build
   cd build
   ```

3. Configure the build with CMake:
   ```powershell
   cmake -G "Visual Studio 17 2022" -A x64 `
     -DRDK_BUILD_PYTHON_WRAPPERS=ON `
     -DRDK_BUILD_CPP_TESTS=ON `
     -DRDK_BUILD_INCHI_SUPPORT=ON ..
   ```

4. Build using Visual Studio or MSBuild:
   - Open the generated solution file in Visual Studio and build the ALL_BUILD and INSTALL targets
   - Or use MSBuild from the command line:
     ```powershell
     msbuild /p:Configuration=Release INSTALL.vcxproj
     ```

### Using a Virtual Environment (Recommended)

It's recommended to use a Python virtual environment when working with custom RDKit builds:

```bash
# Create a virtual environment
python -m venv rdkit_astro_env

# Activate the environment
# On Linux/macOS:
source rdkit_astro_env/bin/activate
# On Windows:
.\rdkit_astro_env\Scripts\activate

# Set environment variables (if needed)
export PYTHONPATH=/path/to/install/lib/python3.x/site-packages:$PYTHONPATH
```

## Usage

### Basic Usage

```python
from rdkit import Chem

# Using the new astrochemical function
mol = Chem.MolFromSmilesAstro("C#CC#CC#CC#CC#[C-]")

# Or use the workaround approach
mol = Chem.MolFromSmiles("C#CC#CC#CC#CC#[C-]", sanitize=False)
mol.UpdatePropertyCache(strict=False)
Chem.GetSymmSSSR(mol)
Chem.SetAromaticity(mol)
Chem.SetConjugation(mol)
Chem.SetHybridization(mol)
```

### Example Scripts

This repository includes several example scripts that demonstrate the usage of the astrochemical extension:

1. `astrochemical_molecules_demo.py` - Basic molecule handling and property calculation
2. `astrochemical_search_demo.py` - Substructure searching and fingerprint generation
3. `astrochemical_3d_visualization.py` - 3D structure generation and visualization

Run these scripts to see examples of how to work with astrochemical molecules:

```bash
python astrochemical_molecules_demo.py
```

## Implementation Details

The key modifications to RDKit include:

1. Added `ASTROCHEMICAL_SANITIZE` flag in MolOps.h
2. Implemented `sanitizeAstroMol` function in MolOps.cpp
3. Added Python wrapper functions (`MolFromSmilesAstro`, `SanitizeAstroMol`)
4. Added C++ literal operator (`_astrochemical`)

The astrochemical sanitization skips the following operations:
- Strict valence checking
- Radical assignment
- Hydrogen adjustment

While retaining other sanitization steps like:
- Aromaticity perception
- Conjugation setting
- Hybridization setting

## License

This extension is covered by the same BSD license as the main RDKit project.

## Acknowledgments

This extension is built upon the excellent work of the RDKit community and is intended to expand RDKit's capabilities for astrochemical research applications. 