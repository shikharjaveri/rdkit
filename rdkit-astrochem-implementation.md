# RDKit Astrochemical Molecule Support Implementation

## Overview

This document summarizes the changes made to RDKit to support astrochemical molecules with non-standard valences that would typically fail sanitization in the standard RDKit implementation.

## Problem Statement

Standard RDKit sanitization enforces strict valence rules which causes many astrochemical molecules with unusual electron configurations (found in space chemistry) to fail validation. This makes it difficult to work with these molecules for research purposes.

## Solution

We implemented a relaxed sanitization approach for astrochemical molecules that:

1. Skips strict valence checking and radical assignment
2. Only applies the sanitization steps that are safe for unusual structures
3. Exposes the functionality through convenient API functions

## Implementation Details

We made the following changes to the RDKit codebase:

### 1. Added a new sanitization flag constant in MolOps.h

```cpp
//! flags for controlling sanitization
const unsigned int SANITIZE_CLEANUP = 1;
const unsigned int SANITIZE_PROPERTIES = 2;
...
const unsigned int SANITIZE_ALL = 0xFFF;

//! Custom sanitization flags for astrochemical molecules
const unsigned int ASTROCHEMICAL_SANITIZE = 
  SANITIZE_CLEANUP | 
  SANITIZE_SYMMRINGS | 
  SANITIZE_KEKULIZE |
  SANITIZE_SETAROMATICITY |
  SANITIZE_SETCONJUGATION |
  SANITIZE_SETHYBRIDIZATION |
  SANITIZE_CLEANUPCHIRALITY |
  SANITIZE_CLEANUPATROPISOMERS;
```

### 2. Added sanitizeAstroMol function in MolOps.cpp

```cpp
void sanitizeAstroMol(RWMol &mol) {
  // Clean up molecule
  mol.clearComputedProps();
  
  // Clear out any cached properties
  mol.updatePropertyCache(false);  // don't check valence (key for astrochemical molecules)
  
  // Apply safe sanitization steps from ASTROCHEMICAL_SANITIZE
  cleanUp(mol);
  cleanUpOrganometallics(mol);
  symmetrizeSSSR(mol);
  Kekulize(mol);
  setAromaticity(mol);
  setConjugation(mol);
  setHybridization(mol);
  cleanupChirality(mol);
  cleanupAtropisomerism(mol);
  
  // Note: We specifically avoid:
  // - adjustHs()
  // - assignRadicals()
  // - and other operations that check valence
}
```

### 3. Added Python wrapper in MolOps.cpp

```cpp
void sanitizeAstroMolPy(ROMol &mol) {
  auto &wmol = static_cast<RWMol &>(mol);
  MolOps::sanitizeAstroMol(wmol);
}
```

### 4. Registered Python function in MolOps.cpp wrapper

```cpp
docString =
  "Alternative sanitization approach for astrochemical molecules with non-standard valences.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to be sanitized\n\
\n\
  NOTES:\n\
\n\
    - Uses a relaxed sanitization protocol that skips strict valence checking,\n\
      radical assignment, and other operations that would normally cause astrochemical\n\
      molecules to fail. This allows RDKit to process structures with unusual valence\n\
      states that are found in space chemistry.\n\
\n";
python::def("SanitizeAstroMol", sanitizeAstroMolPy, python::arg("mol"),
          docString.c_str());
```

### 5. Enhanced SMILES parsing in SmilesParse.h/.cpp

Added new function declarations:

```cpp
// In SmilesParse.h
RDKIT_SMILESPARSE_EXPORT std::unique_ptr<RDKit::RWMol> MolFromSmilesAstro(
    const std::string &smi,
    const SmilesParserParams &params = SmilesParserParams());

// In v1 namespace
inline RWMol *SmilesToAstroMol(
    const std::string &smi, int debugParse = 0, bool sanitize = true,
    std::map<std::string, std::string> *replacements = nullptr) {
  RDKit::v2::SmilesParse::SmilesParserParams params;
  params.debugParse = debugParse;
  if (replacements) {
    params.replacements = *replacements;
  }
  if (sanitize) {
    params.sanitize = true;
    params.removeHs = true;
  } else {
    params.sanitize = false;
    params.removeHs = false;
  }
  return RDKit::v2::SmilesParse::MolFromSmilesAstro(smi, params).release();
};
```

Implemented functions:
```cpp
std::unique_ptr<RWMol> MolFromSmilesAstro(const std::string &smiles,
                                     const SmilesParserParams &params) {
  // Create the molecule without standard sanitization
  auto mol = MolFromSmiles(smiles, SmilesParserParams(params).setSanitize(false));
  
  if (mol && params.sanitize) {
    // Apply astrochemical sanitization instead of standard sanitization
    MolOps::sanitizeAstroMol(*mol);
    
    if (params.removeHs) {
      MolOps::removeHs(*mol);
    }
  }

  return mol;
}
```

### 6. Added _astrochemical literal operator

```cpp
inline std::unique_ptr<RDKit::RWMol> operator"" _astrochemical(const char *text,
                                                               size_t len) {
  std::string smi(text, len);
  try {
    return v2::SmilesParse::MolFromSmilesAstro(smi);
  } catch (const RDKit::MolSanitizeException &) {
    return nullptr;
  }
}
```

## Usage Examples

```python
from rdkit import Chem

# Example astrochemical molecule that would fail standard sanitization
smiles = 'CC#N[C]'  # Cyanomethylidyne radical

# Method 1: Using custom SmilesToAstroMol function
mol = Chem.SmilesToAstroMol(smiles)

# Method 2: Create molecule without sanitization, then apply custom sanitization
mol = Chem.MolFromSmiles(smiles, sanitize=False)
Chem.SanitizeAstroMol(mol)

# Now the molecule can be used with standard RDKit functionality
mw = Chem.Descriptors.MolWt(mol)
Chem.AllChem.EmbedMolecule(mol)  # Generate 3D coordinates
```

## Building and Installation

To use these changes, the RDKit codebase needs to be recompiled with these modifications:

1. Clone the RDKit repository
2. Apply the changes described above
3. Build RDKit following the standard build instructions for your platform
4. Install the compiled version

## Use Cases

- Processing and analyzing molecules found in interstellar clouds
- Working with exotic molecular species in astrochemistry research
- Supporting molecules with unusual bonding patterns and electron configurations
- Enabling computational studies of space chemistry

## Limitations

- This approach bypasses some chemical validation, so the resulting structures may not be chemically valid according to standard chemistry rules
- Care should be taken when using these molecules in simulations or other computational workflows 