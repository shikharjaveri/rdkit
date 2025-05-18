#!/usr/bin/env python
"""
Test script to evaluate astrochemical molecule functionality with a polyyne anion
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import sys

def simulate_astro_sanitize(mol):
    """
    Python implementation that simulates our C++ sanitizeAstroMol function
    Uses only the operations that are safe for astrochemical molecules
    """
    # First, make sure we're working with a writable molecule
    if not isinstance(mol, Chem.RWMol):
        mol = Chem.RWMol(mol)
    
    # Update property cache without checking valence - this is key for unusual molecules
    mol.UpdatePropertyCache(strict=False)
    
    # Symmetrize rings (should be safe for unusual valences)
    try:
        Chem.GetSymmSSSR(mol)
    except Exception as e:
        print(f"  Warning: GetSymmSSSR failed - {str(e)}")
    
    # Set aromaticity (should be safe)
    try:
        Chem.SetAromaticity(mol)
    except Exception as e:
        print(f"  Warning: SetAromaticity failed - {str(e)}")
    
    # Set conjugation (should be safe)
    try:
        Chem.SetConjugation(mol)
    except Exception as e:
        print(f"  Warning: SetConjugation failed - {str(e)}")
    
    # Set hybridization (should be safe)
    try:
        Chem.SetHybridization(mol)
    except Exception as e:
        print(f"  Warning: SetHybridization failed - {str(e)}")
    
    # Update property cache again
    mol.UpdatePropertyCache(strict=False)
    
    return mol

def test_polyyne_anion():
    """Test astrochemical molecule handling with a polyyne anion"""
    
    # The polyyne anion to test
    smiles = "C#CC#CC#CC#CC#[C-]"
    
    print(f"\nTesting polyyne anion: {smiles}")
    print("=" * 60)
    
    # Test 1: Standard sanitization (should fail)
    print("\n1. STANDARD SANITIZATION")
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            print("  Standard sanitization: SUCCEEDED (unexpected!)")
            print("  Molecule details:")
            print(f"    - Num atoms: {mol.GetNumAtoms()}")
            print(f"    - Num bonds: {mol.GetNumBonds()}")
            print(f"    - Formula: {Chem.rdMolDescriptors.CalcMolFormula(mol)}")
        else:
            print("  Standard sanitization: FAILED (mol is None)")
    except Exception as e:
        print(f"  Standard sanitization: FAILED (expected) - {str(e)}")
    
    # Test 2: Simulated astrochemical sanitization
    print("\n2. SIMULATED ASTROCHEMICAL SANITIZATION")
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol:
            mol = simulate_astro_sanitize(mol)
            print("  Simulated astro sanitization: SUCCEEDED")
            print("  Molecule details:")
            print(f"    - Num atoms: {mol.GetNumAtoms()}")
            print(f"    - Num bonds: {mol.GetNumBonds()}")
            print(f"    - Formula: {Chem.rdMolDescriptors.CalcMolFormula(mol)}")
            
            # Try to generate 2D coordinates for visualization
            try:
                AllChem.Compute2DCoords(mol)
                print("  2D coordinate generation: SUCCEEDED")
                # Save molecule image
                img = Draw.MolToImage(mol, size=(400, 200))
                img.save('polyyne_anion_2d.png')
                print("  Saved 2D structure as 'polyyne_anion_2d.png'")
            except Exception as e:
                print(f"  2D coordinate generation: FAILED - {str(e)}")
            
            # Try to generate 3D coordinates
            try:
                status = AllChem.EmbedMolecule(mol)
                if status == 0:
                    print("  3D coordinate generation: SUCCEEDED")
                    # Save molecule as MOL file
                    with open('polyyne_anion_3d.mol', 'w') as f:
                        f.write(Chem.MolToMolBlock(mol))
                    print("  Saved 3D structure as 'polyyne_anion_3d.mol'")
                else:
                    print("  3D coordinate generation: FAILED")
            except Exception as e:
                print(f"  3D coordinate generation: FAILED - {str(e)}")
            
        else:
            print("  Simulated astro sanitization: FAILED (mol is None)")
    except Exception as e:
        print(f"  Simulated astro sanitization: FAILED - {str(e)}")
    
    # Test 3: Direct astrochemical functions
    print("\n3. DIRECT ASTROCHEMICAL FUNCTIONS")
    
    # Test SmilesToAstroMol
    try:
        mol = Chem.SmilesToAstroMol(smiles)
        if mol:
            print("  SmilesToAstroMol: SUCCEEDED")
            print("  Molecule details:")
            print(f"    - Formula: {Chem.rdMolDescriptors.CalcMolFormula(mol)}")
            
            # Try to generate 3D coordinates
            try:
                status = AllChem.EmbedMolecule(mol)
                if status == 0:
                    print("  3D coordinate generation: SUCCEEDED")
                else:
                    print("  3D coordinate generation: FAILED")
            except Exception as e:
                print(f"  3D coordinate generation: FAILED - {str(e)}")
        else:
            print("  SmilesToAstroMol: FAILED (returned None)")
    except AttributeError:
        print("  SmilesToAstroMol: FUNCTION NOT AVAILABLE (needs C++ implementation)")
    except Exception as e:
        print(f"  SmilesToAstroMol: FAILED - {str(e)}")
    
    # Test SanitizeAstroMol
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol:
            Chem.SanitizeAstroMol(mol)
            print("  SanitizeAstroMol: SUCCEEDED")
            print(f"    - Formula: {Chem.rdMolDescriptors.CalcMolFormula(mol)}")
        else:
            print("  SanitizeAstroMol: FAILED (mol is None)")
    except AttributeError:
        print("  SanitizeAstroMol: FUNCTION NOT AVAILABLE (needs C++ implementation)")
    except Exception as e:
        print(f"  SanitizeAstroMol: FAILED - {str(e)}")
    
    print("\nTEST COMPLETED")

if __name__ == "__main__":
    test_polyyne_anion() 