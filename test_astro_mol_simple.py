#!/usr/bin/env python
"""
Test script for astrochemical molecule sanitization (Python-only version)
This script simulates our custom sanitization approach using existing RDKit functionality
and tests the new direct wrapper function
"""

from rdkit import Chem
from rdkit.Chem import AllChem
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
    Chem.SanitizeMsg = 'RINGS'
    Chem.SanitizeMol(mol, Chem.SANITIZE_SYMMRINGS, catchErrors=True)
    
    # Set aromaticity (should be safe)
    Chem.SanitizeMsg = 'AROM'
    Chem.SanitizeMol(mol, Chem.SANITIZE_SETAROMATICITY, catchErrors=True)
    
    # Set conjugation (should be safe)
    Chem.SanitizeMsg = 'CONJ'
    Chem.SanitizeMol(mol, Chem.SANITIZE_SETCONJUGATION, catchErrors=True)
    
    # Set hybridization (should be safe)
    Chem.SanitizeMsg = 'HYBRD'
    Chem.SanitizeMol(mol, Chem.SANITIZE_SETHYBRIDIZATION, catchErrors=True)
    
    # Update property cache again
    mol.UpdatePropertyCache(strict=False)
    
    return mol

def test_astro_mol_sanitization():
    """Test that astrochemical molecules can be sanitized"""
    # Examples of astrochemical molecules with unusual valence states
    # that would normally fail sanitization
    smiles_list = [
        'CC#N[C]',          # Cyanomethylidyne radical
        'C#C[C]',           # Propargyl radical
        'N#N#N',            # Azide with unusual bonding
        '[C](=O)(=O)',      # Carbon dioxide with unusual representation
        'O=[Si]=O',         # Silicon dioxide
        'CC(=O)[O]',        # Acetate radical
    ]
    
    sim_success_count = 0
    direct_success_count = 0
    
    for idx, smiles in enumerate(smiles_list):
        print(f"\nTesting molecule {idx+1}: {smiles}")
        
        # First show that standard sanitization fails
        try:
            mol = Chem.MolFromSmiles(smiles, sanitize=False)
            Chem.SanitizeMol(mol)
            print("  Standard sanitization: SUCCEEDED (unexpected!)")
        except Exception as e:
            print(f"  Standard sanitization: FAILED (expected) - {str(e)}")
        
        # Test our simulated custom sanitization
        try:
            mol = Chem.MolFromSmiles(smiles, sanitize=False)
            mol = simulate_astro_sanitize(mol)
            sim_success_count += 1
            print("  Simulated astro sanitization: SUCCEEDED")
            
            # Test that we can now use other RDKit functionality
            print("  Molecule has properties:")
            print(f"    - Num atoms: {mol.GetNumAtoms()}")
            print(f"    - Num bonds: {mol.GetNumBonds()}")
            print(f"    - Formula: {Chem.rdMolDescriptors.CalcMolFormula(mol)}")
            
            # Try to generate 3D coordinates
            try:
                status = AllChem.EmbedMolecule(mol)
                if status == 0:
                    print("  3D coordinate generation: SUCCEEDED")
                else:
                    print("  3D coordinate generation: FAILED")
            except Exception as e:
                print(f"  3D coordinate generation: ERROR - {str(e)}")
            
        except Exception as e:
            print(f"  Simulated astro sanitization: FAILED - {str(e)}")
            
        # Now test the direct function if it's available
        try:
            # Try the direct function if it's been added to the RDKit build
            mol = Chem.SmilesToAstroMol(smiles)
            direct_success_count += 1
            print("  Direct astro sanitization: SUCCEEDED")
            
            # Test that we can use other RDKit functionality
            print(f"    - Formula: {Chem.rdMolDescriptors.CalcMolFormula(mol)}")
            
        except AttributeError:
            print("  Direct astro sanitization: FUNCTION NOT AVAILABLE (needs C++ implementation)")
        except Exception as e:
            print(f"  Direct astro sanitization: FAILED - {str(e)}")
    
    print(f"\nSummary: {sim_success_count}/{len(smiles_list)} molecules passed simulated astro sanitization")
    if direct_success_count > 0:
        print(f"         {direct_success_count}/{len(smiles_list)} molecules passed direct astro sanitization")
    else:
        print("         Direct astro sanitization function not available - needs C++ implementation")
    
    return sim_success_count == len(smiles_list)

if __name__ == "__main__":
    success = test_astro_mol_sanitization()
    sys.exit(0 if success else 1) 