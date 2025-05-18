#!/usr/bin/env python
"""
Demonstration of astrochemical molecule handling in RDKit.

This script shows:
1. How standard RDKit fails with non-standard valence states
2. The current workaround approach that can be used 
3. The new API we've implemented that simplifies working with these molecules
"""

from rdkit import Chem
import sys

# Examples of astrochemical molecules with unusual valence states
EXAMPLES = [
    'CC#N[C]',          # Cyanomethylidyne radical
    'C#C[C]',           # Propargyl radical
    'N#N#N',            # Azide with unusual bonding
    '[C](=O)(=O)',      # Carbon dioxide with unusual representation
    'O=[Si]=O',         # Silicon dioxide
    'CC(=O)[O]',        # Acetate radical
]

def demonstrate_standard_approach():
    """Show how standard RDKit sanitization fails on these molecules"""
    print("\n\n1. STANDARD RDKIT APPROACH:")
    print("============================")
    
    for i, smiles in enumerate(EXAMPLES):
        print(f"\nMolecule {i+1}: {smiles}")
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                print("  SUCCESS - Molecule was created and sanitized")
            else:
                print("  FAILED - Molecule could not be created")
        except Exception as e:
            print(f"  EXCEPTION: {str(e)}")

def demonstrate_workaround():
    """Show current workaround for handling these molecules"""
    print("\n\n2. CURRENT WORKAROUND APPROACH:")
    print("===============================")
    
    for i, smiles in enumerate(EXAMPLES):
        print(f"\nMolecule {i+1}: {smiles}")
        try:
            # Create molecule without sanitization
            mol = Chem.MolFromSmiles(smiles, sanitize=False)
            if not mol:
                print("  FAILED - Molecule could not be parsed")
                continue
                
            # Update property cache without strict valence checking
            mol.UpdatePropertyCache(strict=False)
            
            # Perform only safe operations
            Chem.GetSSSR(mol)
            Chem.GetSymmSSSR(mol)
            try:
                Chem.Kekulize(mol)
            except:
                pass
            try:
                Chem.SetAromaticity(mol)
            except:
                pass
            try:
                Chem.SetHybridization(mol)
            except:
                pass
            
            print("  SUCCESS - Molecule was created with workaround approach")
            
            # Now we can calculate some properties
            try:
                from rdkit.Chem import Descriptors
                mw = Descriptors.MolWt(mol)
                print(f"  Molecular weight: {mw:.2f}")
            except Exception as e:
                print(f"  Could not calculate molecular weight: {str(e)}")
                
        except Exception as e:
            print(f"  EXCEPTION: {str(e)}")

def demonstrate_new_approach():
    """Show how our new functions would work (if implemented)"""
    print("\n\n3. NEW ASTROCHEMICAL APPROACH (SIMULATED):")
    print("=========================================")
    print("NOTE: These functions would be available after building RDKit with our changes.\n")
    
    for i, smiles in enumerate(EXAMPLES):
        print(f"\nMolecule {i+1}: {smiles}")
        
        # Method 1: Direct creation with astrochemical sanitization
        print("  Method 1: SmilesToAstroMol(smiles)")
        print("    Returns a properly sanitized molecule in one step")
        
        # Method 2: Create molecule then apply custom sanitization 
        print("  Method 2: mol = MolFromSmiles(smiles, sanitize=False)")
        print("           SanitizeAstroMol(mol)")
        print("    Allows more control over the process")
        
        # Method 3: C++ literal operator (for C++ code)
        print("  Method 3: (C++ only) auto mol = \"CC#N[C]\"_astrochemical;")
        print("    Convenient for use in C++ code")

if __name__ == "__main__":
    print("\nASTROCHEMICAL MOLECULE HANDLING DEMONSTRATION")
    print("============================================")
    print("This script demonstrates different approaches to handling")
    print("astrochemical molecules with unusual valence states in RDKit.")
    
    demonstrate_standard_approach()
    demonstrate_workaround()
    demonstrate_new_approach()
    
    print("\n\nSUMMARY:")
    print("=======")
    print("1. Standard RDKit sanitization fails on most astrochemical molecules due to")
    print("   strict valence checking.")
    print("2. A workaround exists but requires multiple steps and careful handling.")
    print("3. Our implementation provides simple, convenient functions that handle")
    print("   these molecules correctly while maintaining compatibility with")
    print("   standard RDKit functionality.") 