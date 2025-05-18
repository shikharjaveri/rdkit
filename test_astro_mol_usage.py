#!/usr/bin/env python
"""
Example usage of RDKit's astrochemical molecule functionality.

This script demonstrates how to use the new astrochemical molecule functions
in RDKit to handle molecules with unusual valence states that would normally
fail standard sanitization.

NOTE: This script requires an RDKit installation that includes the astrochemical
molecule sanitization functionality.
"""

# This example assumes RDKit is properly installed with our modifications
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

def demonstrate_astro_sanitization():
    """Demonstrate the astrochemical molecule sanitization features"""
    
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
    
    print("COMPARING STANDARD VS. ASTROCHEMICAL SANITIZATION\n")
    
    mols = []
    mol_labels = []
    
    for idx, smiles in enumerate(smiles_list):
        print(f"Molecule {idx+1}: {smiles}")
        
        # Try standard sanitization
        try:
            std_mol = Chem.MolFromSmiles(smiles)
            if std_mol:
                print("  Standard sanitization: SUCCESS")
            else:
                print("  Standard sanitization: FAILED")
        except Exception as e:
            print(f"  Standard sanitization: FAILED - {str(e)}")
        
        # Try astrochemical sanitization method 1:
        # Create unsanitized molecule and use custom sanitization
        try:
            mol1 = Chem.MolFromSmiles(smiles, sanitize=False)
            if mol1:
                # This would use the C++ implementation
                Chem.SanitizeAstroMol(mol1)
                print("  Method 1 (SanitizeAstroMol): SUCCESS")
                
                # Calculate molecular descriptors
                mw = Chem.Descriptors.MolWt(mol1)
                print(f"    Molecular weight: {mw:.2f}")
                
                # Generate 3D coordinates 
                AllChem.EmbedMolecule(mol1)
                print("    Generated 3D structure")
                
                mols.append(mol1)
                mol_labels.append(f"{idx+1}: {smiles}")
            else:
                print("  Method 1 (SanitizeAstroMol): FAILED - Could not parse SMILES")
        except AttributeError:
            print("  Method 1 (SanitizeAstroMol): FUNCTION NOT AVAILABLE (implementation needed)")
        except Exception as e:
            print(f"  Method 1 (SanitizeAstroMol): FAILED - {str(e)}")
        
        # Try astrochemical sanitization method 2:
        # Direct creation via SmilesToAstroMol
        try:
            mol2 = Chem.SmilesToAstroMol(smiles)
            if mol2:
                print("  Method 2 (SmilesToAstroMol): SUCCESS")
            else:
                print("  Method 2 (SmilesToAstroMol): FAILED - Returned None")
        except AttributeError:
            print("  Method 2 (SmilesToAstroMol): FUNCTION NOT AVAILABLE (implementation needed)")
        except Exception as e:
            print(f"  Method 2 (SmilesToAstroMol): FAILED - {str(e)}")
            
        # Try astrochemical sanitization method 3: 
        # String literal operator (from C++)
        try:
            # This would be used in C++ code:
            # auto mol = "CC#N[C]"_astrochemical;
            print("  Method 3 (_astrochemical literal): Available in C++ only")
        except Exception as e:
            print(f"  Method 3 (_astrochemical literal): FAILED - {str(e)}")
            
        print("")
    
    # If we have molecules, draw them as a grid
    if mols:
        try:
            img = Draw.MolsToGridImage(mols, molsPerRow=3, subImgSize=(200, 200), 
                                      legends=mol_labels)
            img.save('astrochemical_molecules.png')
            print("Saved molecule grid image to 'astrochemical_molecules.png'")
        except Exception as e:
            print(f"Could not create molecule grid image: {str(e)}")

if __name__ == "__main__":
    demonstrate_astro_sanitization()
    print("\nNOTE: This script demonstrates the expected API for astrochemical molecule")
    print("functionality. If you see 'FUNCTION NOT AVAILABLE' messages, it means the")
    print("implementation has not been compiled and installed in your RDKit version yet.") 