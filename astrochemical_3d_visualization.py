#!/usr/bin/env python
"""
Demonstration of generating and visualizing 3D structures for astrochemical molecules 
using the modified RDKit with astrochemical sanitization features.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import os
import numpy as np
import py3Dmol
from IPython.display import display
import matplotlib.pyplot as plt

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
    
    # Symmetrize rings
    try:
        Chem.GetSymmSSSR(mol)
    except Exception as e:
        print(f"  Warning: GetSymmSSSR failed - {str(e)}")
    
    # Set aromaticity
    try:
        Chem.SetAromaticity(mol)
    except Exception as e:
        print(f"  Warning: SetAromaticity failed - {str(e)}")
    
    # Set conjugation
    try:
        Chem.SetConjugation(mol)
    except Exception as e:
        print(f"  Warning: SetConjugation failed - {str(e)}")
    
    # Set hybridization
    try:
        Chem.SetHybridization(mol)
    except Exception as e:
        print(f"  Warning: SetHybridization failed - {str(e)}")
    
    # Update property cache again
    mol.UpdatePropertyCache(strict=False)
    
    return mol

def create_astro_mol(smiles):
    """
    Create a molecule using our astrochemical sanitization approach
    """
    if hasattr(Chem, 'MolFromSmilesAstro'):
        return Chem.MolFromSmilesAstro(smiles)
    else:
        # Simulate the function if not available
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return None
        return simulate_astro_sanitize(mol)

def generate_3d_structure(mol, name, n_conformers=10):
    """
    Generate 3D conformers and optimize the lowest energy one
    Returns the optimized molecule
    """
    if mol is None:
        return None
    
    try:
        # Make a copy to avoid modifying the original
        mol_3d = Chem.Mol(mol)
        
        # Generate multiple conformers
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        params.numThreads = 0  # Use all available CPUs
        
        print(f"Generating {n_conformers} conformers for {name}...")
        cids = AllChem.EmbedMultipleConfs(mol_3d, numConfs=n_conformers, params=params)
        
        if len(cids) == 0:
            print(f"  ✗ Failed to generate conformers for {name}")
            return None
        
        print(f"  ✓ Generated {len(cids)} conformers")
        
        # Calculate energy for each conformer and optimize them
        print("  Optimizing conformers with UFF...")
        energies = []
        for cid in cids:
            # Optimize the conformer
            try:
                energy = AllChem.UFFOptimizeMolecule(mol_3d, confId=cid)
                energies.append((cid, energy))
                print(f"    Conformer {cid}: Energy = {energy}")
            except Exception as e:
                print(f"    Failed to optimize conformer {cid}: {e}")
        
        if not energies:
            print(f"  ✗ Failed to optimize any conformers for {name}")
            return None
        
        # Find the lowest energy conformer
        energies.sort(key=lambda x: x[1])
        best_cid = energies[0][0]
        print(f"  ✓ Best conformer: {best_cid} with energy {energies[0][1]}")
        
        # Create a molecule with just the best conformer
        best_mol = Chem.Mol(mol_3d)
        best_conf = mol_3d.GetConformer(best_cid)
        best_mol.RemoveAllConformers()
        best_mol.AddConformer(best_conf)
        
        return best_mol
    
    except Exception as e:
        print(f"  ✗ Error generating 3D structure for {name}: {e}")
        return None

def save_molecule_formats(mol, name):
    """
    Save the molecule in various formats useful for visualization
    """
    if mol is None:
        return
    
    base_name = name.lower().replace(' ', '_')
    
    # Save as MOL file
    mol_file = f"{base_name}.mol"
    Chem.MolToMolFile(mol, mol_file)
    print(f"  ✓ Saved MOL file: {mol_file}")
    
    # Save as PDB file
    pdb_file = f"{base_name}.pdb"
    Chem.MolToPDBFile(mol, pdb_file)
    print(f"  ✓ Saved PDB file: {pdb_file}")
    
    # Save as SDF file
    sdf_file = f"{base_name}.sdf"
    writer = Chem.SDWriter(sdf_file)
    writer.write(mol)
    writer.close()
    print(f"  ✓ Saved SDF file: {sdf_file}")

def create_py3dmol_view(mol, width=700, height=400):
    """
    Create a py3Dmol visualization for the molecule
    """
    if mol is None:
        return None
    
    # Convert molecule to MOL block
    mol_block = Chem.MolToMolBlock(mol)
    
    # Set up py3Dmol viewer
    view = py3Dmol.view(width=width, height=height)
    view.addModel(mol_block, 'mol')
    
    # Style the molecule
    view.setStyle({'stick': {'radius': 0.15, 'colorscheme': 'cyanCarbon'}})
    view.addStyle({'elem': 'C'}, {'sphere': {'scale': 0.25, 'colorscheme': 'cyanCarbon'}})
    view.addStyle({'elem': 'H'}, {'sphere': {'scale': 0.2, 'color': 'white'}})
    view.addStyle({'elem': 'O'}, {'sphere': {'scale': 0.25, 'color': 'red'}})
    view.addStyle({'elem': 'N'}, {'sphere': {'scale': 0.25, 'color': 'blue'}})
    
    # Add dashed lines for hydrogen bonds
    view.addUnitCell()
    
    # Zoom to fit the molecule
    view.zoomTo()
    
    return view

def save_py3dmol_html(mol, name):
    """
    Save a py3Dmol visualization as an HTML file
    """
    if mol is None:
        return
    
    base_name = name.lower().replace(' ', '_')
    html_file = f"{base_name}_3d.html"
    
    # Convert molecule to MOL block
    mol_block = Chem.MolToMolBlock(mol)
    
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>3D Structure of {name}</title>
        <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
        <script src="https://3Dmol.org/build/3Dmol.ui-min.js"></script>
        <style>
            .mol-container {{
                width: 100%;
                height: 500px;
                position: relative;
            }}
            body {{
                font-family: Arial, sans-serif;
                margin: 20px;
            }}
            h2 {{
                color: #2c3e50;
            }}
        </style>
    </head>
    <body>
        <h2>3D Structure of {name}</h2>
        <div id="container" class="mol-container"></div>
        <script>
            let element = $('#container');
            let config = {{ backgroundColor: 'white' }};
            let viewer = $3Dmol.createViewer(element, config);
            let molData = `{mol_block}`;
            viewer.addModel(molData, "mol");
            viewer.setStyle({{}}, {{stick: {{radius: 0.15, colorscheme: 'cyanCarbon'}}}});
            viewer.addStyle({{elem: 'C'}}, {{sphere: {{scale: 0.25, colorscheme: 'cyanCarbon'}}}});
            viewer.addStyle({{elem: 'H'}}, {{sphere: {{scale: 0.2, color: 'white'}}}});
            viewer.addStyle({{elem: 'O'}}, {{sphere: {{scale: 0.25, color: 'red'}}}});
            viewer.addStyle({{elem: 'N'}}, {{sphere: {{scale: 0.25, color: 'blue'}}}});
            viewer.zoomTo();
            viewer.render();
        </script>
    </body>
    </html>
    """
    
    with open(html_file, 'w') as f:
        f.write(html_content)
    
    print(f"  ✓ Saved 3D visualization HTML: {html_file}")

def analyze_molecule_geometry(mol, name):
    """
    Analyze geometric properties of the 3D molecule
    """
    if mol is None or mol.GetNumConformers() == 0:
        return
    
    print(f"\nGeometric analysis for {name}:")
    
    # Get the conformer
    conf = mol.GetConformer()
    
    # Calculate bond lengths
    print("  Bond lengths:")
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtomIdx()
        atom2 = bond.GetEndAtomIdx()
        atom1_name = mol.GetAtomWithIdx(atom1).GetSymbol()
        atom2_name = mol.GetAtomWithIdx(atom2).GetSymbol()
        
        # Calculate the distance
        pos1 = conf.GetAtomPosition(atom1)
        pos2 = conf.GetAtomPosition(atom2)
        distance = np.sqrt((pos1.x - pos2.x)**2 + 
                           (pos1.y - pos2.y)**2 + 
                           (pos1.z - pos2.z)**2)
        
        bond_type = str(bond.GetBondType())
        print(f"    {atom1_name}{atom1}-{atom2_name}{atom2} ({bond_type}): {distance:.3f} Å")
    
    # Calculate angles for each atom that has at least 2 connections
    print("\n  Bond angles:")
    for atom_idx in range(mol.GetNumAtoms()):
        atom = mol.GetAtomWithIdx(atom_idx)
        neighbors = [bond.GetOtherAtomIdx(atom_idx) for bond in mol.GetBonds() if bond.GetBeginAtomIdx() == atom_idx or bond.GetEndAtomIdx() == atom_idx]
        
        if len(neighbors) >= 2:
            atom_name = atom.GetSymbol()
            
            # Calculate all angles for this atom
            for i in range(len(neighbors)):
                for j in range(i+1, len(neighbors)):
                    # Get positions
                    pos_center = conf.GetAtomPosition(atom_idx)
                    pos_i = conf.GetAtomPosition(neighbors[i])
                    pos_j = conf.GetAtomPosition(neighbors[j])
                    
                    # Calculate vectors
                    v1 = np.array([pos_i.x - pos_center.x, pos_i.y - pos_center.y, pos_i.z - pos_center.z])
                    v2 = np.array([pos_j.x - pos_center.x, pos_j.y - pos_center.y, pos_j.z - pos_center.z])
                    
                    # Normalize vectors
                    v1_norm = v1 / np.linalg.norm(v1)
                    v2_norm = v2 / np.linalg.norm(v2)
                    
                    # Calculate angle in degrees
                    dot_product = np.dot(v1_norm, v2_norm)
                    # Clip to prevent numerical errors
                    dot_product = max(min(dot_product, 1.0), -1.0)
                    angle = np.degrees(np.arccos(dot_product))
                    
                    neighbor1_name = mol.GetAtomWithIdx(neighbors[i]).GetSymbol()
                    neighbor2_name = mol.GetAtomWithIdx(neighbors[j]).GetSymbol()
                    
                    print(f"    {neighbor1_name}{neighbors[i]}-{atom_name}{atom_idx}-{neighbor2_name}{neighbors[j]}: {angle:.1f}°")

def define_astrochemical_molecules():
    """
    Define a set of astrochemical molecules to visualize in 3D
    """
    molecules = [
        # Linear molecules
        ("Polyyne anion", "C#CC#CC#CC#CC#[C-]"),
        ("Acetylene", "C#C"),
        ("Diacetylene", "C#CC#C"),
        
        # Simple molecules
        ("Carbon monoxide", "[C]=O"),
        ("Formaldehyde", "C=O"),
        
        # Radical species
        ("Cyanomethylidyne radical", "CC#N[C]"),
        ("Methyl radical", "[CH3]"),
        
        # Ions
        ("HCO+", "[CH+]=O"),
        ("N2H+", "[NH+]#N"),
        
        # Astrochemically relevant molecules
        ("Methanol", "CO"),
        ("Hydrogen cyanide", "C#N"),
        ("Ammonia", "N"),
        
        # PAHs
        ("Benzene", "c1ccccc1"),
        ("Naphthalene", "c1ccc2ccccc2c1")
    ]
    return molecules

def main():
    """
    Main function to demonstrate 3D structure generation for astrochemical molecules
    """
    print("="*80)
    print("RDKIT ASTROCHEMICAL 3D STRUCTURE VISUALIZATION")
    print("="*80)
    
    # Create a directory for the output files
    output_dir = "astrochem_3d_structures"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")
    
    os.chdir(output_dir)
    
    # Load molecules and generate 3D structures
    molecules = define_astrochemical_molecules()
    
    for name, smiles in molecules:
        print("\n" + "="*60)
        print(f"PROCESSING: {name} ({smiles})")
        print("="*60)
        
        # Create molecule using astrochemical approach
        mol = create_astro_mol(smiles)
        
        if mol is not None:
            print(f"  ✓ Created molecule successfully")
            print(f"    SMILES: {Chem.MolToSmiles(mol)}")
            
            # Generate 2D coordinates and save an image
            AllChem.Compute2DCoords(mol)
            img_file = f"{name.lower().replace(' ', '_')}_2d.png"
            img = Draw.MolToImage(mol, size=(300, 200))
            img.save(img_file)
            print(f"  ✓ Saved 2D image to {img_file}")
            
            # Generate 3D structure
            mol_3d = generate_3d_structure(mol, name)
            
            if mol_3d is not None:
                # Save molecule in various formats
                save_molecule_formats(mol_3d, name)
                
                # Create HTML visualization
                save_py3dmol_html(mol_3d, name)
                
                # Analyze geometric properties
                analyze_molecule_geometry(mol_3d, name)
            
        else:
            print(f"  ✗ Failed to create molecule")
    
    print("\n" + "="*80)
    print(f"3D structure generation completed!")
    print(f"Output files saved in: {os.path.abspath(output_dir)}")
    print("="*80)

if __name__ == "__main__":
    main() 