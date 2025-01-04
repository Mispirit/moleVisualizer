from flask import Flask, request, jsonify, send_file
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
from io import BytesIO
from flask_cors import CORS

app = Flask(__name__)
CORS(app)

# Helper function to calculate chemical properties
def calculate_properties(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    properties = {
        "molecular_weight": Descriptors.MolWt(mol),
        "logP": Descriptors.MolLogP(mol),
        "TPSA": Descriptors.TPSA(mol),
        "num_h_acceptors": Descriptors.NumHAcceptors(mol),
        "num_h_donors": Descriptors.NumHDonors(mol),
        "num_rotatable_bonds": Descriptors.NumRotatableBonds(mol),
        "formula": Chem.rdMolDescriptors.CalcMolFormula(mol),
        # "name": get_compound_name(smiles),
    }
    return properties

# API route to calculate properties
@app.route('/api/calculate', methods=['POST'])
def calculate():
    data = request.json
    smiles = data.get('smiles')

    if not smiles:
        return jsonify({"error": "No SMILES string provided"}), 400

    # Calculate properties
    properties = calculate_properties(smiles)

    if properties is None:
        return jsonify({"error": "Invalid SMILES input"}), 400

    return jsonify(properties), 200

# API route to generate a 2D image of the molecule
@app.route('/api/render2d', methods=['POST'])
def render2d():
    data = request.json
    smiles = data.get('smiles')

    if not smiles:
        return jsonify({"error": "No SMILES string provided"}), 400

    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        return jsonify({"error": "Invalid SMILES input"}), 400

    # Generate 2D image of the molecule
    img = Draw.MolToImage(mol)

    # Convert image to binary stream
    img_io = BytesIO()
    img.save(img_io, 'PNG')
    img_io.seek(0)

    return send_file(img_io, mimetype='image/png')

# Run the app
if __name__ == "__main__":
    app.run(debug=True, port=5001)
