# Personal Chemoinformatics project by Chikelu Elizabeth #Clizdesigns
import os
from flask import Flask, request, jsonify, send_file
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
from io import BytesIO
from dotenv import load_dotenv
from flask_cors import CORS

# Load environment variables
load_dotenv()

# Get configuration from environment variables
FLASK_ENV = os.getenv("FLASK_ENV", "development")
FLASK_PORT = int(os.getenv("PORT", 5000))  # Render assigns this automatically

app = Flask(__name__)
CORS(app, origins=["https://elizabethmolevisualizer.netlify.app"])  # Allow frontend requests

# Function to calculate chemical properties
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
    }
    return properties

# Properties calculation route
@app.route('/api/calculate', methods=['POST'])
def calculate():
    data = request.json
    smiles = data.get('smiles')

    if not smiles:
        return jsonify({"error": "No SMILES string provided"}), 400

    properties = calculate_properties(smiles)
    if properties is None:
        return jsonify({"error": "Invalid SMILES input"}), 400

    return jsonify(properties), 200

# Route to generate a 2D image of the molecule
@app.route('/api/render2d', methods=['POST'])
def render2d():
    data = request.json
    smiles = data.get('smiles')

    if not smiles:
        return jsonify({"error": "No SMILES string provided"}), 400

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return jsonify({"error": "Invalid SMILES input"}), 400

    img = Draw.MolToImage(mol)

    img_io = BytesIO()
    img.save(img_io, 'PNG')
    img_io.seek(0)

    return send_file(img_io, mimetype='image/png')

# Start the app using environment variables
if __name__ == "__main__":
    app.run(debug=(FLASK_ENV == "development"), port=FLASK_PORT, host="0.0.0.0")
