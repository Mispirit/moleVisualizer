import sqlite3

def create_database():
    conn = sqlite3.connect('compounds.db')
    cursor = conn.cursor()
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS compounds (
        id INTEGER PRIMARY KEY,
        smiles TEXT UNIQUE,
        formula TEXT,
        molecular_weight REAL,
        logP REAL,
        TPSA REAL,
        num_h_acceptors INTEGER,
        num_h_donors INTEGER,
        num_rotatable_bonds INTEGER
    )
    ''')
    conn.commit()
    conn.close()

def insert_compound(smiles, formula, molecular_weight, logP, TPSA, num_h_acceptors, num_h_donors, num_rotatable_bonds):
    conn = sqlite3.connect('compounds.db')
    cursor = conn.cursor()
    cursor.execute('''
    INSERT OR IGNORE INTO compounds (smiles, formula, molecular_weight, logP, TPSA, num_h_acceptors, num_h_donors, num_rotatable_bonds)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?)
    ''', (smiles, formula, molecular_weight, logP, TPSA, num_h_acceptors, num_h_donors, num_rotatable_bonds))
    conn.commit()
    conn.close()

# Call this function to create the database
create_database()
