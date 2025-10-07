# Scientific Data Enrichment Tool for Open WebUI

An Open WebUI tool that automatically enriches chemistry and materials science queries with contextual data from multiple scientific databases.

## Features

- **PubChem Integration**: Chemical structure data, molecular properties, and descriptions
- **ChEMBL Integration**: Drug bioactivity, clinical trial phases, therapeutic targets, and indications
- **Materials Project Integration**: Computational materials properties (band gap, density, formation energy, crystal structure)
- **RDKit Calculations**: Molecular descriptors, Lipinski rule-of-five compliance, drug-likeness metrics

## How It Works

The tool automatically detects and enriches:

1. **Chemical Formulas** (e.g., H2O, H2SO4, LiFePO4)
   - Looks up materials properties from Materials Project
   - Returns band gap, density, formation energy, space group, etc.

2. **Drug/Chemical Names** (e.g., ibuprofen, aspirin, caffeine)
   - Resolves name to structure via PubChem
   - Calculates molecular descriptors with RDKit
   - Retrieves drug information from ChEMBL

3. **Materials Project IDs** (e.g., mp-149)
   - Direct lookup of computational materials data

## Installation

### Option 1: Docker (Recommended)

1. Clone this repository:
```bash
git clone https://github.com/yourusername/scientific-enrichment-tool.git
cd scientific-enrichment-tool
```

2. Create a `.env` file with your Materials Project API key (optional but recommended):
```bash
echo "MATERIALS_PROJECT_API_KEY=your_api_key_here" > .env
```

Get a free API key at: https://next-gen.materialsproject.org/api

3. Start the service:
```bash
docker-compose up -d
```

The API will be available at `http://localhost:8099`

### Option 2: Local Python Installation

1. Install dependencies:
```bash
pip install -r requirements.txt
```

2. Set your Materials Project API key:
```bash
export MATERIALS_PROJECT_API_KEY=your_api_key_here
```

3. Run the API service:
```bash
python enrichment_api.py
```

## Open WebUI Integration

### Installing the Function

1. In Open WebUI, navigate to **Workspace** → **Functions**
2. Click **+ Add Function**
3. Upload the `scientific_enrichment_function.py` file
4. Configure the valves:
   - `ENRICHMENT_API_URL`: Set to `http://localhost:8099` (or your Docker service URL)
   - `MATERIALS_PROJECT_API_KEY`: Your Materials Project API key (optional)

### Using the Tool

1. Start a new chat in Open WebUI
2. Click the **tools icon** next to the message input
3. Enable **Scientific Data Enrichment**
4. Ask chemistry or materials science questions!

**Example queries:**
- "What happens if I mix ibuprofen with H2SO4?"
- "What are the properties of LiFePO4?"
- "Tell me about aspirin's mechanism of action"
- "What is the band gap of mp-149?"

The tool will automatically detect relevant chemicals and materials, fetch data from scientific databases, and provide enriched context to the AI for more informed responses.

## API Endpoints

### Chemistry Tools

#### POST `/tools/pubchem/lookup`
Resolve chemical name or formula via PubChem.

**Request:**
```json
{"name": "ibuprofen"}
```
or
```json
{"formula": "C13H18O2"}
```

**Response:**
```json
{
  "success": true,
  "cid": 3672,
  "description": "Ibuprofen is a nonsteroidal anti-inflammatory drug...",
  "properties": {
    "MolecularFormula": "C13H18O2",
    "MolecularWeight": "206.28",
    "IUPACName": "2-[4-(2-methylpropyl)phenyl]propanoic acid"
  },
  "smiles": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
  "iupac": "2-[4-(2-methylpropyl)phenyl]propanoic acid"
}
```

#### POST `/tools/chembl/lookup`
Lookup bioactivity, targets, and drug information from ChEMBL.

**Request:**
```json
{"name": "aspirin"}
```

**Response:**
```json
{
  "success": true,
  "chembl_id": "CHEMBL25",
  "pref_name": "ASPIRIN",
  "max_phase": 4,
  "properties": {
    "molecular_weight": 180.16,
    "alogp": 1.19
  },
  "indications": [
    {"mesh_heading": "Pain", "max_phase": 4}
  ],
  "bioactivity_targets": [
    {"target_id": "CHEMBL233", "target_type": "SINGLE PROTEIN"}
  ]
}
```

#### POST `/tools/chem/lipinski_pains`
Calculate molecular descriptors and Lipinski rule-of-five compliance.

**Request:**
```json
{"smiles": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"}
```

**Response:**
```json
{
  "success": true,
  "lipinski": {
    "mw": 206.28,
    "hbd": 1,
    "hba": 2,
    "logp": 3.5,
    "tpsa": 37.3,
    "passes": true
  }
}
```

### Materials Tools

#### POST `/tools/materials/lookup`
Lookup materials by formula, Materials Project ID, or elements.

**Request:**
```json
{"formula": "LiFePO4"}
```
or
```json
{"mp_id": "mp-19017"}
```
or
```json
{"elements": ["Li", "Fe", "O"]}
```

**Response:**
```json
{
  "success": true,
  "data": [
    {
      "material_id": "mp-19017",
      "formula_pretty": "LiFePO4",
      "band_gap": 3.41,
      "density": 3.6,
      "formation_energy_per_atom": -2.197,
      "e_above_hull": 0.0,
      "is_metal": false,
      "spacegroup_symbol": "Pnma"
    }
  ]
}
```

## Example Enriched Context

When you ask: **"What happens if I mix ibuprofen with H2SO4?"**

The tool provides:
```
Scientific Context:
- Materials Project (H2SO4) | ID: mp-558129 | Band gap: 0.00 eV | Density: 1.83 g/cm³ | Metal: No
- PubChem (ibuprofen) | IUPAC: 2-[4-(2-methylpropyl)phenyl]propanoic... | Formula: C13H18O2 | MW: 206.28
- RDKit (ibuprofen) | MW: 206.3 | HBD: 1 | HBA: 2 | LogP: 3.5 | TPSA: 37.3 | Lipinski: ✓
- ChEMBL (ibuprofen) | Drug: IBUPROFEN | Status: Approved | Treats: Pain, Inflammation

User Query: What happens if I mix ibuprofen with H2SO4?
```

This enriched context helps the AI provide accurate, informed responses about chemical reactions, drug properties, and safety considerations.

## Data Sources

- **PubChem**: Free public database of chemical structures and properties (NIH)
- **ChEMBL**: Free database of bioactive drug-like molecules (EMBL-EBI)
- **Materials Project**: Computational materials science database (requires free API key)
- **RDKit**: Open-source cheminformatics library

## Requirements

- Python 3.11+
- Docker (recommended)
- Materials Project API key (optional but recommended for best results)

## License

MIT License - feel free to use and modify for your projects.

## Contributing

Contributions welcome! Please open an issue or pull request.

## Related Projects

- [AB-MCTS Reasoning Engine](https://github.com/johnsonfarmsus/openwebui-ab-mcts-pipeline) - Advanced tree search reasoning for Open WebUI
