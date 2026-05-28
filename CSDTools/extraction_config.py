ENTRY_PROPS = [
    "identifier",
    "is_organic",
    "has_disorder",
    "disorder_details",
    "color",
    "polymorph",
    "phase_transition",
    "radiation_source",
    "temperature",
    "calculated_density",
    "melting_point",
    "melting_point_default_units",
#    "heat_capacity",        # Solubility Platform not available
#    "heat_capacity_notes",
#    "heat_of_fusion",
#    "heat_of_fusion_notes",
#    "solubility_data",
    "bioactivity",
    "ccdc_number",
    "doi",
#   "publications",         # Citation tuple → extracted manually in chunk_extraction
    "solvent",
]

CRYSTAL_PROPS = [
#   "spacegroup_number_and_setting",  # (number, setting) tuple → split into 2 columns in chunk_extraction
    "spacegroup_symbol",
#   "symmetry_operators",             # str tuple → extracted manually in chunk_extraction
    "z_prime",
    "z_value",
    "cell_volume",
#   "cell_lengths",                   # non-native object → split into 3 columns in chunk_extraction
#   "cell_angles",                    # non-native object → split into 3 columns in chunk_extraction
]

MOLECULE_PROPS = [
    "smiles",
    "molecular_volume",
    "molecular_weight",
]

COLUMN_ORDER = [
    # ENTRY
    "identifier", "is_organic", "has_metal", "num_components",
    "has_disorder", "disorder_details", "color",
    "polymorph", "phase_transition", "radiation_source", "temperature",
    "calculated_density", "melting_point", "melting_point_default_units",
    "bioactivity", "ccdc_number", "doi", "publications", "publication_years", "solvent",
    # MOLECULE
    "smiles", "molecular_volume", "molecular_weight", "inchi",
    # CRYSTAL
    "spacegroup_number", "spacegroup_setting", "spacegroup_symbol",
    "symmetry_operators", "z_prime", "z_value",
    "cell_length_a", "cell_length_b", "cell_length_c",
    "cell_angle_alpha", "cell_angle_beta", "cell_angle_gamma",
    "cell_volume",
]
