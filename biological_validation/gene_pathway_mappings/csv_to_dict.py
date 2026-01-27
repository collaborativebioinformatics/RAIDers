import pandas as pd

file = "ALS Gene-Pathway Validation - Gene-to-pathway mapping table.csv"
df = pd.read_csv(file, header=0)

#print(df.head())
#print(df.shape) # (32, 5)

#print(df.columns)

df = df.rename(columns={'x': 'Variant_Gene',
                         'Primary Pathway': 'Primary_Pathway',
                         'Secondary Pathway(s)': 'Secondary_Pathways',
                         'Mechanism Summary': 'Mechanism_Summary'})

df_wo_citations = df.drop(columns=['Citation #s'])



def create_metadata_dict(df):
    metadata = {}

    for _, row in df.iterrows():
        # Variant Gene would be the key
        gene = row['Variant_Gene'].strip()

        # Parse secondary pathways - split by comma into list
        secondary_raw = row['Secondary_Pathways']
        if str(secondary_raw).lower() == 'none' or secondary_raw == '' or pd.isna(secondary_raw):
            # If 'none', empty, NaN, make it an empty list
            secondary = []
        else:
            # Split by comma and strip whitespace
            secondary = [s.strip() for s in str(secondary_raw).split(',')]

        metadata[gene] = {
            'primary': str(row['Primary_Pathway']).strip(),
            'secondary': secondary,
            'mechanism': row['Mechanism_Summary'].strip()
        }
    return metadata


def create_pathway_dict(metadata):
    pathway_dict = {
        gene: [data["primary"]] + data["secondary"]
        for gene, data in metadata.items()
    }
    return pathway_dict

def get_all_pathways(gene_pathway):
    all_pathways = set()
    for pathway_list in gene_pathway.values():
        for pathway in pathway_list:
            all_pathways.add(pathway)
    return sorted(all_pathways)


def write_to_file(metadata, gene_to_pathway, all_pathways):
    """Write dictionaries to als_master_list.py"""

    with open('als_master_list.py', 'w') as f:
        f.write('"""\n')
        f.write('ALS Gene-to-Pathway Mappings\n')
        f.write('Auto-generated from gene_pathway_raw.csv\n')
        f.write('Last updated: January 2026\n')
        f.write('"""\n\n')

        # Write GENE_PATHWAY_METADATA
        f.write('GENE_PATHWAY_METADATA = {\n')
        for gene, data in metadata.items():
            f.write(f'    "{gene}": {{\n')
            f.write(f'        "primary": "{data["primary"]}",\n')
            f.write(f'        "secondary": {data["secondary"]},\n')
            f.write(f'        "mechanism": "{data["mechanism"]}"\n')
            f.write(f'    }},\n')
        f.write('}\n\n')

        # Write GENE_TO_PATHWAY
        f.write('GENE_TO_PATHWAY = {\n')
        for gene, pathways in gene_to_pathway.items():
            f.write(f'    "{gene}": {pathways},\n')
        f.write('}\n\n')

        # Write ALL_PATHWAYS
        f.write(f'ALL_PATHWAYS = {all_pathways}\n')

metadata = create_metadata_dict(df)
gene_to_pathway = create_pathway_dict(metadata)
pathways = get_all_pathways(gene_to_pathway)
print(metadata)
print(gene_to_pathway)
print(pathways)

# Write
write_to_file(metadata, gene_to_pathway, pathways)

print(f"✓ Generated mappings for {len(metadata)} genes")
print(f"✓ Found {len(pathways)} unique pathways")
