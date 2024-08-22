# coding: utf-8

import pandas as pd
import re
import numpy as np

# File URL
url = '/scratch/Funcoup5/data_used/10090.mitab.01-22-2018.txt'

def iRefIndex_parser(file_path):
    """
    Parses the iRefIndex dataset to extract interaction types, methods, tax IDs,
    and handles complex interactions.

    Args:
        file_path (str): Path to the iRefIndex file.

    Returns:
        None
    """
    # Load the dataset
    data = pd.read_csv(file_path, sep='\t')
    data = data.rename(columns={'#uidA': 'uidA'})
    
    # Select relevant columns
    filtered_data = data[['uidA', 'uidB', 'interactionType', 'method', 'Host_organism_taxid']]

    # Extract interaction type
    filtered_data['interactionType'] = filtered_data['interactionType'].apply(
        lambda x: re.findall(r"\((.*)\)", x)
    )
    filtered_data['interactionType'] = filtered_data['interactionType'].apply(
        lambda x: x[0] if x else 'No_info'
    )

    # Extract experiment method
    filtered_data['method'] = filtered_data['method'].apply(
        lambda x: re.findall(r"\((.*)\)", x)
    )
    filtered_data['method'] = filtered_data['method'].apply(
        lambda x: x[0] if x else 'No_info'
    )

    # Extract tax IDs
    filtered_data['Host_organism_taxid'] = filtered_data['Host_organism_taxid'].apply(
        lambda x: re.findall(r":(.*)\(", x)
    )
    filtered_data['Host_organism_taxid'] = filtered_data['Host_organism_taxid'].apply(
        lambda x: int(x[0]) if x and int(x[0]) > 0 else 0
    )

    # Handle complex interactions
    complex_data = filtered_data[filtered_data['uidA'].str.contains("complex")]
    complex_dict = {
        key: list(set(value)) 
        for key, value in complex_data.groupby('uidA')['uidB'].apply(list).to_dict().items()
    }
    
    complex_df = pd.DataFrame()
    for key, proteins in complex_dict.items():
        protein_pairs = list(itertools.combinations(proteins, 2))
        temp_df = pd.DataFrame(protein_pairs, columns=['uidA', 'uidB'])
        temp_values = complex_data[complex_data['uidA'] == key].iloc[1][['interactionType', 'method', 'Host_organism_taxid']]
        repeated_attributes = pd.DataFrame([temp_values] * len(temp_df))
        result = pd.concat([temp_df.reset_index(drop=True), repeated_attributes.reset_index(drop=True)], axis=1)
        complex_df = pd.concat([complex_df, result], ignore_index=True)

    # Mark complex interactions
    complex_df['Gold_Standard'] = 'Complex'

    # Handle non-complex interactions
    no_complex_data = filtered_data[~filtered_data['uidA'].str.contains('complex')]
    no_complex_data['Gold_Standard'] = 'PPI'

    # Combine complex and non-complex interactions
    final_dataset = pd.concat([no_complex_data, complex_df], ignore_index=True)

    # Clean up the uidA and uidB columns
    final_dataset['uidA'] = final_dataset['uidA'].str.split(':').str[1]
    final_dataset['uidB'] = final_dataset['uidB'].str.split(':').str[1]

    # Filter out unwanted rows
    to_drop = ['No_info', 'unspecified method', 0, '0']
    final_dataset = final_dataset[~final_dataset.isin(to_drop).any(axis=1)]

    # Save the cleaned dataset
    final_dataset.to_csv('/scratch/Funcoup5/Clean_iRefIndex2.csv', sep='\t', index=False)

# Run the parser function
iRefIndex_parser(url)
