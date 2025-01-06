import os
import json
import warnings
from io import StringIO

import joblib
import requests
import numpy as np
import pandas as pd
from bioservices import UniProt
from sklearn.exceptions import InconsistentVersionWarning

# Suppress InconsistentVersionWarning from scikit-learn
warnings.filterwarnings(action='ignore', category=InconsistentVersionWarning)

def fetch_uniprot_data(taxon_id, gene_name):
    """
    Fetch UniProt data for a given taxon ID and gene name.

    If data is not already cached locally, it will download and save the JSON file.

    Parameters
    ----------
    taxon_id : int or str
        The NCBI Taxonomy ID of the organism.
    gene_name : str
        The name of the gene to query.

    Returns
    -------
    str
        The UniProt entry ID if data is found.
        Returns "No data found for the specified query." if no data is found.
    """
    # Adjust gene name if necessary
    gene = "ORF1ab" if gene_name == "ORF1b" else gene_name

    # Define the query for UniProt
    query = f"{taxon_id} AND {gene} AND reviewed:true"
    uniprot = UniProt(verbose=False)
    data = uniprot.search(query=query, frmt="tsv", limit=5, columns="")
    df = pd.read_csv(StringIO(data), sep='\t')

    if not df.empty and 'Entry' in df.columns:
        entry_id = df["Entry"].iloc[0]
        current_dir = os.path.dirname(os.path.abspath(__file__))
        json_file_path = os.path.join(current_dir, 'uniprot', f'{entry_id}.json')

        # Download the JSON data if it does not exist
        if not os.path.exists(json_file_path):
            endpoint = f"https://www.uniprot.org/uniprot/{entry_id}.json"
            response = requests.get(endpoint)
            response.raise_for_status()
            os.makedirs(os.path.dirname(json_file_path), exist_ok=True)
            with open(json_file_path, "w") as outfile:
                outfile.write(response.text)
        return entry_id
    else:
        return "No data found for the specified query."

def count_entries_in_cross_references(data, database, properties_value=None):
    """
    Count the number of entries in UniProt cross-references for a specified database.

    Parameters
    ----------
    data : dict
        The UniProt entry data in JSON format.
    database : str
        The name of the database to count entries for (e.g., "GO", "InterPro").
    properties_value : str, optional
        Specific property value to match in the cross-reference entries.

    Returns
    -------
    int
        The count of entries matching the specified database and property value.
    """
    return sum(
        d["database"] == database and
        (properties_value is None or d["properties"][0]["value"][0] == properties_value)
        for d in data.get("uniProtKBCrossReferences", [])
    )

def extract_keywords(data):
    """
    Extract keyword counts and categories from UniProt data.

    Parameters
    ----------
    data : dict
        The UniProt entry data in JSON format.

    Returns
    -------
    tuple
        - keyword_category_counts : dict
            A dictionary with keyword categories as keys and their counts as values.
        - keyword_name_counts : dict
            A dictionary with keyword names as keys and True as values to indicate presence.
    """
    keyword_name_counts = {}
    keyword_category_counts = {}
    for d in data.get("keywords", []):
        name = d["name"].replace(" ", "_")
        category = d["category"].replace(" ", "_")
        keyword_name_counts[name] = True
        keyword_category_counts[category] = keyword_category_counts.get(category, 0) + 1
    return keyword_category_counts, keyword_name_counts

def extract_features_from_json(data):
    """
    Extract features from UniProt JSON data for use in the prediction model.

    Parameters
    ----------
    data : dict
        The UniProt entry data in JSON format.

    Returns
    -------
    dict
        A dictionary containing extracted features for model input.
    """
    protein = data.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "")
    x_data = {
        "is_polymerase": int("polymerase" in protein.lower()),
        "is_uncharacterized": int("uncharacterized" in protein.lower()),
        "is_probable": int("probable" in protein.lower() or "putative" in protein.lower()),
        "Annotation_score": data.get("annotationScore", 0),
        "Protein_existence": int(data.get("proteinExistence", [0])[0]),
        "Number_keywords": len(data.get("keywords", [])),
        "Number_references": len(data.get("references", [])),
        "uniProtKBCrossReferences": len(data.get("uniProtKBCrossReferences", []))
    }

    # Add specific counting functions here
    x_data["n_GO_molecular_functions"] = count_entries_in_cross_references(data, "GO", "F")
    x_data["n_GO_biological_processes"] = count_entries_in_cross_references(data, "GO", "P")
    x_data["n_GO_cellular_components"] = count_entries_in_cross_references(data, "GO", "C")

    # Extracting categories and keywords
    keyword_category_counts, keyword_name_counts = extract_keywords(data)
    x_data["Number_keywords_names"] = len(keyword_name_counts)
    x_data.update(keyword_name_counts)
    x_data["Number_keywords_categories"] = len(keyword_category_counts)
    x_data.update({f"Number_keywords_{category}": count for category, count in keyword_category_counts.items()})
    x_data["Number_family_and_domains"] = count_entries_in_cross_references(data, "InterPro")
    x_data["PDB_3D_structures"] = count_entries_in_cross_references(data, "PDB")
    x_data["SMR_3D_structures"] = count_entries_in_cross_references(data, "SMR")
    x_data["Number_biological_pathways"] = count_entries_in_cross_references(data, "Reactome")
    x_data["Number_drug_targets"] = count_entries_in_cross_references(data, "DrugBank")
    x_data["Number_protein_protein_interactions_BioGRID"] = count_entries_in_cross_references(data, "BioGRID")
    x_data["Number_protein_protein_interactions_IntAct"] = count_entries_in_cross_references(data, "IntAct")

    return x_data

def refine_predictions(y_pred, df):
    """
    Refine the predictions based on specific rules.

    - If 'Protein_existence' is 5, set the prediction to 1.
    - For predictions equal to 5, check if any of the GO terms are missing.
      If so, downgrade the prediction to 4.

    Parameters
    ----------
    y_pred : list or np.array
        List or array of predicted class labels.
    df : pd.DataFrame
        DataFrame containing feature data corresponding to the predictions.

    Returns
    -------
    np.array
        Array of refined predictions.
    """
    y_pred = np.array(y_pred)

    for i in range(len(y_pred)):
        if df.iloc[i]["Protein_existence"] == 5:
            y_pred[i] = 1

        if y_pred[i] == 5:
            mol_func_sum = df.iloc[i]["Number_keywords_Molecular_function"] + df.iloc[i]["n_GO_molecular_functions"]
            bio_proc_sum = df.iloc[i]["Number_keywords_Biological_process"] + df.iloc[i]["n_GO_biological_processes"]
            cell_comp_sum = df.iloc[i]["Number_keywords_Cellular_component"] + df.iloc[i]["n_GO_cellular_components"]

            if mol_func_sum == 0 or bio_proc_sum == 0 or cell_comp_sum == 0:
                y_pred[i] = 4

    return y_pred

def generate_data_from_json(annotations, json_folder_path):
    """
    Generate data for model prediction from UniProt JSON files.

    Parameters
    ----------
    annotations : dict
        Dictionary where keys are UniProt entry IDs and values are labels.
    json_folder_path : str
        Path to the folder containing UniProt JSON files.

    Returns
    -------
    tuple
        - X : np.ndarray
            Feature matrix for model prediction.
        - y : np.ndarray
            Array of labels.
        - ids : np.ndarray
            Array of UniProt entry IDs.
        - df : pd.DataFrame
            DataFrame containing the features.
    """
    X, y, ids = [], [], []
    for annotation, label in annotations.items():
        json_file_path = os.path.join(json_folder_path, f"{annotation}.json")

        # Verify the existence of the file before attempting to open it
        if not os.path.exists(json_file_path):
            raise FileNotFoundError(f"The file {json_file_path} does not exist.")

        with open(json_file_path, 'r') as f:
            data = json.load(f)

        x_data = extract_features_from_json(data)
        X.append(x_data)
        ids.append(annotation)
        y.append(int(label))

    df = pd.DataFrame(X).fillna(0)
    return df.to_numpy(), np.array(y), np.array(ids), df

def adjust_test_data(df_train, df_test):
    """
    Adjust the columns of the test DataFrame to match those of the training DataFrame.

    Any missing columns in the test set are filled with zeros.
    NaN values are also replaced with zeros.

    Parameters
    ----------
    df_train : pd.DataFrame
        The training DataFrame with the target column structure.
    df_test : pd.DataFrame
        The test DataFrame to be adjusted.

    Returns
    -------
    np.ndarray
        A numpy array of the adjusted test DataFrame.
    """
    # Initialize a DataFrame with the same columns as df_train
    adjusted_df_test = pd.DataFrame(columns=df_train.columns)

    # Match test DataFrame columns to training DataFrame
    for col in df_train.columns:
        if col in df_test.columns:
            adjusted_df_test[col] = df_test[col]
        else:
            adjusted_df_test[col] = 0

    # Fill any missing values with 0
    adjusted_df_test.fillna(0, inplace=True)

    # Convert the DataFrame to a numpy array
    X_test_adjusted = adjusted_df_test.to_numpy()

    return X_test_adjusted

def get_protein_score(taxon_ids, gene):
    """
    Get the protein score for a given taxon ID(s) and gene.

    Accepts either a single taxon ID or a list of taxon IDs.

    Parameters
    ----------
    taxon_ids : int or list of int
        A single taxon ID or a list of taxon IDs.
    gene : str
        The name of the gene.

    Returns
    -------
    int or list of int
        The protein score(s) corresponding to the given taxon ID(s).
        Returns a single score if a single taxon ID is provided,
        or a list of scores if multiple taxon IDs are provided.

    Notes
    -----
    If no data is found for a taxon ID and gene, a default score of 1 is assigned.
    """
    try:
        # Ensure taxon_ids is a list
        if not isinstance(taxon_ids, list):
            taxon_ids = [taxon_ids]

        # Initialize a list to store protein scores
        protein_scores = []

        # Load model components outside the loop for efficiency
        current_dir = os.path.dirname(os.path.abspath(__file__))
        df_ref = pd.read_csv(os.path.join(current_dir, 'ml', 'df_ref.csv'))
        scaler = joblib.load(os.path.join(current_dir, 'ml', 'scaler.pkl'))
        selector = joblib.load(os.path.join(current_dir, 'ml', 'selector.pkl'))
        classifier = joblib.load(os.path.join(current_dir, 'ml', 'classifier.pkl'))

        for taxon_id in taxon_ids:
            # Fetching UniProt data for the gene
            entry_id = fetch_uniprot_data(taxon_id, gene)
            if entry_id == "No data found for the specified query.":
                protein_scores.append(1)  # Default score when data is not found
                continue  # Skip to the next taxon_id

            # Processing test data
            X_test, y_test, ids_test, df_test = generate_data_from_json(
                {entry_id: 0},
                os.path.join(current_dir, 'uniprot')
            )
            X_test_adjusted = adjust_test_data(df_ref, df_test)
            X_test_scaled = scaler.transform(X_test_adjusted)
            X_test_selected = selector.transform(X_test_scaled)
            y_test_pred = classifier.predict(X_test_selected)
            protein_score = refine_predictions(y_test_pred, df_test)[0]
            protein_scores.append(protein_score)

        return protein_scores if len(protein_scores) > 1 else protein_scores[0]
    except Exception as e:
        print(f"Error in get_protein_score: {e}")
        return 1