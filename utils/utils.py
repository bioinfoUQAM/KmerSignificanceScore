import json
import numpy as np
import pandas as pd

def load_data_from_json(filename):
    """
    Load data from a JSON file.

    Parameters:
        filename (str): Path to the JSON file.

    Returns:
        dict: Data loaded from the JSON file.

    Raises:
        FileNotFoundError: If the file does not exist.
        json.JSONDecodeError: If the file contains invalid JSON.
    """
    try:
        with open(filename, 'r') as file:
            return json.load(file)
    except json.JSONDecodeError as e:
        print(f"Error decoding JSON from {filename}: {e}")
        raise
    except FileNotFoundError as e:
        print(f"File not found: {filename}")
        raise

def save_data_as_json(data, filename):
    """
    Save data to a JSON file, converting non-serializable keys and values to
    appropriate JSON-compatible types.

    Parameters:
    - data: The data to save.
    - filename: The filename for the JSON file.
    """
    def convert_keys(obj):
        if isinstance(obj, dict):
            return {str(key): convert_keys(value) for key, value in obj.items()}
        elif isinstance(obj, list):
            return [convert_keys(element) for element in obj]
        elif isinstance(obj, set):
            return list(obj)  # Convert sets to lists for JSON compatibility
        elif isinstance(obj, np.ndarray):
            return obj.tolist()  # Convert NumPy arrays to lists
        elif isinstance(obj, (np.integer, np.floating)):
            return obj.item()  # Convert NumPy numbers to native Python types
        elif isinstance(obj, np.bool_):
            return bool(obj)  # Convert NumPy booleans to Python bool
        else:
            return obj  # Return the object if it's already serializable

    data_serializable = convert_keys(data)

    with open(filename, 'w') as file:
        json.dump(data_serializable, file, indent=4)

def load_annotations(csv_filepath: str) -> dict:
    """
    Load annotations from a CSV file and return them as a dictionary.

    This function reads a CSV file with two columns, where the first column is used as the keys
    and the second column as the values for the dictionary. The first column is set as the index,
    and the second column's values are converted to a dictionary.

    Parameters:
    -----------
    csv_filepath : str
        The file path to the CSV file to be loaded.

    Returns:
    --------
    dict
        A dictionary where keys are from the first column and values from the second column
        of the CSV file.
    
    Example:
    --------
    If the CSV file looks like this:
    
        id;annotation
        1;label1
        2;label2

    The output will be:
    
        {1: 'label1', 2: 'label2'}
    """
    df = pd.read_csv(csv_filepath, sep=';', skipinitialspace=True)
    return df.set_index(df.columns[0]).to_dict()[df.columns[1]]
