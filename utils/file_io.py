import os
import numpy as np

def load_numpy_file(file_path):
    """
    Load a numpy file and return the data.
    
    Parameters:
    file_path (str): The path to the numpy file.
    
    Returns:
    np.ndarray: The loaded data.
    """
    if os.path.exists(file_path):
        data = np.load(file_path)
        return data
    else:
        raise FileNotFoundError(f"The file {file_path} does not exist.")
    
    