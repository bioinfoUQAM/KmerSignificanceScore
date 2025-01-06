import numpy as np
from sklearn.svm import SVC
from sklearn.metrics import f1_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict

def compute_rms_rank_weighted_f1_score(y_true, y_pred, classes, verbose=False):
    """
    Compute the Root Mean Square (RMS) Rank Weighted F1 Score.

    Parameters:
    - y_true (array-like): True labels.
    - y_pred (array-like): Predicted labels.
    - classes (array-like): Unique class labels.
    - verbose (bool): If True, prints detailed computation info.

    Returns:
    - float: The computed RMS Rank Weighted F1 Score.
    """
    y_true = np.array(y_true)
    y_pred = np.array(y_pred)
    
    # Calculate F1 scores for each class
    f1_scores_dict = {c: f1_score(y_true == c, y_pred == c) for c in classes}

    # Sort F1 scores in descending order
    sorted_scores = sorted(f1_scores_dict.items(), key=lambda item: item[1], reverse=True)
    f1_scores_sorted = dict(sorted_scores)

    # Calculate weights based on ranks
    n = len(f1_scores_sorted)
    ranks = np.arange(1, n + 1)
    weights = n / ranks ** (np.log(n) / 2)

    # Optional verbose output
    if verbose:
        for i, (class_id, score) in enumerate(f1_scores_sorted.items()):
            print(f"Class: {class_id}, Weight: {weights[i]:.2f}, F1-score: {score:.3f}")
        print()

    # Compute the RMS Rank Weighted F1 Score
    f1_scores_array = np.array(list(f1_scores_sorted.values()))
    numerator = np.sum(weights * f1_scores_array ** 2)
    denominator = np.sum(weights)
    rms_rank_weighted_f1_score = np.sqrt(numerator / denominator)

    return rms_rank_weighted_f1_score

def evaluate_model(estimator, X, y, n_splits=5, n_jobs=1):
    """
    Evaluate the model using cross-validation and return predictions.

    Parameters:
    - estimator: The classifier to use.
    - X (array-like): Feature matrix.
    - y (array-like): Target labels.
    - n_splits (int): Number of cross-validation folds.
    - n_jobs (int): Number of jobs to run in parallel.

    Returns:
    - array-like: Cross-validated predictions.
    """
    cv = StratifiedKFold(n_splits=n_splits, random_state=0, shuffle=True)
    y_pred = cross_val_predict(estimator, X, y, cv=cv, n_jobs=n_jobs)
    return y_pred

def categorize_scores(scores, thresholds):
    """
    Categorize scores based on provided thresholds.

    Parameters:
    - scores (list or array-like): Scores to categorize.
    - thresholds (list or array-like): Threshold values for categorization.

    Returns:
    - list: Categories corresponding to each score.
    """
    categories = [np.digitize(score, thresholds) + 1 for score in scores]
    return categories

def get_discriminative_scores(X, y):
    """
    Compute discriminative scores for the dataset.

    Parameters:
    - X (array-like): Feature matrix.
    - y (array-like): Target labels.

    Returns:
    - tuple: (categorized_scores, uncategorized_scores)
        - categorized_scores (list): Scores categorized into bins.
        - uncategorized_scores (list): Raw score values.
    """
    y = np.array(y)
    y_pred = evaluate_model(SVC(kernel="linear"), X, y, n_splits=5, n_jobs=-1)
    n_classes = len(np.unique(y))
    
    # Define thresholds based on the number of classes
    if n_classes <= 4:
        thresholds = [0.60, 0.70, 0.80, 0.90]
    elif n_classes <= 7:
        thresholds = [0.575, 0.675, 0.775, 0.875]
    else:
        thresholds = [0.55, 0.65, 0.75, 0.85]


    # Compute various F1 scores
    classes = np.unique(y)
    rms_rank_weighted_f1 = compute_rms_rank_weighted_f1_score(y, y_pred, classes)
    f1_macro = f1_score(y, y_pred, average='macro')
    f1_micro = f1_score(y, y_pred, average='micro')
    f1_weighted = f1_score(y, y_pred, average='weighted')

    f1_scores = [f1_macro, f1_micro, f1_weighted, rms_rank_weighted_f1]
    
    # Categorize the scores based on thresholds
    categorized_scores = categorize_scores(f1_scores, thresholds)
    uncategorized_scores = f1_scores

    return categorized_scores, uncategorized_scores
