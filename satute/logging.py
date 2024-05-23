def format_matrix(matrix, precision: int = 4):
    """Format a matrix for pretty printing."""
    formatted_matrix = "\n".join(
        ["\t".join([f"{item:.{precision}f}" for item in row]) for row in matrix]
    )
    return formatted_matrix

def format_array(array, precision=4):
    """Format a 1D array for pretty printing."""
    formatted_array = "\t".join([f"{item:.{precision}f}" for item in array])
    return formatted_array
