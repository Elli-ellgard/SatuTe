import re
from typing import List, Dict

class BaseParser:
    def __init__(self, file_content: List[str]):
        self.file_content = file_content
        
    def find_line_index(self,lines: List[str], search_string: str) -> int:
        """
        Returns the index of the first line that contains the given search string.
        Raises an exception if the search string is not found.

        Args:
        - lines (List[str]): The list of lines to search through.
        - search_string (str): The string to search for.

        Returns:
        - int: The index of the first line containing the search string.

        Raises:
        - ValueError: If the search string is not found in any line.
        """
        for idx, line in enumerate(lines):
            if search_string in line:
                return idx
        raise ValueError(
            f"Search string '{search_string}' not found in the provided lines."
        )
        
    def find_dimension_by_rate_matrix_parsing(self,start_index: int, file_content: str) -> int:
        # Detect the number of matrix rows based on numeric entries
        n = 0
        current_idx = start_index + 2  # Adjusting to start from matrix values
        while current_idx < len(file_content) and re.search(
            r"(\s*-?\d+\.\d+\s*)+", file_content[current_idx]
        ):
            n += 1
            current_idx += 1
        return n
