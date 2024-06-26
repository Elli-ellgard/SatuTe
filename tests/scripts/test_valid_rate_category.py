import unittest
from unittest.mock import MagicMock, call
from logging import Logger
from satute.valid_data_input import validate_and_check_rate_categories

class TestValidateRateCategories(unittest.TestCase):
    def setUp(self):
        # Setup a mock logger
        self.logger = MagicMock(spec=Logger)

    def test_valid_category(self):
        categorized_sites = {
            '1': ['site1', 'site2'],
            '2': [],
            '3': ['site3']
        }
        chosen_category = 1
        
        # This should pass without exception
        validate_and_check_rate_categories(categorized_sites, chosen_category, self.logger)

        # Ensure no warning is called for the chosen category
        self.logger.warning.assert_called_once_with("Skipping empty rate category '2'.")

    def test_empty_category_warning(self):
        categorized_sites = {
            '1': ['site1', 'site2'],
            '2': [],
            '3': ['site3']
        }
        chosen_category = 3

        validate_and_check_rate_categories(categorized_sites, chosen_category, self.logger)
        
        # Check if warning was called for empty category 2
        self.logger.warning.assert_called_once_with("Skipping empty rate category '2'.")

    def test_invalid_category_error(self):
        categorized_sites = {
            '1': ['site1', 'site2'],
            '2': [],
            '3': ['site3']
        }
        chosen_category = 2

        # Expect ValueError to be raised for empty chosen category
        with self.assertRaises(ValueError) as context:
            validate_and_check_rate_categories(categorized_sites, chosen_category, self.logger)
        
        self.assertIn(
            f"Chosen category rate '{chosen_category}' is empty. Choose a different category with assigned sites.",
            str(context.exception)
        )

    def test_non_existent_category_no_error(self):
        categorized_sites = {
            '1': ['site1', 'site2'],
            '2': [],
            '3': ['site3']
        }
        chosen_category = 4

        # Should not raise an error for non-existent category
        validate_and_check_rate_categories(categorized_sites, chosen_category, self.logger)

        # Ensure no warning is logged for non-existent chosen category
        self.logger.warning.assert_called_once_with("Skipping empty rate category '2'.")

    def test_valid_empty_category_no_error(self):
        categorized_sites = {
            '1': [],
            '2': ['site1', 'site2'],
            '3': []
        }
        chosen_category = 2

        # This should pass without exception because category 2 is valid and not empty
        validate_and_check_rate_categories(categorized_sites, chosen_category, self.logger)

        # Ensure warnings were called for empty categories
        expected_warnings = [
            call("Skipping empty rate category '1'."),
            call("Skipping empty rate category '3'.")
        ]
        self.logger.warning.assert_has_calls(expected_warnings, any_order=True)

if __name__ == '__main__':
    unittest.main()