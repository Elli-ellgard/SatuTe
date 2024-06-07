import logging
import unittest
from unittest.mock import MagicMock
from satute.valid_data_input import validate_and_check_rate_categories

class TestValidateRateCategories(unittest.TestCase):
    def setUp(self):
        # Setup a mock logger
        self.logger = logging.getLogger('TestLogger')
        self.logger.addHandler(logging.NullHandler())

    def test_valid_category(self):
        categorized_sites = {
            'p1': ['site1', 'site2'],
            'p2': [],
            'p3': ['site3']
        }
        chosen_category = 1
        
        # This should pass without exception
        validate_and_check_rate_categories(categorized_sites, chosen_category, self.logger)

    def test_empty_category_warning(self):
        categorized_sites = {
            'p1': ['site1', 'site2'],
            'p2': [],
            'p3': ['site3']
        }
        chosen_category = 3

        # Mock logger method
        self.logger.warning = MagicMock()

        validate_and_check_rate_categories(categorized_sites, chosen_category, self.logger)
        
        # Check if warning was called for empty category 2
        self.logger.warning.assert_called_with("Will be skipping Rate category p2")

    def test_invalid_category_error(self):
        categorized_sites = {
            'p1': ['site1', 'site2'],
            'p2': [],
            'p3': ['site3']
        }
        chosen_category = 2

        # Expect ValueError to be raised for empty chosen category
        with self.assertRaises(ValueError) as context:
            validate_and_check_rate_categories(categorized_sites, chosen_category, self.logger)
        
        self.assertIn(
            f"Chosen category rate {chosen_category}",
            str(context.exception)
        )

if __name__ == '__main__':
    unittest.main()