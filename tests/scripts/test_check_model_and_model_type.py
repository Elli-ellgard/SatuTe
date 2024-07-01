import unittest

class TestCheckModelAndType(unittest.TestCase):

    def setUp(self):
        # Initialize any resources needed for the tests
        pass

    def tearDown(self):
        # Clean up any resources used in the tests
        pass

    def test_check_model(self):
        from satute.iqtree_parser import IqTreeParser
        test_instance = IqTreeParser()

        # Test with a not accepted Protein model
        aa_models = ["NQ.BIRD+F", "NQ.PLANT"]
        for model in aa_models:
            with self.assertRaises(ValueError):
                test_instance.check_model(model)

        # Test with a not accepted DNA model
        dna_models = ["12.12+R4", "6.8a"]
        for model in dna_models:
            with self.assertRaises(ValueError):
                test_instance.check_model(model)

        # Test with an accepted model (neither not accepted DNA nor not accepted protein)
        models = ["GTR+G4", "FLU+FQ", "000000"]
        for model in models: 
            result = test_instance.check_model(model)
            self.assertIsNone(result)  

 
    def test_get_model_type(self):
        from satute.iqtree_parser import (
            IqTreeParser,
            ModelType
        )
        test_instance = IqTreeParser()

        # Test with a protein model
        result = test_instance.get_model_type("FLU+FQ")
        self.assertEqual(result, ModelType.PROTEIN)

        # Test with a DNA model
        result = test_instance.get_model_type("GTR+I+G4")
        self.assertEqual(result,  ModelType.DNA)

if __name__ == '__main__':
    unittest.main()
