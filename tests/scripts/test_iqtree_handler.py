import unittest
from unittest.mock import Mock, patch
from satute.handler.iqtree_handler import IqTreeHandler
import subprocess

class TestIqTreeHandler(unittest.TestCase):
    def setUp(self):
        self.iqtree_path = "iqtree2"
        self.logger = Mock()
        self.iqtree_handler = IqTreeHandler(self.iqtree_path, self.logger)

    @patch("subprocess.run")
    def test_run_iqtree_with_arguments_logs_output_and_error(self, mock_run):
        arguments = ["-s", "tests/data/data_aa/cbl3.5_A8+B8_100bp_LG/ali-len100nt_cbl3.5.phy", "-m", "GTR+G"]
        extra_arguments = ["--ufboot", "1000"]
        expected_command = f"{self.iqtree_path} {' '.join(arguments)} {' '.join(extra_arguments)}"

        mock_run.return_value.stdout = "IQ-TREE output"
        mock_run.return_value.stderr = "IQ-TREE error"

        self.iqtree_handler.run_iqtree_with_arguments(arguments, extra_arguments)

        self.logger.info.assert_any_call(f"Running IQ-TREE command: {expected_command}")
        self.logger.error.assert_any_call("IQ-TREE Error: IQ-TREE error")

    @patch("subprocess.run")
    def test_run_iqtree_with_arguments_raises_runtime_error_on_failure(self, mock_run):
        arguments = ["-s","tests/data/data_aa/cbl3.5_A8+B8_100bp_LG/ali-len100nt_cbl3.5.phy" ,"-m", "GTR+G"]
        extra_arguments = ["--ufboot", "1000"]
        expected_command = f"{self.iqtree_path} {' '.join(arguments)} {' '.join(extra_arguments)}"

        mock_run.side_effect = subprocess.CalledProcessError(1, expected_command, output="IQ-TREE output", stderr="IQ-TREE error")

        with self.assertRaises(RuntimeError) as context:
            self.iqtree_handler.run_iqtree_with_arguments(arguments, extra_arguments)

        self.assertEqual(
            str(context.exception),
            f"IQ-TREE execution failed with error code 1.\nCommand: {expected_command}\nOutput: IQ-TREE output\nError: IQ-TREE error\nPlease check the command and ensure that IQ-TREE is installed and accessible."
        )
        self.logger.error.assert_any_call(
            f"IQ-TREE execution failed with error code 1.\nCommand: {expected_command}\nOutput: IQ-TREE output\nError: IQ-TREE error\nPlease check the command and ensure that IQ-TREE is installed and accessible."
        )

if __name__ == "__main__":
    unittest.main()