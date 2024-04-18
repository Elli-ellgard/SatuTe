from typing import Dict, Any, List
import pandas as pd


class TestStatisticComponents:
    def __init__(self, coefficients: List[float], variances: List[float]):
        """
        Initializes the TestStatisticComponents with coefficients and variances.

        Args:
            coefficients (list[float]): A list of coefficients for the test statistic components.
            variances (list[float]): A list of variances corresponding to the test statistic components.
        """
        self.coefficients = coefficients
        self.variances = variances

    def to_dataframe(self) -> pd.DataFrame:
        """
        Converts the test statistic components into a pandas DataFrame.

        Returns:
            pd.DataFrame: A DataFrame where each row represents a test statistic component, with columns for 'Coefficient' and 'Variance'.
        """
        # Ensure the lists are of equal length to avoid data misalignment
        if len(self.coefficients) != len(self.variances):
            raise ValueError(
                "Coefficients and variances lists must be of the same length."
            )

        # Create a DataFrame from the components
        data = {
            "coefficient": self.coefficients,
            "sample_variance": self.variances,
        }
        df = pd.DataFrame(data)

        return df


class TestStatisticComponentsContainer:
    """
    Manages multiple instances of TestStatisticComponents.

    This class acts as a container for managing a collection of TestStatisticComponents instances, facilitating
    the addition, retrieval, and collective handling of multiple sets of test statistic components across different statistical tests or scenarios.

    Attributes:
        components (Dict[str, TestStatisticComponents]): A dictionary storing TestStatisticComponents instances, keyed by an identifier.
    """

    def __init__(self) -> None:
        """
        Initializes a new instance of TestStatisticComponentsContainer with an empty dictionary to store TestStatisticComponents instances.
        """
        self.components: Dict[str, TestStatisticComponents] = {}

    def add_component(
        self, identifier: str, test_statistic_component: TestStatisticComponents
    ) -> None:
        """
        Adds a TestStatisticComponents instance to the collection under a specified identifier.

        Args:
            identifier: The identifier to associate with the TestStatisticComponents instance.
            test_statistic_component: The TestStatisticComponents instance to be added.
        """
        self.components[identifier] = test_statistic_component

    def get_component(self, identifier: str) -> TestStatisticComponents:
        """
        Retrieves a TestStatisticComponents instance by its identifier.

        Args:
            identifier: The identifier of the component to retrieve.

        Returns:
            The TestStatisticComponents instance associated with the given identifier, or None if no such component exists.
        """
        return self.components.get(identifier)

    def get_all_components(self) -> Dict[str, TestStatisticComponents]:
        """
        Returns all stored TestStatisticComponents instances.

        Returns:
            A dictionary of all TestStatisticComponents instances, keyed by their associated identifiers.
        """
        return self.components

    def to_dataframe(self) -> pd.DataFrame:
        """
        Converts the stored TestStatisticComponents instances into a pandas DataFrame.

        Returns:
            pd.DataFrame: A DataFrame containing all test statistic components, with each row representing a component from an instance,
            including an 'Identifier' column to distinguish between different TestStatisticComponents instances.
        """
        data: List[Dict[str, Any]] = []

        for identifier, component in self.components.items():
            df = component.to_dataframe()
            for _, row in df.iterrows():
                data_row = {"Edge": identifier}
                data_row.update(row.to_dict())
                data.append(data_row)

        return pd.DataFrame(data)


class TestResultBranch:
    """
    Represents individual test results for saturation tests.

    This class provides a structured way to store and manage test results, where each test is identified by a unique name
    and associated with a numerical score indicating the test outcome.

    Attributes:
        results (Dict[str, float]): A dictionary mapping test names to their scores.
    """

    def __init__(self, **results: float) -> None:
        """
        Initializes a new instance of TestResultBranch, optionally with initial test results.

        Args:
            **results: Variable length keyword arguments, where each key-value pair represents a test name and its score.
        """
        self.results: Dict[str, Any] = results

    def add_result(self, result_name: str, score: Any) -> None:
        """
        Adds or updates a test result.

        If the test result already exists, its score will be updated; otherwise, a new test result entry is created.

        Args:
            result_name: The unique name of the test.
            score: The score of the test.
        """
        self.results[result_name] = score

    def get_results(self) -> Dict[str, Dict[str, Any]]:
        """
        Retrieves all test results.

        Returns:
            A dictionary of all stored test results, with test names as keys and their scores as values.
        """
        return self.results


class TestResultsBranches:
    """
    Manages multiple instances of TestResultBranch.

    This class acts as a container for managing a collection of TestResultBranch instances, facilitating
    the addition, retrieval, and collective handling of multiple sets of test results across different branches.

    Attributes:
        branches (Dict[str, TestResultBranch]): A dictionary storing TestResultBranch instances, keyed by branch name.
    """

    def __init__(self) -> None:
        """
        Initializes a new instance of TestResultsBranches with an empty dictionary to store TestResultBranch instances.
        """
        self.branches: Dict[str, TestResultBranch] = {}

    def add_branch(
        self, branch_name: str, test_result_branch: TestResultBranch
    ) -> None:
        """
        Adds a TestResultBranch instance to the collection under a specified branch name.

        Ensures that the added instance is indeed a TestResultBranch to maintain the integrity of the container.

        Args:
            branch_name: The name to associate with the TestResultBranch instance.
            test_result_branch: The TestResultBranch instance to be added.

        Raises:
            ValueError: If the provided instance is not an instance of TestResultBranch.
        """
        if not isinstance(test_result_branch, TestResultBranch):
            raise ValueError(
                "test_result_branch must be an instance of TestResultBranch"
            )
        self.branches[branch_name] = test_result_branch

    def get_branch(self, branch_name: str) -> TestResultBranch:
        """
        Retrieves a TestResultBranch instance by its branch name.

        Args:
            branch_name: The name of the branch to retrieve.

        Returns:
            The TestResultBranch instance associated with the given branch name, or None if no such branch exists.
        """
        return self.branches.get(branch_name)

    def get_all_branches(self) -> Dict[str, TestResultBranch]:
        """
        Returns all stored TestResultBranch instances.

        Returns:
            A dictionary of all TestResultBranch instances, keyed by their associated branch names.
        """
        return self.branches

    def to_dataframe(self) -> pd.DataFrame:
        """
        Converts the stored TestResultBranch instances into a pandas DataFrame.

        This method dynamically creates a DataFrame structure. Each row in the DataFrame represents the results
        from a single branch. The DataFrame will have columns dynamically adjusted based on the content of the results,
        with a fixed column for 'Branch Name'.

        Returns:
            pd.DataFrame: A DataFrame containing all test results, structured dynamically based on the results content.
        """
        # List to accumulate row data for the DataFrame
        data: List[Dict[str, Any]] = []

        # Iterate through each branch
        for branch_name, test_result_branch in self.branches.items():
            # Retrieve the results for the current branch
            results = test_result_branch.get_results()

            # Assuming each branch's results are a dict that you want to turn into a single row
            # Add the branch name to the row
            row = {"edge": branch_name}

            # Assuming results are structured in a way that can directly fill into row data
            row.update(results)

            # Append the constructed row to the data list
            data.append(row)

        # Convert the data list into a DataFrame
        return pd.DataFrame(data)
