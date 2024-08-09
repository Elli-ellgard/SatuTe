import pandas as pd
from typing import Dict, Any, List


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
            "coherence": self.coefficients,
            "category_variance": self.variances,
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
                data_row = {"branch": identifier}
                data_row.update(row.to_dict())
                data.append(data_row)

        return pd.DataFrame(data)

