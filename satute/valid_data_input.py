
def validate_dir_conflicts(input_args):
    if input_args.dir and any(
        getattr(input_args, opt)
        for opt in ["msa", "tree", "model", "ufboot", "boot", "add_iqtree_options"]
    ):
        raise ValueError(
            "Cannot run Satute with -dir option combined with -'add_iqtree_options', -msa, -tree, -model, -ufboot, or -boot."
        )

def validate_msa_presence(input_args):
    if not input_args.dir and not input_args.msa:
        raise ValueError("MSA file is required when -dir is not used.")
    if input_args.dir and input_args.msa:
        raise ValueError("MSA and dir cannot be used together.")
    
def validate_tree_and_model(input_args):
    if input_args.tree and not input_args.model:
        raise ValueError("Model must be specified when using a tree file.")
    
def validate_boot_options(input_args):
    if input_args.tree and (input_args.ufboot or input_args.boot):
        raise ValueError("Cannot use -ufboot or -boot with -tree.")

def validate_category_range(input_args, number_rates):
    if input_args.model and input_args.category:
        if not (1 <= input_args.category <= number_rates):
            raise ValueError("Chosen category of interest is out of range.")

def validate_and_set_rate_category(input_category: int, number_rates: int):
    """
    Validates the input category against the number of rates and sets the rate category.
    Args:
    - input_category (int): The category input from the user.
    - number_rates (int): The number of rates from the substitution model.
    Returns:
    - str: The validated and set rate category.
    Raises:
    - ValueError: If the input category is out of the valid range.
    """
    if not 1 <= input_category <= number_rates:
        self.logger.error("Chosen category of interest is out of range.")
        raise ValueError("Chosen category of interest is out of range.")
    return str(input_category)
