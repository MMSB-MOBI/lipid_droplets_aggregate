class MoleculeNotFound(Exception):
    """Raised when molecule not found"""
    pass

class InvalidAxis(Exception):
    """Raised when axis has not valid value"""
    pass

class ArgumentError(Exception):
    """Raised when error in command line argument"""
    pass

class NotComputedError(Exception):
    """Raised when something is not computed"""
    pass
