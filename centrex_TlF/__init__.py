import logging

from rich.logging import RichHandler

from . import states
from . import hamiltonian
from . import couplings
from . import constants
from . import lindblad
from . import utils
from . import transitions

from .states.states import State, UncoupledBasisState, CoupledBasisState

FORMAT = "%(message)s"
logging.basicConfig(
    level="WARNING", format=FORMAT, datefmt="[%X]", handlers=[RichHandler()]
)
