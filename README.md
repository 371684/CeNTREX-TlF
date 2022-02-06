# CeNTREX-TlF
 
Collection of code tools written for the CeNTREX TlF Schiff moment experiment.

Consists of three main components:
* `states`
* `hamiltonians`
* `couplings`
* `lindblad`

# Installation
To install run `python setup.py install`, which will install the `centrex_TlF` package.

Dependencies are:
* tqdm
* scipy
* numpy
* sympy

# Description

## `states`
`states` contains the functions and classes to represent the TlF states:  
`CoupledBasisState` is a class representing a TlF state with coupled quantum numbers, i.e. F, mF, F1, J, I1, I2, Ω, P.  
`UncoupledBasisState` is a class representing a TlF state with uncoupled quantum numbers, i.e. J, mJ, I1, m1, I2, m2, Ω, P.  
Finally `State` is a class representing a collection of states, since in most cases the TlF molecules are in a superposition state.

The states can be initialized manually, using 
```Python
import centrex_TlF as centrex
centrex.CoupledBasisState(F=1, mF=0, F1 = 1/2, J = 0, I1 = 1/2, I2 = 1/2, Omega = 0, P = 1)
```
or using some of the functions to generate all hyperfine substates in a given J level: 
```Python
import centrex_TlF as centrex
QN = centrex.states.generate_uncoupled_states_ground(Js = [0,1])
```
which returns an array containing the UncoupledBasisStates
```python
array([|X, J = 0, mJ = 0, I₁ = 1/2, m₁ = -1/2, I₂ = 1/2, m₂ = -1/2, P = +, Ω = 0>,
       |X, J = 0, mJ = 0, I₁ = 1/2, m₁ = -1/2, I₂ = 1/2, m₂ = 1/2, P = +, Ω = 0>,
       |X, J = 0, mJ = 0, I₁ = 1/2, m₁ = 1/2, I₂ = 1/2, m₂ = -1/2, P = +, Ω = 0>,
       |X, J = 0, mJ = 0, I₁ = 1/2, m₁ = 1/2, I₂ = 1/2, m₂ = 1/2, P = +, Ω = 0>,
       |X, J = 1, mJ = -1, I₁ = 1/2, m₁ = -1/2, I₂ = 1/2, m₂ = -1/2, P = -, Ω = 0>,
       |X, J = 1, mJ = -1, I₁ = 1/2, m₁ = -1/2, I₂ = 1/2, m₂ = 1/2, P = -, Ω = 0>,
       |X, J = 1, mJ = -1, I₁ = 1/2, m₁ = 1/2, I₂ = 1/2, m₂ = -1/2, P = -, Ω = 0>,
       |X, J = 1, mJ = -1, I₁ = 1/2, m₁ = 1/2, I₂ = 1/2, m₂ = 1/2, P = -, Ω = 0>,
       |X, J = 1, mJ = 0, I₁ = 1/2, m₁ = -1/2, I₂ = 1/2, m₂ = -1/2, P = -, Ω = 0>,
       |X, J = 1, mJ = 0, I₁ = 1/2, m₁ = -1/2, I₂ = 1/2, m₂ = 1/2, P = -, Ω = 0>,
       |X, J = 1, mJ = 0, I₁ = 1/2, m₁ = 1/2, I₂ = 1/2, m₂ = -1/2, P = -, Ω = 0>,
       |X, J = 1, mJ = 0, I₁ = 1/2, m₁ = 1/2, I₂ = 1/2, m₂ = 1/2, P = -, Ω = 0>,
       |X, J = 1, mJ = 1, I₁ = 1/2, m₁ = -1/2, I₂ = 1/2, m₂ = -1/2, P = -, Ω = 0>,
       |X, J = 1, mJ = 1, I₁ = 1/2, m₁ = -1/2, I₂ = 1/2, m₂ = 1/2, P = -, Ω = 0>,
       |X, J = 1, mJ = 1, I₁ = 1/2, m₁ = 1/2, I₂ = 1/2, m₂ = -1/2, P = -, Ω = 0>,
       |X, J = 1, mJ = 1, I₁ = 1/2, m₁ = 1/2, I₂ = 1/2, m₂ = 1/2, P = -, Ω = 0>],
      dtype=object)
```
State objects, which are superpositions of BasisStates are also generated easily:
```Python
superposition = 1*QN[0] + 0.1j*QN[1]
```
which returns
```Python
1.00 x |X, J = 0, mJ = 0, I₁ = 1/2, m₁ = -1/2, I₂ = 1/2, m₂ = -1/2, P = +, Ω = 0>
0.00+0.10j x |X, J = 0, mJ = 0, I₁ = 1/2, m₁ = -1/2, I₂ = 1/2, m₂ = 1/2, P = +, Ω = 0>
```
A subset of `State`, `CoupledBasisStates` can be selected with the `QuantumSelector` as follows:
```Python
QN = centrex.states.generate_coupled_states_ground(Js = [0,1])
qn_select = centrex.states.QuantumSelector(J = 1, mF = 0, electronic = 'X')
qn_select.get_indices(QN)
```
which returns all the indices with `J=1` and `mJ=0`:
```python
array([ 4,  6,  9, 13], dtype=int64)
```
## `hamiltonian`
`hamiltonian` contains the functions to generate TlF hamiltonians in the X and B state in either coupled or uncoupled form.  
Generating a ground state X hamiltonian can be accomplished easily using some convenience functions:
```Python
import centrex_TlF as centrex

# generate the hyperfine sublevels in J=0 and J=1
QN = centrex.states.generate_uncoupled_states_ground(Js = [0,1])

# generate a dictionary with X hamiltonian terms
H = centrex.hamiltonian.generate_uncoupled_hamiltonian_X(QN)

# create a function outputting the hamiltonian as a function of E and B
Hfunc = centrex.hamiltonian.generate_uncoupled_hamiltonian_X_function(H)
```
All functions generating hamiltonians only require a list or array of TlF states. Generating the hamiltonian only for certain hyperfine sublevels is hence also straightforward. The function `calculate_uncoupled_hamiltonian_X` calculates the hamiltonians from scratch, whereas `generate_uncoupled_hamiltonian_X` pulls the non-zero elements from an sqlite database.

To convert a hamiltonian from one basis to another transformation matrices can be generated or calculated
(`generate_transform_matrix` pulls non-zero matrix elements from an sqlite database, `calculate_transform_matrix` does the full element wise calculation):
```Python
import centrex_TlF as centrex

# generate the hyperfine sublevels in J=0 and J=1
QN = centrex.states.generate_uncoupled_states_ground(Js = [0,1])
# generate the coupled hyperfine sublevels in J=0 and J=1
QNc = centrex.states.generate_coupled_states_ground(Js = [0,1])

# generate a dictionary with X hamiltonian terms
H = centrex.hamiltonian.generate_uncoupled_hamiltonian_X(QN)
Hfunc = centrex.hamiltonian.generate_uncoupled_hamiltonian_X_function(H)
H0 = Hfunc(E = [0,0,0], B = [0,0,1e-3])

# generate the transformation matrix
transform = centrex.hamiltonian.generate_transform_matrix(QN, QNc)

# calculate the transformed matrix
H0c = transform.conj().T@H0@transform
```
This is mostly used for optical bloch simulations where the coupled states representation is more convenient.

## `couplings`
`couplings` contains functions to generate coupling matrices and collapse matrices for optical bloch simulations, plus some other helpful tools.  

Transitions are defined with `TransitionSelectors`, which specify the ground states, excited states, polarizations, and symbols used to show the states:
```python
import sympy as smp

gnd = centrex.states.QuantumSelector(J=[1,2], electronic = 'X')
gnd_laser = centrex.states.QuantumSelector(J=1, electronic = 'X')
exc = centrex.states.QuantumSelector(J=1, F=1, F1=1/2, electronic = 'B', P = +1)

transitions = [
    centrex.couplings.TransitionSelector(
        ground = 1*centrex.states.generate_coupled_states_ground_X(gnd_laser),
        excited = 1*centrex.states.generate_coupled_states_excited_B(exc),
        polarizations = [[1,0,0], [0,0,1]],
        polarization_symbols = smp.symbols("Plx Plz"),
        Ω = smp.Symbol('Ωl', complex = True),
        δ = smp.Symbol('δl'),
        description = "laser transition"
    )
]
```
This specifies a transition between the `X` states of `J=1` to the excited `B` state of `|J=1, F=1, F1=1/2>`. The polarizations included are `x` (`[1,0,0]`) and `z` (`[0,0,1]`), with symbolic representation `Plx` and `Plz`. The transition Rabi rate is `Ωl` and the detuning is `δl`. These representations are used for the symbolic matrices and for the optical-bloch equation solver in `linblad`, where all these symbolic parameters can be set to a fixed value or scanned over.

`generate_coupling_field_automatic` generates the matrices that specify the couplings between the states given in the `TransitionSelector` object. To calculate this you need to now the exact states included in the system, calculated with `generate_total_reduced_hamiltonian` (from `hamiltonian`):
```python
ground_states, excited_states, QN, H_int, V_ref_int = \
        centrex.hamiltonian.generate_total_reduced_hamiltonian(
            ground_states_approx  = \
                    centrex.states.generate_coupled_states_ground_X(gnd),
            excited_states_approx = \
                    centrex.states.generate_coupled_states_excited_B(exc)
        )
```
Note that this by default calculates the Hamiltonian in zero field.  
Then the couplings are calculated with:
```python
couplings = [
        centrex.couplings.generate_coupling_field_automatic(
            transition.ground,
            transition.excited,
            H_int, QN, V_ref_int, 
            pol_vec = transition.polarizations,
            nprocs = 4)
        for transition in transitions
    ]
```
`couplings` is a list of dictionaries, with a dictionary for each transition (in this case only 1). The dictionary contains the following keys:
```
ground_main       : main ground state (i.e. defining 0 detuning)
excited_main      : main excited state (i.e. defining 0 detuning)
main_coupling     : coupling element size between ground_main and excited_main
ground_states     : ground states involved in the coupling
excited_states    : excited states involved in the coupling
fields            : list of dictionaries, with a single dictionary per polarization.
                    Each dict contains the polarization vector ("pol") and coupling field ("field")
```
The indices of the coupling fields correspond to equivalent indices of the supplied `QN`. 

The collapse matrices are generated with `collapse_matrices`
```Python
import numpy as np
Γ = 2*np.pi*1.56e6
C_array = centrex.couplings.collapse_matrices(
            QN, ground_states, excited_states, gamma = Γ
        )
```
Collapse matrices returns an `n x m x m` array, where `m` is the number of states involved, and `n` is the total number of decays; e.g. `C_array[0]` represents the decay for a single state, `C_array[1]` for the next state, etc.
## `lindblad`
This contains all the function for generating and solving the Optical Bloch Equation (OBE) system. The actual solving of the systems of equations is done with `Julia`, specifically `DifferentialEquations.jl`, which is many times faster than either `SciPy` or equivalent `Python` ODE solver packages. 

For ease of use a unified interface is defined, which we'll go through step-by-step
```Python
import numpy as np
import centrex_TlF as centrex
from sympy import symbols, Symbol

syspars = centrex.states.SystemParameters(
    nprocs  = 6,
    Γ       = 2*np.pi*1.56e6,
    X  = [centrex.states.QuantumSelector(J=1, electronic = 'X'),
               centrex.states.QuantumSelector(J=3, electronic = 'X')],
    B = centrex.states.QuantumSelector(J=1, F=1, F1=3/2, P=+1, electronic = 'B')
)
```
`SystemParameters` defines the general system parameters used in the simulations  
* nprocs -> the number of processes to use in generating and solving the system; generally set equal to the number of processor cores
* Γ -> decay rate in 2π⋅Hz
* X -> involved X states
* B -> involved B states
  

```Python
transitions = [
    centrex.couplings.TransitionSelector(
        ground               = 1*centrex.states.generate_coupled_states_ground_X(syspars.B[0]),
        excited              = 1*centrex.states.generate_coupled_states_excited_B(syspars.X),
        polarizations        =        [[1,0,0], [0,0,1]],
        polarization_symbols = symbols("Plx     Plz"),
        Ω                    = Symbol('Ωl', complex = True),
        δ                    = Symbol('δl'),
        description          = "laser transition"
    )
]
```
A laser transition between `J=1` and the excited state is defined here, see `transitions` for more detail. 
```Python
odepars = centrex.lindblad.odeParameters(
    Ωl = "Ωl0 * phase_modulation(t, β, ωphase)",
    Ωl0     = 1*syspars.Γ,    # Rabi frequency of the laser [rad/s]
    δl      = 0,              # detuning of the laser [rad/s]
    
    # laser phase modulation
    ωphase = syspars.Γ,       # laser phase modulation frequency [rad/s]
    β      = 3.8,             # laser phase modulation depth [rad]

    # laser polarization switching
    ωp = 2*np.pi*1.56e6,           # polarization switching frequency [rad/s]
    φp = 0.0,                 # polarization switching phase [rad]
    Pl  = "sin(ωp*t + φp)",
    Plz = "Pl>0",
    Plx = "Pl<0",
    
    # molecules
    y0 = 0,                   # molecule start y position [m]
    vz = 184,                 # longitudinal molecular velocity [m/s]
    vy = 0,                   # molecule vy [m/s]
)
```
`odeParameters` defines the parameters used in the OBE system. I.e. the Rabi rate, detuning and whatever else you'd like to define. Note that compound expressions are allowed as well, here the polarization switches between `x` and `z` with a frequency of `ωp`, and the laser is phase modulated.
Upon initialization, `odeParameters` automatically checks whether all defined variables are specified, including variables defined in compound expressions such as the Rabi rate and polarizations. 

Note that the symbols defined in the `TransitionSelector` objects have to be defined here; the function generating the system of equations checks if these symbols are defined and throws an error otherwise.

Finally `setup_OBE_system` generates the system, transcribes it into `julia` code and initializes all necessary `julia` packages and processes in the background.
```Python
obe_system = centrex.lindblad.setup_OBE_system_julia(
                  syspars, odepars, transitions, verbose=True, full_output=True, 
                  qn_compact=centrex.states.QuantumSelector(J=3, electronic = 'X')
                )
```
Note specifically `qn_compact`; for some simulations certain decays have to be included, but you might not care about which specific hyperfine sublevel the decays go to. `qn_compact` specifies which quantum numbers to combine into a single level (in this case `J=3`) in order to speed up the simulations.
Additionally other decay paths (such as vibrational decays, for instance) can be added with the optional keyword argument `decay_channels`, which expects a `DecayChannel` or list/array of `DecayChannel`. `DecayChannel` is a dataclass of format:
```Python
@dataclass
class DecayChannel:
    ground: CoupledBasisState
    excited: np.ndarray
    branching: float
    description: str = ""
```
`ground`, a `CoupledBasisState`, indicates the decay path, `excited` specifies which excited states decay to `ground`.
`branching` is the branching ration to the specified ground state, e.g. 0.04 for 4%. `description` is just an optional string to describe the decay channel.

To then run a single trajectory the integration time then has to be specified and initial density:
```Python
tspan = [0,110e-6]
ρ = centrex.utils.generate_population_states([0], len(obe_system.QN))
```
Here all initial population is put into
```python
1.00+0.00j x |X, J = 1, F₁ = 1/2, F = 0, mF = 0, I₁ = 1/2, I₂ = 1/2, P = -, Ω = 0>
```
and the integration is run from `0 -> 110 us`.  
The simulation is run with
```Python
t_array_compact, pop_results_compact = centrex.lindblad.do_simulation_single(odepars, tspan, ρ)
```
`t_array` `(n,)` contains the timesteps corresponding to the solutions `pop_results_compact` `(m x n)`. Here `m` is axis corresponding to the included states `obe_system.QN`.