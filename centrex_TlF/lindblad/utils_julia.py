import numpy as np
import sympy as smp
from julia import Main
from pathlib import Path
from sympy import Symbol
from sympy.utilities.lambdify import lambdify

__all__ = [
    "initialize_julia",
    "generate_ode_fun_julia",
    "setup_variables_julia",
    "odeParameters",
    "setup_parameter_scan_1D",
    "setup_ratio_calculation",
    "setup_initial_condition_scan",
    "setup_state_integral_calculation",
    "get_indices_diag_flattened",
    "setup_parameter_scan_ND",
    "handle_randomized_ensemble_solution",
    "setup_discrete_callback_terminate",
    "setup_problem",
    "solve_problem",
    "get_results",
    "do_simulation_single",
    "setup_problem_parameter_scan",
    "solve_problem_parameter_scan",
    "get_results_parameter_scan",
]


def initialize_julia(nprocs):
    """
    Function to initialize Julia over nprocs processes.
    Creates nprocs processes and loads the necessary Julia
    packages.

    Args:
        nprocs (int): number of Julia processes to initialize.
    """
    Main.eval(
        """
        using Logging: global_logger
        using TerminalLoggers: TerminalLogger
        global_logger(TerminalLogger())

        using Distributed
        using ProgressMeter
    """
    )

    if Main.eval("nprocs()") < nprocs:
        Main.eval(f"addprocs({nprocs}-nprocs())")

    if Main.eval("nprocs()") > nprocs:
        procs = Main.eval("procs()")
        procs = procs[nprocs:]
        Main.eval(f"rmprocs({procs})")

    Main.eval(
        """
        @everywhere begin
            using LinearAlgebra
            using Trapz
            using DifferentialEquations
        end
    """
    )
    # loading common julia functions from julia_common.jl
    path = Path(__file__).parent / "julia_common.jl"
    Main.eval(f'include(raw"{path}")')

    print(f"Initialized Julia with {nprocs} processes")


def generate_ode_fun_julia(preamble, code_lines):
    """
    Generate the ODE function from the preamble and code lines
    generated in Python.

    Args:
        preamble (str): preamble of the ODE function initializing the 
                        function variable definitions.
        code_lines (list): list of strings, each line is a generated
                            line of Julia code for part of the ODE.

    Returns:
        str : function definition of the ODE 
    """
    ode_fun = preamble
    for cline in code_lines:
        ode_fun += "\t\t" + cline + "\n"
    ode_fun += "\t end \n \t nothing \n end"
    Main.eval(f"@everywhere {ode_fun}")
    return ode_fun


def setup_variables_julia(Γ, ρ, vars=None):
    """
    Convenience function for loading Γ, ρ and specified vars 
    in Julia memory.

    Args:
        Γ (float): Γ of the simulated system
        ρ (np.ndarray): complex density matrix of the simulated
                        system
        vars (dict): dict of variables (str key) to load in
                    Julia memory
    """
    Main.Γ = Γ
    Main.ρ = ρ
    Main.eval(
        """
        @everywhere begin
            Γ = $Γ
            ρ = $ρ
        end
    """
    )
    if vars:
        for key, val in vars.items():
            Main.eval(f"@everywhere {key} = {val}")


type_conv = {
    "int": "Int64",
    "float": "Float64",
    "complex": "ComplexF64",
    "float64": "Float64",
    "int32": "Int32",
    "complex128": "ComplexF64",
    "list": "Array",
    "ndarray": "Array",
}
julia_funcs = [
    "gaussian_2d",
    "phase_modulation",
    "square_wave",
    "multipass_2d_intensity",
]


class odeParameters:
    def __init__(self, *args, **kwargs):
        # if elif statement is for legacy support, where a list of parameters was supplied
        if len(args) > 1:
            raise AssertionError(
                "For legacy support supply a list of strings, one for each parameter"
            )
        elif len(args) == 1:
            assert isinstance(
                args[0][0], str
            ), "For legacy support supply a list of strings, one for each parameter"
            if "ρ" not in args[0]:
                args[0].append("ρ")
            kwargs = {par: 0.0 for par in args[0]}
            odeParameters(**kwargs)

        # kwargs = self._check_for_density(kwargs)
        # kwargs = self._check_for_states(kwargs)
        self._parameters = [
            key for key, val in kwargs.items() if not isinstance(val, str)
        ]
        self._compound_vars = [
            key for key, val in kwargs.items() if isinstance(val, str)
        ]

        for key, val in kwargs.items():
            # ϕ = 'ϕ') results in different unicode representation
            # replace the rhs with the rep. used on the lhs
            if isinstance(val, str):
                val = val.replace("\u03d5", "\u03c6")
            setattr(self, key, val)

        self._check_symbols_defined()
        self._order_compound_vars()

        # checking types, necessary if the ODE parameters contain arrays or list
        # Julia can't do type inference then and this tanks performance
        # storing the input types here for use in generate_preamble, but this is
        # only used if one of the input parameters is an array
        self._parameter_types = [
            type_conv.get(type(getattr(self, par)).__name__) for par in self._parameters
        ]
        self._array_types = dict(
            (par, type_conv.get(type(getattr(self, par)[0]).__name__))
            for par in self._parameters
            if type_conv.get(type(getattr(self, par)).__name__) == "Array"
        )

    def _check_for_density(self, kwargs):
        assert "ρ" in kwargs.keys(), "Supply an initial density ρ to odeParameters"
        self.ρ = kwargs.get("ρ")
        del kwargs["ρ"]
        return kwargs

    def _check_for_states(self, kwargs):
        assert (
            "ground" in kwargs.keys()
        ), "Supply ground states `ground` to odeParameters"
        self.ground = kwargs.get("ground")
        del kwargs["ground"]
        assert (
            "excited" in kwargs.keys()
        ), "Supply excited states `excited` to odeParameters"
        self.excited = kwargs.get("excited")
        del kwargs["excited"]
        return kwargs

    def __setattr__(self, name, value):
        if name in [
            "_parameters",
            "_compound_vars",
            "_parameter_types",
            "_array_types",
        ]:
            super(odeParameters, self).__setattr__(name, value)
        elif name == "ρ":
            super(odeParameters, self).__setattr__(name, value)
        elif name in ["ground", "excited"]:
            super(odeParameters, self).__setattr__(name, value)
        elif name in self._parameters:
            assert not isinstance(
                value, str
            ), "Cannot change parameter from numeric to str"
            super(odeParameters, self).__setattr__(name, value)
        elif name in self._compound_vars:
            assert isinstance(value, str), "Cannot change parameter from str to numeric"
            super(odeParameters, self).__setattr__(name, value)
        else:
            raise AssertionError(
                "Cannot instantiate new parameter on initialized OdeParameters object"
            )

    def _get_defined_symbols(self):
        symbols_defined = self._parameters + self._compound_vars
        symbols_defined += ["t"]
        symbols_defined = set([Symbol(s) for s in symbols_defined])
        return symbols_defined

    def _get_numerical_symbols(self):
        symbols_numerical = self._parameters[:]
        symbols_numerical += ["t"]
        symbols_numerical = set([Symbol(s) for s in symbols_numerical])
        return symbols_numerical

    def _get_expression_symbols(self):
        symbols_expressions = [
            smp.parsing.sympy_parser.parse_expr(getattr(self, s))
            for s in self._compound_vars
        ]
        symbols_expressions = set().union(
            *[s.free_symbols for s in symbols_expressions]
        )
        return symbols_expressions

    def _check_symbols_defined(self):
        symbols_defined = self._get_defined_symbols()
        symbols_expressions = self._get_expression_symbols()

        warn_flag = False
        warn_string = f"Symbol(s) not defined: "
        for se in symbols_expressions:
            if se not in symbols_defined:
                warn_flag = True
                warn_string += f"{se}, "
        if warn_flag:
            raise AssertionError(warn_string.strip(" ,"))

    def check_symbols_in_parameters(self, symbols_other):
        if not isinstance(symbols_other, (list, np.ndarray, tuple, set)):
            symbols_other = [symbols_other]
        elif isinstance(symbols_other, set):
            symbols_other = list(symbols_other)
        if isinstance(symbols_other[0], smp.Symbol):
            symbols_other = [str(sym) for sym in symbols_other]

        parameters = self._parameters[:]
        parameters += ["t"]

        warn_flag = False
        warn_string = f"Symbol(s) not defined: "
        for se in symbols_other:
            if se not in parameters:
                warn_flag = True
                warn_string += f"{se}, "
        if warn_flag:
            raise AssertionError(warn_string.strip(" ,"))

    def _order_compound_vars(self):
        symbols_num = self._get_numerical_symbols()
        unordered = list(self._compound_vars)
        ordered = []

        while len(unordered) != 0:
            for compound in unordered:
                if compound not in ordered:
                    symbols = smp.parsing.sympy_parser.parse_expr(
                        getattr(self, compound)
                    ).free_symbols
                    m = [
                        True if (s in symbols_num) or (str(s) in ordered) else False
                        for s in symbols
                    ]
                    if all(m):
                        ordered.append(compound)
            unordered = [val for val in unordered if val not in ordered]
        self._compound_vars = ordered

    def _get_index_parameter(self, parameter, mode="python"):
        # OdeParameter(ϕ = 'ϕ') results in different unicode representation
        # replace the rhs with the rep. used on the lhs
        parameter = parameter.replace("\u03d5", "\u03c6")
        if mode == "python":
            return self._parameters.index(parameter)
        elif mode == "julia":
            return self._parameters.index(parameter) + 1

    @property
    def p(self):
        return [getattr(self, p) for p in self._parameters]

    def get_index_parameter(self, parameter, mode="python"):
        if isinstance(parameter, str):
            return self._get_index_parameter(parameter, mode)
        elif isinstance(parameter, (list, np.ndarray)):
            return [self._get_index_parameter(par, mode) for par in parameter]

    def check_transition_symbols(self, transitions):
        # check Rabi rate and detuning symbols
        to_check = ["Ω", "δ"]
        symbols_defined = [str(s) for s in self._get_defined_symbols()]
        not_defined = []
        for transition in transitions:
            for var in to_check:
                var = str(getattr(transition, var))
                if var is not None:
                    if var not in symbols_defined:
                        not_defined.append(var)
        if len(not_defined) > 0:
            not_defined = set([str(s) for s in not_defined])
            raise AssertionError(
                f"Symbol(s) from transitions not defined: {', '.join(not_defined)}"
            )

        # check polarization switching symbols
        to_check = []
        for transition in transitions:
            to_check.extend(transition.polarization_symbols)

        symbols_defined = self._get_defined_symbols()

        warn_flag = False
        warn_string = f"Symbol(s) in transition polarization switching not defined: "
        for ch in to_check:
            if ch not in symbols_defined:
                warn_flag = True
                warn_string += f"{ch}, "
        if warn_flag:
            raise AssertionError(warn_string.strip(" ,"))
        return True

    def generate_p_julia(self):
        Main.eval(f"p = {self.p}")

    def __repr__(self):
        rep = f"OdeParameters("
        for par in self._parameters:
            rep += f"{par}={getattr(self, par)}, "
        return rep.strip(", ") + ")"

    def get_parameter_evolution(self, t, parameter):
        """Get the time evolution of parameters in odeParameters.
        Evaluates expressions in python if possible, otherwise calls julia to 
        evaluate the expressions

        Args:
            t (np.ndarray[float]): array of timestamps
            parameter (str): parameter to evaluate over t

        Returns:
            np.ndarray : array with values of parameter corresponding to t

        """
        # get a list of all parameters in odeParameters
        parameters = self._parameters + self._compound_vars
        # check if `parameter` is defined in odeParameters
        assert parameter in parameters, f"{parameter} not defined in odeParameters"
        # if `parameter` is not a compound variable, an array of size t of
        # parameter
        if parameter in self._parameters:
            return np.ones(len(t)) * getattr(self, str(parameter))
        elif parameter in self._compound_vars:
            expression = getattr(self, str(parameter))
            # parse expression string to sympy equation
            expression = smp.parsing.sympy_parser.parse_expr(expression)
            while True:
                symbols_in_expression = [
                    sym for sym in expression.free_symbols if sym is not smp.Symbol("t")
                ]
                symbols_in_expression = [str(sym) for sym in symbols_in_expression]
                compound_bool = [
                    sym in self._compound_vars for sym in symbols_in_expression
                ]
                # if any of the symbols in the expression are compound variables
                # fetch the related compound expression and insert it in the
                # expression
                if np.any(compound_bool):
                    for idx in np.where(compound_bool)[0]:
                        compound_var = smp.parsing.sympy_parser.parse_expr(
                            getattr(self, symbols_in_expression[idx])
                        )
                        expression = expression.subs(
                            symbols_in_expression[idx], compound_var
                        )
                else:
                    # break if none of the symbols in the expression are compound
                    # variables
                    break

            # substitute the numerical variables for the expressions
            symbols_in_expression = [
                sym for sym in expression.free_symbols if sym is not smp.Symbol("t")
            ]
            array_symbols = []
            for sym in symbols_in_expression:
                # can't subs a symbol with a list, tuple or array
                if isinstance(getattr(self, str(sym)), (list, tuple, np.ndarray)):
                    array_symbols.append(sym)
                expression = expression.subs(sym, getattr(self, str(sym)))
            functions_in_expression = [
                str(fn).split("(")[0] for fn in expression.atoms(smp.Function)
            ]
            # check if any of the functions are special julia defined functions
            if np.any([fn in julia_funcs for fn in functions_in_expression]):
                expression = str(expression)
                # broadcast the function, allows for input of an array of t
                for fn in julia_funcs:
                    expression = expression.replace(fn, f"{fn}.")
                for sym in array_symbols:
                    expression = expression.replace(
                        str(sym), f"Ref({getattr(self, str(sym))})"
                    )
                # quick workaround to convert int to float for function inputs
                # LOOK AT THIS LATER
                expression = expression.replace(", 0,", ", 0.0,")
                # evaluate the specified parameter expression in julia
                Main.eval(f"_tmp_func(t) = {str(expression)}")
                # can't get broadcasting to work if some variables are of array or list type, use map
                Main.tmp_t = t
                return Main.eval(f"map(_tmp_func, tmp_t)")
            else:
                # evaluate the specified parameter expression in python
                func = lambdify(smp.Symbol("t"), expression, modules=["numpy", "scipy"])
                res = func(t)
                if np.shape(res) == ():
                    return np.ones(len(t)) * res
                else:
                    return res


def setup_parameter_scan_1D(odePar, parameter, values):
    """
    Convenience function for setting up a 1D parameter scan.
    Scan can be performed over multiple parameters simultaneously,
    but only for the same value for each parameter. For different
    values per parameter see setup_parameter_scan_zipped.

    Args:
        odePar (odeParameters): object containing all the parameters
                                for the ODE system
        parameters (str, list): parameter or list of parameters to
                                scan over
        values (list, np.ndarray): values to scan the parameter(s)
                                    over.
    """
    # check if parameters is a list, get indices of the parameters
    # as defined in odePar
    if isinstance(parameter, (list, tuple)):
        indices = [odePar.get_index_parameter(par) for par in parameter]
    else:
        indices = [odePar.get_index_parameter(parameter)]

    # generate the parameter sequence for the prob_func function
    pars = list(odePar.p)
    for idx in indices:
        pars[idx] = "params[i]"
    pars = "[" + ",".join([str(p) for p in pars]) + "]"

    # generate prob_func which remakes the ODE problem
    # for each different parameter set
    Main.params = values
    Main.eval(
        f"""
    @everywhere params = $params
    @everywhere function prob_func(prob,i,repeat)
        remake(prob, p = {pars})
    end
    """
    )


def setup_parameter_scan_zipped(odePar, parameters, values):
    """
    Convenience function for initializing a 1D parameter scan over
    multiple parameters, with each parameter scanning over a different
    set of parameters.

    Args:
        odePar (odeParameters): object containing all the parameters
                                for the OBE system.
        parameters (list): list of parameters to scan over
        values (list, np.ndarray): list/array of values to scan over.
    """
    # get the indices of each parameter that is scanned over,
    # as defined in odePars. If a parameter is not scanned over,
    # use the variable definition
    pars = list(odePar.p)
    for idN, parameter in enumerate(parameters):
        if isinstance(parameter, (list, tuple)):
            indices = [odePar.get_index_parameter(par) for par in parameter]
        else:
            indices = [odePar.get_index_parameter(parameter)]
        for idx in indices:
            pars[idx] = f"params[i,{idN+1}]"
    params = np.array(list(zip(*values)))

    pars = "[" + ",".join([str(p) for p in pars]) + "]"

    # generate prob_func which remakes the ODE problem for
    # each different parameter set
    Main.params = params
    Main.eval(
        f"""
    @everywhere params = $params
    @everywhere function prob_func(prob, i, repeat)
        remake(prob, p = {pars})
    end
    """
    )


def setup_parameter_scan_ND(odePar, parameters, values, randomize=False):
    """
    Convenience function for generating an ND parameter scan.
    For each parameter a list or np.ndarray of values is supplied,
    and each possible combination between all parameters is simulated.

    Args:
        odePar (odeParameters): object containing all the parameters for
                                the OBE system.
        parameters (list, np.ndarray): strs of parameters to scan over.
        values (list, np.ndarray): list or np.ndarray of values to scan over
                                    for each parameter
        randomize (bool, False): randomize the parameter scan if True,
                                might increase performance when some parameter
                                combinations are more computatationally
                                expensive than others.
    Returns:
        np.ndarray: return the indices that randomized the parameter
                    combinations, if randomize=True
    """
    pars = str(odePar.p)[1:-1].split(",")

    for idN, parameter in enumerate(parameters):
        if isinstance(parameter, (list, tuple)):
            indices = [odePar.get_index_parameter(par) for par in parameter]
        else:
            indices = [odePar.get_index_parameter(parameter)]
        for idx in indices:
            pars[idx] = f"params[i,{idN+1}]"
    pars = "[" + ",".join(pars) + "]"
    # create all possible combinations between parameter values with
    # meshgrid
    params = np.array(np.meshgrid(*values)).T.reshape(-1, len(values))
    Main.params = params
    if randomize:
        # randomize value sequence if randomize=True
        ind_random = np.random.permutation(len(params))
        Main.params = params[ind_random]

    # generate the prob_func that remakes the ODE problem for each
    # new parameter value set
    Main.eval(
        f"""
    @everywhere params = $params
    @everywhere function prob_func(prob, i, repeat)
        remake(prob, p = {pars})
    end
    """
    )
    if randomize:
        return ind_random


def setup_ratio_calculation(states, output_func=None):
    if output_func is None:
        output_func = "output_func"
    cmd = ""
    if isinstance(states[0], (list, np.ndarray, tuple)):
        for state in states:
            cmd += f"sum(real(diag(sol.u[end])[{state}]))/sum(real(diag(sol.u[1])[{state}])), "
        cmd = cmd.strip(", ")
        cmd = "[" + cmd + "]"
    else:
        cmd = (
            f"sum(real(diag(sol.u[end])[{states}]))/sum(real(diag(sol.u[1])[{states}]))"
        )

    Main.eval(
        f"""
    @everywhere function {output_func}(sol,i)
        if size(sol.u)[1] == 1
            return NaN, false
        else
            val = {cmd}
            return val, false
        end
    end"""
    )
    return output_func


def handle_randomized_ensemble_solution(ind_random, sol_name="sim"):
    order_restored = np.argsort(ind_random)
    # Main.params = Main.params[order_restored]
    return (
        np.asarray(Main.params[order_restored]),
        np.array(Main.eval(f"{sol_name}.u"))[order_restored],
    )


def setup_initial_condition_scan(values):
    Main.params = values
    Main.eval("@everywhere params = $params")
    Main.eval(
        """
    @everywhere function prob_func(prob,i,repeat)
        remake(prob,u0=params[i])
    end
    """
    )


def setup_state_integral_calculation(states, output_func=None, nphotons=False, Γ=None):
    """Setup an integration output_function for an EnsembleProblem. 
    Uses trapezoidal integration to integrate the states.
    
    Args:
        states (list): list of state indices to integrate
        nphotons (bool, optional): flag to calculate the number of photons, 
                                    e.g. normalize with Γ
        Γ (float, optional): decay rate in 2π Hz (rad/s), not necessary if already
                                loaded into Julia globals
    """
    if output_func is None:
        output_func = "output_func"
    if nphotons & Main.eval("@isdefined Γ"):
        Main.eval(
            f"""
        @everywhere function {output_func}(sol,i)
            return Γ.*trapz(sol.t, [real(sum(diag(sol.u[j])[{states}])) for j in 1:size(sol)[3]]), false
        end"""
        )
    else:
        if nphotons:
            assert (
                Γ is not None
            ), "Γ not defined as a global in Julia and not supplied to function"
            Main.eval(f"@everywhere Γ = {Γ}")
            Main.eval(
                f"""
            @everywhere function {output_func}(sol,i)
                return {Γ}.*trapz(sol.t, [real(sum(diag(sol.u[j])[{states}])) for j in 1:size(sol)[3]]), false
            end"""
            )
        else:
            Main.eval(
                f"""
            @everywhere function {output_func}(sol,i)
                return trapz(sol.t, [real(sum(diag(sol.u[j])[{states}])) for j in 1:size(sol)[3]]), false
            end"""
            )
    return output_func


def setup_state_integral_map(
    states, Δx, xmax, vx_index, integral_steps=10, output_func=None
):
    if output_func is None:
        output_func = "output_func"
    Main.eval(
        f"""
    @everywhere function {output_func}(sol, i)
        results = zeros({int(xmax/Δx)},2)
        vx = sol.prob.p[{vx_index}]
        Δt = {Δx}/vx
        tarray = collect(range(0,stop=Δt, length = {integral_steps}))
        @inbounds begin
            for i=1:{int(xmax/Δx)}
                t = Δt*(i-1) + Δt/2
                results[i,1] = t*vx
                results[i,2] = Γ*trapz(tarray .+ (t-Δt/2), 
                                    [sum(diag(real(sol(ti + (t-Δt/2))))[{states}]) for ti in tarray]
                                )
            end
        end
        return results, false
    end
    """
    )

    return output_func


def setup_discrete_callback_terminate(
    odepars: odeParameters, stop_expression: str, callback_name=None
):
    # parse expression string to sympy equation
    expression = smp.parsing.sympy_parser.parse_expr(stop_expression)
    # extract symbols in expression and convert to a list of strings
    symbols_in_expression = list(expression.free_symbols)
    symbols_in_expression = [str(sym) for sym in symbols_in_expression]
    # check if all symbols are parameters of the ODE
    odepars.check_symbols_in_parameters(symbols_in_expression)

    # remove t
    symbols_in_expression.remove("t")
    # get indices of symbols
    indices = [
        odepars.get_index_parameter(sym, mode="julia") for sym in symbols_in_expression
    ]
    for idx, sym in zip(indices, symbols_in_expression):
        stop_expression = stop_expression.replace(str(sym), f"integrator.p[{idx}]")
    if callback_name is None:
        callback_name = "cb"
    Main.eval(
        f"""
        @everywhere condition(u,t,integrator) = {stop_expression}
        @everywhere affect!(integrator) = terminate!(integrator)
        {callback_name} = DiscreteCallback(condition, affect!)
    """
    )
    return callback_name


def setup_problem(
    odepars: odeParameters, tspan: list, ρ: np.ndarray, problem_name="prob"
):
    odepars.generate_p_julia()
    Main.ρ = ρ
    Main.tspan = tspan
    assert Main.eval(
        "@isdefined Lindblad_rhs!"
    ), "Lindblad function is not defined in Julia"
    Main.eval(
        f"""
        {problem_name} = ODEProblem(Lindblad_rhs!,ρ,tspan,p)
    """
    )


def setup_problem_parameter_scan(
    odepars: odeParameters,
    tspan: list,
    ρ: np.ndarray,
    parameters: list,
    values: np.ndarray,
    dimensions: int = 1,
    problem_name="prob",
    output_func=None,
    zipped=False,
):
    setup_problem(odepars, tspan, ρ, problem_name)
    if dimensions == 1:
        if zipped:
            setup_parameter_scan_zipped(odepars, parameters, values)
        else:
            setup_parameter_scan_1D(odepars, parameters, values)
    else:
        setup_parameter_scan_ND(odepars, parameters, values)
    if output_func is not None:
        Main.eval(
            f"""
            ens_{problem_name} = EnsembleProblem({problem_name}, 
                                                    prob_func = prob_func,
                                                    output_func = {output_func}
                                                )
        """
        )
    else:
        Main.eval(
            f"""
            ens_{problem_name} = EnsembleProblem({problem_name}, 
                                                    prob_func = prob_func)
        """
        )
    return f"ens_{problem_name}"


def solve_problem(
    method="Tsit5()",
    abstol=1e-7,
    reltol=1e-4,
    dt=1e-8,
    callback=None,
    problem_name="prob",
    progress=False,
    saveat=None,
    dtmin=None,
    maxiters=None,
):
    if saveat is None:
        saveat = "[]"
    if dtmin is None:
        dtmin = "nothing"
        force_dtmin = "false"
    else:
        force_dtmin = "true"
    if maxiters is None:
        maxiters = 1e5
    if callback is not None:
        Main.eval(
            f"""
            sol = solve({problem_name}, {method}, abstol = {abstol}, 
                        reltol = {reltol}, dt = {dt}, 
                        progress = {str(progress).lower()}, 
                        callback = {callback}, saveat = {saveat}, 
                        dtmin = {dtmin}, maxiters = {maxiters}, 
                        force_dtmin = {force_dtmin}
                    )
        """
        )
    else:
        Main.eval(
            f"""
            sol = solve({problem_name}, {method}, abstol = {abstol}, 
                        reltol = {reltol}, dt = {dt},
                        progress = {str(progress).lower()}, saveat = {saveat},
                        dtmin = {dtmin}, maxiters = {maxiters},
                        force_dtmin = {force_dtmin}
                    )
        """
        )


def solve_problem_parameter_scan(
    method="Tsit5()",
    distributed_method="EnsembleDistributed()",
    abstol=1e-7,
    reltol=1e-4,
    dt=1e-8,
    save_everystep=True,
    callback=None,
    ensemble_problem_name="ens_prob",
    trajectories=None,
    saveat=None,
):
    if trajectories is None:
        trajectories = "size(params)[1]"
    if saveat is None:
        saveat = "[]"
    if callback is not None:
        Main.eval(
            f"""
            sol = solve({ensemble_problem_name}, {method}, {distributed_method}, 
                        abstol = {abstol}, reltol = {reltol}, dt = {dt},
                        trajectories = {trajectories}, callback = {callback},
                        save_everystep = {str(save_everystep).lower()},
                        saveat = {saveat}
                    )
        """
        )
    else:
        Main.eval(
            f"""
            sol = solve({ensemble_problem_name}, {method}, {distributed_method}, 
                        abstol = {abstol}, reltol = {reltol}, dt = {dt},
                        trajectories = {trajectories},
                        save_everystep = {str(save_everystep).lower()},
                        saveat = {saveat}
                    )
        """
        )


def get_results_parameter_scan(scan_values=None):
    results = np.array(Main.eval("sol.u"))
    if scan_values is not None:
        if isinstance(scan_values, list) or (scan_values.ndim > 1):
            results = results.reshape([len(val) for val in scan_values])
            X, Y = np.meshgrid(*scan_values)
            return X, Y, results.T
        else:
            return scan_values, results
    return results


def get_results():
    """Retrieve the results of a single trajectory OBE simulation solution.

    Returns:
        tuple: tuple containing the timestamps and an n x m numpy arra, where
                n is the number of states, and m the number of timesteps
    """
    results = np.real(np.einsum("jji->ji", np.array(Main.eval("sol[:]")).T))
    t = Main.eval("sol.t")
    return t, results


def do_simulation_single(
    odepars,
    tspan,
    ρ,
    terminate_expression=None,
    dt=1e-8,
    saveat=None,
    dtmin=None,
    maxiters=None,
):
    """Perform a single trajectory solve of the OBE equations for a specified 
    TlF system.

    Args:
        odepars (odeParameters): object containing the ODE parameters used in
        the solver
        tspan (list, tuple): time range to solve for
        terminate_expression (str, optional): Expression that determines when to 
                                            stop integration. Defaults to None.
        saveat (array or float, optional): save solution at timesteps given by 
                                            saveat, either a list or every 
                                            saveat
        dtmin (float, optional): minimum dt allowed for adaptive timestepping
        maxiters (float, optional): maximum number of steps allowed

    Returns:
        tuple: tuple containing the timestamps and an n x m numpy array, where
                n is the number of states, and m the number of timesteps
    """
    callback_flag = False
    if terminate_expression is not None:
        setup_discrete_callback_terminate(odepars, terminate_expression)
        callback_flag = True
    setup_problem(odepars, tspan, ρ)
    if callback_flag:
        solve_problem(
            callback="cb", saveat=saveat, dtmin=dtmin, dt=dt, maxiters=maxiters
        )
    else:
        solve_problem(saveat=saveat, dtmin=dtmin, dt=dt, maxiters=maxiters)
    return get_results()


def get_indices_diag_flattened(n):
    return np.diag(np.arange(0, n * n).reshape(-1, n))
