from julia import Main

__all__ = ["solve_problem_parameter_scan_progress"]


def solve_problem_parameter_scan_progress(
    method="Tsit5()",
    distributed_method="EnsembleDistributed()",
    abstol=1e-7,
    reltol=1e-4,
    save_everystep=True,
    callback=None,
    problem_name="prob",
    ensemble_problem_name="ens_prob",
    trajectories=None,
    output_func=None,
    saveat=None,
):
    if trajectories is None:
        trajectories = "size(params)[1]"
    if callback is None:
        callback = "nothing"
    if saveat is None:
        saveat = "[]"
    if output_func is None:
        Main.eval(
            """
            @everywhere function output_func_progress(sol, i)
                put!(channel, 1)
                sol, false
            end
        """
        )
    else:
        Main.eval(
            f"""
            @everywhere function output_func_progress(sol, i)
                put!(channel, 1)
                a,b = {output_func}(sol, i)
                return a,b
            end
        """
        )
    Main.eval(
        f"""
        {ensemble_problem_name} = EnsembleProblem({problem_name}, 
                                                prob_func = prob_func,
                                                output_func = output_func_progress
                                            )
    """
    )

    Main.eval(
        """
        if !@isdefined channel
            const channel = RemoteChannel(()->Channel{Int}(1))
            @everywhere const channel = $channel
        end
    """
    )

    Main.eval(
        f"""
        progress = Progress({trajectories}, showspeed = true)
        @sync sol = begin
            @async begin
                tasksdone = 0
                while tasksdone < {trajectories}
                    tasksdone += take!(channel)
                    update!(progress, tasksdone)
                end
            end
            @async begin
                @time global sol = solve({ensemble_problem_name}, {method}, 
                            {distributed_method}, trajectories={trajectories},
                            abstol = {abstol}, reltol = {reltol}, 
                            callback = {callback}, 
                            save_everystep = {str(save_everystep).lower()},
                            saveat = {saveat})
            end
    end
    """
    )
