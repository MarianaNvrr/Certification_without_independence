
include("main.jl")

"""
Run the certification protocol for a Hamiltonian observable `H`.
The protocol estimates the average energy of the use states
using only the randomly selected test states.

Input
- `N`: total number of experimental rounds.
- `H`: Hamiltonian observable, assumed to be an n-qubit Hermitian matrix.
- (optional) `p`: probability of testing a state.
- (optional) `δ`: failure probability. 
- (optional) `verbose`: if true, print protocol diagnostics and results.

Output
- `empirical_deviation_bound`: certification error (see Theorem 3).
- `avg_tested_estimator`: estimator inferred from the test rounds.
- `true_avg_use_gs_energy`: exact average energy of the use states,
  computed only for numerical validation.
"""
function hamiltonian_certification(N,H; p=0.7, δ=0.05,verbose=true)
    d = size(H, 1) ; n = Int(log2(d))    
    N = Int(N)

    # We create a combination of all the indices
    relevant_indices = pauli_indices(n)

    # We calculate the coefficients excluding the identity
    αP = Float64[]
    for ind in relevant_indices
        P = pauli(collect(ind)) 
        push!(αP, real(tr(H * P)) / d)
    end

    L1 = sum(abs.(αP))
    πP = pweights(abs.(αP)/(L1) )

    # Martingale Bounds   
    Wmin = eigen_sparse(H)[1][1][1]
    Wmax = eigen_sparse(H)[1][2][1]
    A = min(Wmin, -(1-p)/p * L1)
    B = max(Wmax,(1-p)/p * L1 ) 
    Δp = B - A

    S_U_omega = Float64[] # True avg energy of use states
    S_T_X     = Float64[] # Estimator values of test states

    ϵ = Δp*sqrt(log(2/δ)/(2N))

    if verbose
        println("--- PRIOR BOUND SETUP ---")
        println("Total rounds N=$N required for target ϵ=$ϵ \n")
    end

    # Source model. We need to assume a state for the simulations
    gs = ground_state(H)
    trajectory = generate_drift_trajectory(N, gs)

    # MAIN LOOP
    loop_time = @elapsed begin
        for t in 1:N
            # Sequence of state at round t
            ρ_t = trajectory[t]

            # Coin flip
            c_t = rand() < p ? 1 : 0
            
            if c_t == 1 #-----------  Test states
                # We sample from the distribution πP
                sample_idx = sample(1:length(relevant_indices), πP)
                Pt_ind = collect(relevant_indices[sample_idx])
                Pt = pauli(Pt_ind)
                
                p_plus = real(1 + tr(Pt * ρ_t)) / 2
                yt = rand() < p_plus ? 1.0 : -1.0
                
                # Single-shot estimator 
                Xt =  yt * L1 * sign(αP[sample_idx])
                push!(S_T_X, Xt)
            else #----------- Use states
                true_gs_energy = real(tr(H * ρ_t))
                push!(S_U_omega, true_gs_energy)
            end
        end
    end

    num_use = length(S_U_omega)
    num_tested = length(S_T_X)

    avg_tested_estimator = (1-p)/p *(sum(S_T_X)/num_use)
    empirical_deviation_bound = (N * Δp / num_use) * sqrt(log(2/δ) / (2N))
    
    certified_lower_bound = avg_tested_estimator - empirical_deviation_bound
    certified_upper_bound = avg_tested_estimator + empirical_deviation_bound
    true_avg_use_gs_energy = sum(S_U_omega) / num_use

    if verbose
        println("--- PROTOCOL RESULTS ---")
        println("States Tested (|S_T| $(p*100)%) : $num_tested")
        println("States Used (|S_U| $((1-p)*100)%):   $num_use")
        println("The certification took $(round(loop_time, digits=2)) seconds.")
        println("-------------------------------------------------")
        println("Certified Upper Bound             :    ≤ $(round(certified_upper_bound, digits=5))")
        println("TRUE Avg energy of use states     :    = $(round(true_avg_use_gs_energy, digits=5))")
        println("Certified Lower Bound             :    ≥ $(round(certified_lower_bound, digits=5))")
        println("Empirical Confidence Width        :    ±$(round(empirical_deviation_bound, digits=5))")
        println("TRUE Avg energy of use states     :    = $(round(true_avg_use_gs_energy, digits=5))")
        println("_________________________________________________")

        if true_avg_use_gs_energy>= certified_lower_bound
            println("\n SUCCESS: Entanglement bounds held correctly. \n")
        else
            println("\n FAILURE: True entanglement fell below the bound. \n")
        end
    end 
   
    return empirical_deviation_bound,avg_tested_estimator, true_avg_use_gs_energy
end

