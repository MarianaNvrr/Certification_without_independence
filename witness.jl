include("main.jl")

"""
Run the certification protocol for an observable `W`.
The protocol estimates the average witness value of `W` on the use states
using only the randomly selected test states.

Input
- `st`: target state used to generate the simulated drifting trajectory.
- `N`: total number of experimental rounds.
- `W`: observable to certify, assumed to be an n-qubit Hermitian matrix.
- (optional) `p`: probability of testing a state.
- (optional) `δ`: confidence value.
- (optional) `verbose`: if true, print protocol diagnostics and results.

Output
- `empirical_deviation_bound`: certification error (see Theorem 3).
- `avg_tested_estimator`: estimator inferred from the test rounds.
- `true_avg_use_val`: exact average value of the use states,
  computed only for numerical validation.
"""
function witness_certified(st,N,W; p=0.7, δ=0.05, verbose=true)
    d = size(W, 1); n = Int(log2(d))
    N = Int(N)

    relevant_indices = pauli_indices(n)

    αP = Float64[]; αI = tr(W)/d
    for ind in relevant_indices
        P = pauli(collect(ind)) 
        push!(αP, real(tr(W * P))/d)
    end

    L1 = sum(abs.(αP)) 
    πP = pweights(abs.(αP)/(L1) )

    # Martingale Bounds 
    W = (W-αI*I(d))  
    Wmin = minimum(real.(eigvals(W)))
    Wmax = maximum(real.(eigvals(W)))
    A = min(Wmin, -(1-p)/p * L1) 
    B = max(Wmax,(1-p)/p * L1 )
    Δp = B - A

    S_U_rho = Float64[] 
    S_T_X   = Float64[] 

    ϵ = Δp*sqrt(log(2/δ)/(2N))

    # Source model. We need to assume a state for the simulations
    trajectory = generate_drift_trajectory(N, st)

    if verbose
        println("--- PRIOR BOUND SETUP ---")
        println("Total rounds N=$N required for target ϵ=$ϵ \n")
    end

    # MAIN LOOP
    loop_time = @elapsed begin
        for t in 1:N
            ρ_t  = trajectory[t]

            # Coin Flip
            c_t = rand() < p ? 1 : 0
            
            if c_t == 1 #-----------  Test states
                # We sample from the distribution πP
                sample_idx = sample(1:length(relevant_indices), πP)
                Pt_ind = collect(relevant_indices[sample_idx])
                Pt = pauli(Pt_ind)
                
                p_plus = real(1 + tr(Pt * rho_t)) / 2
                yt = rand() < p_plus ? 1.0 : -1.0
                
                # Single-shot estimator for fidelity
                Xt =  αI +  yt * L1 * sign(αP[sample_idx])
                push!(S_T_X, Xt)
            else #----------- Use states
                true_val = real(tr(W * rho_t))
                push!(S_U_rho, true_val)
            end
        
        end
    end

    num_use = length(S_U_rho)
    num_tested = length(S_T_X)

    avg_tested_estimator =(1-p)/p *(sum(S_T_X)/num_use)
    empirical_deviation_bound = (N * Δp / num_use) * sqrt(log(2/δ) / (2N))
    
    certified_lower_bound = avg_tested_estimator - empirical_deviation_bound
    certified_upper_bound = avg_tested_estimator + empirical_deviation_bound
    true_avg_use_val = sum(S_U_rho) / num_use

    if verbose
        println("--- PROTOCOL RESULTS ---")
        println("  For p = $(p*100)%")
        println("States Tested (|S_T|): $num_tested")
        println("States use (|S_U|):   $num_use")
        println("The certification took $(round(loop_time, digits=2)) seconds.")
        println("-------------------------------------------------")
        println("Certified upper Bound           :  ≤ $(round(certified_upper_bound, digits=5))")
        println("TRUE Avg Witness of use States  :  = $(round(true_avg_use_val, digits=5))")
        println("Certified Lower Bound           :  ≥ $(round(certified_lower_bound, digits=5))")
        println("Empirical Confidence Width      :  ± $(round(empirical_deviation_bound, digits=5))")
        println("_________________________________________________")

        if true_avg_use_val  >= certified_lower_bound
            println("SUCCESS: The true witness is bounded correctly by the protocol.")
        else
            println("FAILURE: The true fidelity fell below the bound (should happen < $δ of the time).")
        end
    end
    return empirical_deviation_bound, avg_tested_estimator, true_avg_use_val
end


"
Computes the sample complexity according to Theorem 2.

Input 
- `trajectory`: sequence of states ρ_t. 
- `W`: observable to estimate. 
- `ϵ`: target additive error. 
- `δ`: failure probability.

Output 
- Estimated average value of `W`. 
- Number of samples used by the verification protocol.
"
function complexity_verification(trajectory,W, ϵ, δ)
    d = size(W, 1); n = Int(log2(d))

    # We create a combination of all the indices
    relevant_indices = pauli_indices(n)
    
    # We calculate the coefficients excluding the identity
    αP = Float64[]
    for ind in relevant_indices
        P = pauli(collect(ind)) 
        push!(αP, real(tr(W * P))/d)
    end

    L1 = sum(abs.(αP) )
    πP = pweights(abs.(αP)/L1)

    # Sample complexity from Theorem 2.
    N = Int(ceil(2 * (L1 / ϵ)^2 * log(2 / δ)))

    X_hat_sum=0.

    # MAIN LOOP
    for t in 1:N
        sample_idx = sample(1:length(relevant_indices), πP)
        Pt_ind = collect(relevant_indices[sample_idx])
        Pt = pauli(Pt_ind)

        ρ_t = trajectory[t]

        p_plus = real(1 + tr(Pt * ρ_t)) / 2
        yt = rand() < p_plus ? 1.0 : -1.0
        
        # Single-shot estimator 
        Xt = yt * L1 * sign(αP[sample_idx])
        X_hat_sum += Xt
    end
    return (real(tr(W))/d) + (X_hat_sum / N), N
end

"
Run a static certification protocol for an observable `W`. 
The protocol first consumes an initial block of states to 
estimate the value of `W` using the verification protocol. 
We also estimate the average witness of the remainin states for
bechmarking.

Input 
- `st`: state needed for the simulations. 
- `N`: total number of rounds. 
- `W`: witness to estimate. 
- (optional) `ϵ`: target additive error for the initial verification block. 
- (optional) `δ`: failure probability for the verification estimate. 
- (optional) `verbose`: if true, print protocol diagnostics and results. 

Output 
- `W_estimated`: estimate obtained from the initial verified block. 
- `avg_true_fid`: exact average value of `W` on the remaining use states, 
computed only for numerical validation.
"
function witness_static(st,N, W; ϵ=0.03, δ=0.05, verbose=true)
    N = Int(N)
    S_U_ρ = Float64[]

    # Source model. We need to assume a state for the simulations
    trajectory = generate_drift_trajectory(N, st)

    W_estimated, N_ver = complexity(trajectory,W, ϵ, δ)
    N_rest = N - N_ver

    if verbose
        println("----- VERIFICATION PROTOCOL RESULTS --------")
        println("Verified states   : $N_ver ($(N_ver*100/N)%)")
        println("Teleported states : $N_rest ($(N_rest*100/N)%)")
    end

    loop_time = @elapsed begin
        for t in 1:N_rest
            ρ_t = trajectory[N_ver + t]

            true_val = real(tr(W * ρ_t))
            push!(S_U_ρ, true_val)
        end
    end

    num_use = length(S_U_ρ)

    avg_true = sum(S_U_ρ)/ num_use

    if verbose
        println("The verification took $(round(loop_time, digits=2)) seconds.")
        println("-------------------------------------------")
        println("Verified witness  : $W_estimated")
        println("True witness      : $avg_true")
        println("_____________________________________________")
    end

    return W_estimated, avg_true
end

