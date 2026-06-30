using JLD2
using Plots


############################# HAMILTONIAN ############################# 

# Parameters used for the plots
# g = 0.5; J=1; Δ= 1
# N = 1e6

residuals = []
real_deviation = []
p_values = [0.2,0.5]

H_type="H_XXZ" #"H_TFIM"  # select Hamiltonian to be plotted
for p_val in p_values
    filepath = "data/Hamiltonian_$(H_type)_p$(p_val)_n3_N1e6_coupling1.jld2"
    data = load(filepath)

    push!(residuals, data["residuals"])
    push!(real_deviation, data["real_deviation"])
end

default(fontfamily="Computer Modern", grid=false, framestyle=:box)

plt= plot(
    xlabel="\$ h_{T} -  h_{U} \$",
    ylabel="Frequency (1000 trajectories)",
    framestyle = :box,
    tickfontsize = 10,
    guidefontsize = 12,
    legendfontsize = 10,
    legend=:topleft,
    fontfamily = "Computer Modern"
)


histogram!(plt, residuals[2], 
        bins=30, 
        normalize=:probability,
        linecolor=:black,
        alpha=0.9,           # CRÍTICO: Transparencia para ver los gráficos de atrás
        linewidth=1,         # Quitar los bordes hace que se vea mucho más limpio superpuesto
        color=:RoyalBlue,
        label="\$ p = $(p_values[2]) \$"
    )

histogram!(plt, residuals[1], 
        bins=40, 
        normalize=:probability,
        linecolor=:black,
        alpha=0.7,           
        linewidth=1,         
        color=:salmon,
        label="\$ p = $(p_values[1]) \$"
    )

colorss = [:salmon,:RoyalBlue]
for i in [2,1]
    vline!(plt, [mean(real_deviation[i])], color=colorss[i], linewidth=2, linestyle=:dash,alpha = 0.8, label=false)
    vline!(plt, [-mean(real_deviation[i])], color=colorss[i], linewidth=2, linestyle=:dash,alpha = 0.8, label=false)
end

display(plt)


############################# WITNESS ############################# 


N_fixed = 1e5

# Parameters used for the plots
# p_val = 0.026 
# ϵ_target = 0.03 
# st = state("weak_entangled", 2; α=3/4)
# W = 3*I(4)/4 - st

filepath = "data/witness_partialEntg_N1e$(Int(log10(N_fixed)))_drift.001.jld2"
data = load(filepath)
residuals_cert = data["residuals_cert"]
epsilons_cert = data["epsilons_cert"]
residuals_verif = data["residuals_verif"]

# epsilon certification 
mean_eps_cert = mean(epsilons_cert)

plt = plot(xlabel="\$  w_{U} - w_{T} \$",
    ylabel="Frequency (1000 trajectories)",
    normalize=:probability,
    tickfontsize=10,
    guidefontsize=12,
    legendfontsize=10,
    legend=:topright
)


histogram!(plt, residuals_verif, 
    bins=20, 
    normalize=:probability,
    linecolor=:black,
    alpha=0.7, 
    linewidth=1, 
    color=:salmon,
    label="Static certification"
)

histogram!(plt, residuals_cert, 
    bins=15, 
    normalize=:probability,
    linecolor=:black,
    alpha=0.9, 
    linewidth=1, 
    color=:RoyalBlue,
    label="Dynamic certification (this work)"
)

display(plt)


############################# IBM EXPERIMENT ############################# 
using NPZ

data = npzread("programs/data/fidelity_p0.1.npz")

# Parameters used for the plots
# p = 0.1
F_est = data["F_est"]
F_exp = data["F_exp"]
epsilon = data["epsilon"]

residuals = F_est - F_exp

plt= plot(xlabel="\$ F_{T} - F_{U} \$",
    ylabel="Frequency (64 trajectories)",
    bins = 30,
    tickfontsize = 10,
    guidefontsize = 12,
    legendfontsize = 10,
    legend=false
)

histogram!(plt, residuals, 
        bins=20, 
        normalize=:probability,
        linecolor=:black,
        alpha=0.9,           
        linewidth=1,         
        color=:RoyalBlue,
        ylims=(0,0.13)
    )

real_deviation = mean(epsilon)
vline!(plt, [real_deviation], color=:RoyalBlue, linewidth=2, linestyle=:dash,alpha = 0.8, label=false)
vline!(plt, [-real_deviation], color=:RoyalBlue, linewidth=2, linestyle=:dash,alpha = 0.8, label=false)
