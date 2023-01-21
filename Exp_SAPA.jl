include("LLR.jl")
include("Approx_Proj.jl")
include("ADMM_LP_SAPA.jl")

###
# This code generates multiple Monte-Carlo simulations to estimate the frame/bit
# error rate and iteration complexity of ADMM-LP decoding with different parity
# polytope projection algorithms. To estimate the error rate, we require sufficient
# error patterns to average over. To do so, we allocate a parameter called
# "total_error_sim" to identify how many error patterns should be observed to obtain
# a quite precise estimation of error rates. The default value is set to 100.
###

# Choose your LDPC code in the following by uncommenting the corresponding lines
# However, make sure to change N_Dec and the arguments for different experiments accordingly
# Uncomment to use the Margulis code:
include("LDPC Codes/Margulis/getH_Margulis_ADMM.jl")
x_inp = readdlm("LDPC Codes/Margulis/codeword2640.txt")[1,:] # input codeword

# Uncomment to use the Quasi-cyclic code:
# include("LDPC Codes/QuasiCyclic/getH_QC_ADMM.jl")
# x_inp = zeros(N)

# Uncomment to use the WiGig code:
# include("LDPC Codes/WiGig/getH_WiGig_ADMM.jl")
# x_inp = readdlm("LDPC Codes/WiGig/wifiCodewords.txt")[2,:]

# Uncomment to use the IEEE [576, 288] code:
# include("LDPC Codes/IEEE_576x288/getH_576x288_ADMM.jl")
# x_inp = readdlm("LDPC Codes/MIEEE_576x288/IEEE_576x288_codeword.txt")[1,:]


# decoder parameters
μ = 3
ρ = 1.9
ϵ = 10^(-5)
method = "L2"
α = 0.8
N_iter = 20
N_ET = 1
γ = zeros(N)

# experiment parameters
SNR_vec = 1:0.2:6
N_SNR = length(SNR_vec)
N_Dec = 8 # 1) EXACT, 2) EVA, 3) LSA, 4) 2-SAPA, 5) 3-SAPA, 6) 4-SAPA, 7) 5-SAPA 8) APA(which is equal to 6-SAPA for Margulis with d = 8)
BER_vec = zeros(N_SNR, N_Dec)
FER_vec = zeros(N_SNR, N_Dec)
Iter_vec = zeros(N_SNR, N_Dec)
Nsim_vec = zeros(N_SNR, N_Dec)

total_error_sim = 100 # number of error patters to get per SNR

for j = 1:N_SNR

    total_sim = zeros(N_Dec)
    correct_sim = zeros(N_Dec)
    error_sim = zeros(N_Dec)
    bit_error = zeros(N_Dec)
    correct_iter = zeros(N_Dec)

    SNR = SNR_vec[j]
    println("SNR = $SNR:")

    while minimum(error_sim) < total_error_sim

        UpdateLLR(x_inp, γ, "AWGN", SNR)

        # EXACT
        if error_sim[1] < total_error_sim
            x_dec, iter = ADMMdecoder(method, F_list, γ, N_iter, ϵ, μ, ρ, α, N_ET, 0, 6)
            if sum( Int.( x_dec .!= x_inp ) ) != 0
                error_sim[1] += 1
                bit_error[1] += sum( Int.( x_dec .!= x_inp ) )
            else
                correct_sim[1] += 1
                correct_iter[1] += iter
            end
            total_sim[1] += 1
        end


        # EVA
        if error_sim[2] < total_error_sim
            x_dec, iter = ADMMdecoder(method, F_list, γ, N_iter, ϵ, μ, ρ, α, N_ET, 3, 6)
            if sum( Int.( x_dec .!= x_inp ) ) != 0
                error_sim[2] += 1
                bit_error[2] += sum( Int.( x_dec .!= x_inp ) )
            else
                correct_sim[2] += 1
                correct_iter[2] += iter
            end
            total_sim[2] += 1
        end


        # LSA
        if error_sim[3] < total_error_sim
            x_dec, iter = ADMMdecoder(method, F_list, γ, N_iter, ϵ, μ, ρ, α, N_ET, 2, 6)
            if sum( Int.( x_dec .!= x_inp ) ) != 0
                error_sim[3] += 1
                bit_error[3] += sum( Int.( x_dec .!= x_inp ) )
            else
                correct_sim[3] += 1
                correct_iter[3] += iter
            end
            total_sim[3] += 1
        end


        # 2-SAPA
        if error_sim[4] < total_error_sim
            x_dec, iter = ADMMdecoder(method, F_list, γ, N_iter, ϵ, μ, ρ, α, N_ET, 4, 2)
            if sum( Int.( x_dec .!= x_inp ) ) != 0
                error_sim[4] += 1
                bit_error[4] += sum( Int.( x_dec .!= x_inp ) )
            else
                correct_sim[4] += 1
                correct_iter[4] += iter
            end
            total_sim[4] += 1
        end


        # 3-SAPA
        if error_sim[5] < total_error_sim
            x_dec, iter = ADMMdecoder(method, F_list, γ, N_iter, ϵ, μ, ρ, α, N_ET, 4, 3)
            if sum( Int.( x_dec .!= x_inp ) ) != 0
                error_sim[5] += 1
                bit_error[5] += sum( Int.( x_dec .!= x_inp ) )
            else
                correct_sim[5] += 1
                correct_iter[5] += iter
            end
            total_sim[5] += 1
        end


        # 4-SAPA
        if error_sim[6] < total_error_sim
            x_dec, iter = ADMMdecoder(method, F_list, γ, N_iter, ϵ, μ, ρ, α, N_ET, 4, 4)
            if sum( Int.( x_dec .!= x_inp ) ) != 0
                error_sim[6] += 1
                bit_error[6] += sum( Int.( x_dec .!= x_inp ) )
            else
                correct_sim[6] += 1
                correct_iter[6] += iter
            end
            total_sim[6] += 1
        end


        # 5-SAPA
        if error_sim[7] < total_error_sim
            x_dec, iter = ADMMdecoder(method, F_list, γ, N_iter, ϵ, μ, ρ, α, N_ET, 4, 5)
            if sum( Int.( x_dec .!= x_inp ) ) != 0
                error_sim[7] += 1
                bit_error[7] += sum( Int.( x_dec .!= x_inp ) )
            else
                correct_sim[7] += 1
                correct_iter[7] += iter
            end
            total_sim[7] += 1
        end


        # APA
        if error_sim[8] < total_error_sim
            x_dec, iter = ADMMdecoder(method, F_list, γ, N_iter, ϵ, μ, ρ, α, N_ET, 1, 6)
            if sum( Int.( x_dec .!= x_inp ) ) != 0
                error_sim[8] += 1
                bit_error[8] += sum( Int.( x_dec .!= x_inp ) )
            else
                correct_sim[8] += 1
                correct_iter[8] += iter
            end
            total_sim[8] += 1
        end


        if maximum(total_sim) % 10000 == 0 # save the results into .txt file each 10000 simulations
            for k = 1:8
                FER_vec[j, k] = error_sim[k] / total_sim[k]
                BER_vec[j, k] = bit_error[k] / total_sim[k] / N
                Iter_vec[j, k] = correct_iter[k] / correct_sim[k]
                Nsim_vec[j, k] = total_sim[k]
            end
            println("SNR = ", SNR, ", # errors = ", error_sim, ", FER = ", FER_vec[j,:])
            open("SNR_FER_BER_CoIter_Nsim_Margulis_ADMM_Niter20_Nerr100.txt", "w") do io
                writedlm(io, [SNR_vec FER_vec BER_vec Iter_vec Nsim_vec], ',')
            end
        end
    end

    for k = 1:8
        FER_vec[j, k] = error_sim[k] / total_sim[k]
        BER_vec[j, k] = bit_error[k] / total_sim[k] / N
        Iter_vec[j, k] = correct_iter[k] / correct_sim[k]
        Nsim_vec[j, k] = total_sim[k]
    end
    println("SNR = ", SNR, ", # errors = ", error_sim, ", FER = ", FER_vec[j,:])
    open("SNR_FER_BER_CoIter_Nsim_Margulis_ADMM_Niter20_Nerr100.txt", "w") do io
        writedlm(io, [SNR_vec FER_vec BER_vec Iter_vec Nsim_vec], ',')
    end

end
