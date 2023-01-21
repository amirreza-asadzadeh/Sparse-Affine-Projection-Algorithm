using LinearAlgebra

function ADMMdecoder(DecoderMethod, F_list, γ, maxIteration, ϵ, μ, ρ, α, N_ET, Proj, χ)

    ###
    # This function applies non-penalized/l2-penalized ADMM-LP decoding.
    #
    # DecoderMethod: "LP" for non-penalized, "L2" for l2-penalized decoder
    # F_list: parity check matrix given in the format of a F_list
    # γ: the Log Likelihood Ratio (LLR) vector of size N
    # maxIteration: the maximum allowed number of iterations
    # ϵ: stopping precision, usually set to 10^(-5)
    # μ: the Lagrangian penalty parameter, usually set to 3
    # ρ: the over-relaxation parameter, generally 1 < ρ < 2, ρ = 1 for non-over-relaxed case, recommended to set to 1.9
    # α: the parameter for l2-penalty term, usually set to 0.8
    # N_ET: the number of iteration to check early termination
    # Proj: determines the type of parity polytope projection. The options are:
    # --- Exact parity polytope projection -------> Proj = 0
    # --- Affine Projeciotn Algorithm ------------> Proj = 1
    # --- Line-Segment Algorithm -----------------> Proj = 2
    # --- Even-Vertex Algorithm ------------------> Proj = 3
    # --- χ-Sparse Affine Projection Algorithm ---> Proj = 4
    # χ: the hyperparameter for χ-SAPA, only effective when Proj is set to 4
    #
    # Returns:
    # x_dec: a binary (0/1) array of size N, the decoded output
    # iter: the number of iterations it took for the decoder to decode
    # feas1_vec: the primal residual vector at convergence (the last iteration)
    # feas2_vec: the dual residual vector at convergence (the last iteration)
    #
    ###

    λ = Array{Array{Float64, 1}, 1}(undef, M)
    zReplica = Array{Array{Float64, 1}, 1}(undef, M)
    zOld = Array{Array{Float64, 1}, 1}(undef, M)
    Pjx = Array{Array{Float64, 1}, 1}(undef, M)

    for j = 1:M
        λ[j] = zeros(CheckDegree[j])
        zReplica[j] = ones(CheckDegree[j]).*0.5
        zOld[j] = ones(CheckDegree[j]).*0.5
        Pjx[j] = ones(CheckDegree[j]).*0.5
    end

    x_dec = ones(N).*0.5
    x_hard = similar(x_dec)
    temp = zeros(N)

    iter = 1;
    feas_1 = 1e3; # holds for sum over j of abs(P_jx^k - z_j^k)
    feas_2 = 1e3; # holds for sum over j of abs(z_j^k - z_j^(k-1))
    feas_tol = ϵ^2 * M * maximum(CheckDegree);

    while( iter <= maxIteration && (feas_1 >= feas_tol || feas_2 >= feas_tol) )

        # x-Update
        fill!(temp, 0.0)
        for j = 1:M
            for k = 1:CheckDegree[j]
                IDX = F_list[j][k]
                temp[IDX] += zReplica[j][k] - λ[j][k] ./ μ
            end
        end
        temp .-= γ ./ μ

        #select decoder
        if occursin("LP", DecoderMethod)
            for i = 1:N
                x_dec[i] = min( max( ( temp[i] ) / (VariableDegree[i]), 0.0), 1.0)
            end
        elseif occursin("L2", DecoderMethod)
            for i = 1:N
                x_dec[i] = min( max( ( temp[i] - (α / μ) ) / ( VariableDegree[i] - ( 2 * α / μ ) ), 0.0), 1.0)
            end
        end


        # z-Update, lambda-Update, feasibility check
        feas_1 = 0
        feas_2 = 0
        for j = 1:M
            zOld[j] .= zReplica[j]
            for k = 1:CheckDegree[j]
                IDX = F_list[j][k]
                Pjx[j][k] = x_dec[IDX]
            end

            if Proj == 0
                zReplica[j] = EXACT( ( ρ .* Pjx[j] ) .+ ( ( 1 - ρ ) .* zOld[j] ) .+ ( λ[j]./ μ ) )
            elseif Proj == 1
                zReplica[j] = APA( ( ρ .* Pjx[j] ) .+ ( ( 1 - ρ ) .* zOld[j] ) .+ ( λ[j]./ μ ) )
            elseif Proj == 2
                zReplica[j] = LSA( ( ρ .* Pjx[j] ) .+ ( ( 1 - ρ ) .* zOld[j] ) .+ ( λ[j]./ μ ) )
            elseif Proj == 3
                zReplica[j] = EVA( ( ρ .* Pjx[j] ) .+ ( ( 1 - ρ ) .* zOld[j] ) .+ ( λ[j]./ μ ) )
            elseif Proj == 4
                zReplica[j] =χ_SAPA( ( ρ .* Pjx[j] ) .+ ( ( 1 - ρ ) .* zOld[j] ) .+ ( λ[j]./ μ ), χ )
            end

            difference1 = ( ρ .* Pjx[j] ) .+ ( ( 1 - ρ ) .* zOld[j] ) .- zReplica[j]
            λ[j] .+= μ .* ( difference1 )
            feas_1 += norm( difference1 )^2
            difference2 = zReplica[j] - zOld[j]
            feas_2 += norm( difference2 )^2
        end

        # Early Termination
        if iter % N_ET == 0
            HardDecode(x_dec, x_hard)
            sum_syn = 0
            for j = 1:M
                sum_syn += sum( x_hard[ F_list[j] ] ) % 2
            end
            if sum_syn == 0
                x_dec = x_hard
                iter += 1
                break
            end
        end
        iter += 1
    end

    iter -= 1
    HardDecode(x_dec, x_hard)
    x_dec = x_hard
    return x_dec, iter
end


function HardDecode(x, x_hard)

    ###
    # This function assigns binary 0/1 values to a vector of continous values. If a value is above 0.5, it is set to 1,
    # otherwise, it is set to 0. The function updates x_hard in-place to enhance performance.
    #
    # x: input vector with continous variables
    # x_hard: the output binary vector, it should be allocated with the same size as x
    #
    # Returns:
    # nothing
    #
    ###

    for i =1:N
        if x[i] <= 0.5
            x_hard[i] = 0;
        else
            x_hard[i] = 1;
        end
    end
    return nothing
end
