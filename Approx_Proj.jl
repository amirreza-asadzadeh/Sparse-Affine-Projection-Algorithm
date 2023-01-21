using LinearAlgebra, Statistics

function probsimplexproj(w, v)

    ###
    # A function for projection onto probability simplex in the centered framework.
    # Function projects the point and assigns the output value to input for performance
    # enhancement.
    #
    # v: input point to be projected
    # w: output vector, the same size as v
    #
    # Returns:
    # nothing
    ###
    
    rho = sort(v, rev = true);
    d = length(v);
    u = Vector{Float64}(undef, d);
    sum = 0;
    for i = 1:d
        sum += rho[i];
        u[i] = (1/i) * (sum - 1);
    end
    i_star = 1;
    while( i_star <= d )
        if rho[i_star] > u[i_star]
            i_star += 1;
        else
            break
        end
    end
    i_star -= 1;
    for i = 1:d
        w[i] = max( v[i] - u[i_star] - 1/2, -1/2 );
    end
    return nothing
end


function simplexhyperproj(w, v)

    ###
    # A simplified function for projection onto probability simplex in the centered framework.
    # Function projects the point onto the supporting hyperplane of the probability simplex
    # and assigns the output value to input for performance enhancement.
    #
    # v: input point to be projected
    # w: output vector, the same size as v
    #
    # Returns:
    # nothing
    ###

    d = length(v);
#     for i = 1:d
#         v[i] = max( v[i], -1/2 );
#     end
    v_shift = mean(v) - 1/d;
    for i = 1:d
        w[i] = v[i] - v_shift - 1/2;
    end
    return nothing
end

function EXACT(v)

    ###
    # This function projects a point onto the parity polytope EXACTLY.
    # The input vector should be given in the non-centered framework. Output is given in the non-centered framework
    # as well. However, the computations inside the function takes place in the centered
    # framework, as it is practically more efficient. The function updates the input vector
    # as the output vector to enhance performance of the code.
    #
    # v: input point to be projected
    #
    # Returns:
    # the projected point
    ###

    d = length(v);
    v -= ones(d) * (1/2); #shifting to fit in the centered projection framework
    f = zeros(d);
    for i = 1:d
        if v[i] >= 0
            f[i] = 1; #facet identifictation
        end
    end
    if sum(f)%2 == 0
        i_star = argmin( abs.(v) );
        f[i_star] = 1 - f[i_star];
    end
    v_tilde = similar(v);
    for i = 1:d
        v_tilde[i] = v[i] * (-1)^f[i]; #similarity transform
    end
    u_tilde = similar(v);

    probsimplexproj(u_tilde, v_tilde); #simplex projection

    v_tildeclipped = similar(v);
    for i = 1:d
        v_tildeclipped[i] = min( max( v_tilde[i], -1/2 ), 1/2 );
    end
    w = similar(v);
    if sum(v_tildeclipped) >= 1 - d/2 #membership test
      for i = 1:d
            w[i] = min( max( v[i], -1/2 ), 1/2 );
        end
    else
        for i = 1:d
            w[i] = u_tilde[i] * (-1)^f[i]; #similarity transform
        end
    end
    #return nothing
    return w + ones(d) * (1/2) #invert shifting to fit in the uncentered projection framework
end

function APA(v)

    ###
    # This function is an approximation of the parity polytope projection, known as Affine Projeciton Algorithm (APA).
    # It approximates the projection by projecting onto the supporing hyperplane of all d
    # even-parity vertices which neighbor the closest odd-parity hypercube vertex to the point to be projected.
    # The input vector should be given in the non-centered framework. Output is given in the non-centered framework
    # as well. However, the computations inside the function takes place in the centered
    # framework, as it is practically more efficient. The function updates the input vector
    # as the output vector to enhance performance of the code.
    #
    # v: input point to be projected
    #
    # Returns:
    # the projected point
    ###

    d = length(v);
    v -= ones(d) * (1/2); #shifting to fit in the centered projection framework
    f = zeros(d);
    for i = 1:d
        if v[i] >= 0
            f[i] = 1; #facet identifictation
        end
    end
    if sum(f)%2 == 0
        i_star = argmin( abs.(v) );
        f[i_star] = 1 - f[i_star];
    end
    v_tilde = similar(v);
    for i = 1:d
        v_tilde[i] = v[i] * (-1)^f[i]; #similarity transform
    end
    u_tilde = similar(v);

    simplexhyperproj(u_tilde, v_tilde); #simple simplex hyperplane projeciton (without sorting needed)

    v_tildeclipped = similar(v);
    for i = 1:d
        v_tildeclipped[i] = min( max( v_tilde[i], -1/2 ), 1/2 );
    end
    w = similar(v);
    if sum(v_tildeclipped) >= 1 - d/2 #membership test
      for i = 1:d
            w[i] = min( max( v[i], -1/2 ), 1/2 );
        end
    else
        for i = 1:d
            w[i] = u_tilde[i] * (-1)^f[i]; #similarity transform
        end
    end
    #return nothing
    return w + ones(d) * (1/2) #invert shifting to fit in the uncentered projection framework
end


function EVA(v)

    ###
    # This function is an approximation of the parity polytope projection, namely as Even-Vertex-Algorithm (EVA).
    # It approximates the projection by projecting onto the closest even-parity vertex of the polytope to the point to be projected.
    # The input vector should be given in the non-centered framework. Output is given in the non-centered framework
    # as well. However, the computations inside the function takes place in the centered
    # framework, as it is practically more efficient.
    #
    # v: input point to be projected
    #
    # Returns:
    # f: the projected point
    ###

    d = length(v);
    v -= ones(d) * (1/2);
    f = zeros(d);
    for i = 1:d
        if v[i] >= 0;
            f[i] = 1;
        end
    end
    if sum(f)%2 == 1
        i_star = argmin( abs.(v) );
        f[i_star] = 1 - f[i_star];
    end
    return f
end


function LSA(v)

    ###
    # This function is an approximation of the parity polytope projection, namely as Line-Segment-Algorithm (LSA).
    # It approximates the projection by projecting onto the line-segment connecting the two closest even-parity vertices neighboring
    # the nearest odd-parity hypercube vertex to the point to be projected.
    # The input vector should be given in the non-centered framework. Output is given in the non-centered framework
    # as well. However, the computations inside the function takes place in the centered
    # framework, as it is practically more efficient.
    #
    # v: input point to be projected
    #
    # Returns:
    # z: the projected point
    ###

    d = length(v);
    u = similar(v);
    theta = -ones(d);
    f = zeros(d);
    for i = 1:d
        u[i] = min( max(v[i], 0.0), 1.0);
        if v[i] >= 0.5
            f[i] = 1
            theta[i] = 1;
        end
    end
    if sum(f)%2 == 0
        i_star = argmin( abs.(v.-0.5) );
        theta[i_star] *= -1;
        f[i_star] = 1 - f[i_star];
    end
    s = sum(f);
    if sum(theta.*u) <= s-1
        z = u;
    else #Line-Segment Projection
        sort_idx = sortperm(abs.(v.-0.5));
        p = sort_idx[1];
        q = sort_idx[2];
        A = copy(f);
        B = copy(f);
        A[p] = 1 - f[p];
        B[q] = 1 - f[q];
        AB = zeros(d);
        AB[p] = B[p] - A[p];
        AB[q] = B[q] - A[q];
        t = 0.5 * ( (B[p]-A[p])*(v[p]-A[p]) + (B[q]-A[q])*(v[q]-A[q]) );
        t = min( max(t, 0.0), 1.0);
        z = A .+ t*AB;
    end
    return z
end

# Limited-Hyperplane Projection Algorithm (χ-HYPPA)
function χ_SAPA(v, χ)

    ###
    # This function is an approximation of the parity polytope projection, namely as χ-Sparse Affine Projection Algorithm (χ-SAPA).
    # It approximates the projection by projecting onto the supporting hyperplance of the χ closest even-parity vertices neighboring
    # the nearest odd-parity hypercube vertex to the point to be projected.
    # The input vector should be given in the non-centered framework. Output is given in the non-centered framework
    # as well. However, the computations inside the function takes place in the centered
    # framework, as it is practically more efficient. The function updates the input vector
    # as the output vector to enhance performance of the code.
    #
    # v: input point to be projected
    #
    # Returns:
    # w_non: the projected point
    ###

    d = length(v)
    v -= ones(d) * (1/2) #shifting to fit in the centered projection framework
    f = zeros(d);
    for i = 1:d
        if v[i] >= 0
            f[i] = 1; #facet identifictation
        end
    end

    if sum(f)%2 == 0
        i_star = argmin( abs.(v) );
        f[i_star] = 1 - f[i_star];
    end

    I_σ = sortperm( abs.( v ) )
    I_χ = I_σ[1:χ] # choosing the indices of χ closest even-parity vertices

    ṽ = similar(v)
    for i = 1:d
        ṽ[i] = v[i] * (-1)^f[i]; #similarity transform
    end
    ũ = similar(v)

    v_par_sum = 0
    for i in I_χ
        v_par_sum += ṽ[i]
    end
    v_shift = ( v_par_sum - 1 ) / χ

    for i = 1:d
        if i ∈ I_χ
            ũ[i] = ṽ[i] - v_shift - 1/2
        else
            ũ[i] = - 1/2
        end
    end

    ṽ_clipped = similar(v);
    for i = 1:d
        ṽ_clipped[i] = min( max( ṽ[i], -1/2 ), 1/2 );
    end

    w = similar(v);
    if sum(ṽ_clipped) >= 1 - d/2 #membership test
      for i = 1:d
            w[i] = min( max( v[i], -1/2 ), 1/2 );
        end
    else
        for i = 1:d
            w[i] = ũ[i] * (-1)^f[i]; #similarity transform
        end
    end
    w_non = w + ones(d).*0.5
    return w_non
end
