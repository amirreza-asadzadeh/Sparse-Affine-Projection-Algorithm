function UpdateLLR(input, gamma, ChannelName, ChannelParameter)

###
# This function generates the log-likelihood ratio (LLR) vector based on the input vector, which is
# actually the noisy recieved vector. Depending on the channel, Binary Symmetric Channel (BSC)
# Additive White Gaussian Noise Channel (AWGN), the llr vector is generated and will be updated to
# the gamma vector. So it is important to allocate a vector of proper size for gamma before calling
# this function. This function updates the vector instead of returning it to enhance the performance
# the code.
#
# input: Noisy received vector
# gamma: an empty array with the same size as input vector
# ChannelName: name of the channel, either BSC or AWGN is allowed
# ChannelParameter: SNR(dB) for AWGN channel, croos-over probability for BSC
#
# Returns:
# nothing
###

    if occursin("BSC", ChannelName)
        p = ChannelParameter;
        gamma = fill!(gamma, log( p / (1 - p) ));
        RandomVector = rand(N);
        for i = 1:N
            if ( RandomVector[i] <= p && input[i] == 1 ) || ( RandomVector[i] > p && input[i] == 0)
                gamma[i] = - gamma[i];
            end
        end
    elseif occursin("AWGN", ChannelName)
        var = 10^( - ChannelParameter / 10 ) ; #ChannelParameter is SNR in dB
        RandomVector = randn(N) * sqrt(var);
        for i= 1:N
            gamma[i] = ( ( 2*input[i] - 2 + RandomVector[i] )^2 - ( 2*input[i] + RandomVector[i] )^2 )/( 2 * var );
        end
    else
        println("Error found in Channel Name.")
    end
    return nothing
end
