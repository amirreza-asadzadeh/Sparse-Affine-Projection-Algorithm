using DelimitedFiles

F = readdlm("LDPC Codes/WiGig/WifiCode.txt")

const M = size(F)[1]
F_list_temp = Array{Array{Int, 1}, 1}(undef, M)
CheckDegree_temp = zeros(Int, M)

for j = 1:M
    vec = F[j,:]
    L = count(i -> (i != ""), vec)
    CheckDegree_temp[j] = L
    vec1 = zeros(L)
    for i = 1:L
        temp = vec[i]
        vec1[i] = temp + 1
    end
    F_list_temp[j] = vec1
end

const CheckDegree = CheckDegree_temp
const F_list = F_list_temp

# computing number of variable nodes
max_vec = zeros(Int, M)
for j = 1:M
    max_vec[j] = findmax(F_list[j])[1]
end
const N = findmax(max_vec)[1]

# computing degree of each variable node
# VariableDegree = Array{Int, 1}(undef, N)
VariableDegree_temp = zeros(Int, N)
H_temp = zeros(Int, M, N)
G_list_temp = Array{Array{Int, 1}, 1}(undef, N)
for i = 1:N
    G_list_temp[i] = []
end
for j = 1:M
    for i = 1:CheckDegree[j]
        VN_temp = F_list[j][i]
        H_temp[j, VN_temp] = 1
        VariableDegree_temp[VN_temp] += 1
        G_list_temp[VN_temp] = push!(G_list_temp[VN_temp], j)
    end
end
const H = H_temp
const VariableDegree = VariableDegree_temp
const G_list = G_list_temp

println("N = $N")
println("M = $M")
