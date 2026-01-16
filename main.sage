#=======================================================================================#
# Supporting implementation of Code-Based Homomorphic Encryption for Secure Aggregation #
#=======================================================================================#

load('LPNAgg.sage')

# parameters

x = 2^16 # range of each user
N = 50 # number of users
total_range = N*(x-1)
k_info = 10^12 # data dimension
lam = 128 # security level

#solver = 0 # Thm 1
solver = 1 # Thm 2

step = max(floor( (N-2)/20) ,1)
zs = [i for i in range(5)] + [i for i in range(5,N-20,step)] + [i for i in range(max(N-20,5),N-1)] 

for z in zs:
    total_comm, qs, ps, ks = optimize_user_comm_crt(k_info, N, z, total_range, lam, solver)
    ps_plot = [numerical_approx(p,16) for p in ps]
    ks_log10 = [numerical_approx(math.log(k,10),16) for k in ks]
    print(f'% z = {z}: qs = {qs}, ps = {ps_plot}, ks = 10^{ks_log10}')
    print(f'{numerical_approx(z/N)} {ceil(total_comm)}\\\\')


print(f'\n% SwiftAgg+:')
for z in range(0,N-1,step):
    swiftagg_comm = swiftagg_plus_comm(k_info, total_range, N, z)
    print(f'{numerical_approx(z/N)} {ceil(swiftagg_comm)}\\\\')


