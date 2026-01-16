#=======================================================================================#
# Supporting implementation of Code-Based Homomorphic Encryption for Secure Aggregation #
#=======================================================================================#

from copy import deepcopy
import math

# high-precision arithmetic
RR = RealField(1000) 

def get_piling(N, q, p_lpn):
    """
    Piling-up lemma: computes error probability of sum of Bernoullis
    
    :param N: number of additions
    :param q: field size
    :param p_lpn: individual error probability
    """
    p = RR(p_lpn)
    return (q-1)/q - (q-1)/q * (1 - q/(q-1)*p)^N

def get_p_hint(N, z, q, p_lpn, solver):
    """
    Compute error probability given hint = 0

    :param N: number of parties
    :param z: collusion parameter
    :param q: field size
    :param p_lpn: individual error probability
    """
    p_lpn = RR(p_lpn)

    assert N - z >= 2, "invalid number of collusions"

    if solver:
        # Thm 2: concrete attack
        p1 = p_lpn
        p2 = get_piling(N-z-1, q, p_lpn)
    
        return p1*p2/((q-1)^2*(1-p1)*(1-p2) + p1*p2)
    
    else: # Thm 1: reduction
        p = p_lpn    
        return p^2/((q-1)^2*(1-p)^2 + p^2)


def capacity(q, p):
    """
    Capacity of q-ary symmetric channel,
    used to compute rate of error-correcting code.
    
    :param q: field size
    :param p: error probability
    """
    return 1 - q_entropy(p, q)


def q_entropy(p, q):
    """
    q-ary entropy fucntion
    
    :param p: error probability
    :param q: field size
    """
    if p == 0. or p == 1.:
        return 0.
    
    if p < 0. or p > 1.:
        return -1000
    
    return RR(p*log(q-1, q) - p*log(p, q) - (1-p)*log(1-p, q))
    


def prange_k(p, lam):
    """
    Somputes LPN dimension that achieves lambda bits of security against Prange/Gauss attack.
    See Esser, Kübler, and May, "LPN Decoded" for details.

    :param p: LPN noise rate
    :param lam: security parameter
    """
    if p <= 0:
        return float('inf')
    
    return -lam / log(1-RR(p),2)

def user_comm(k_info, p_lpn, N, z, q, lam, solver):
    """
    compute per-user communication cost for aggregation with
    * LPN-based KAHE
    *  committee-based decryptor

    :param k_info: information dimension
    :param p_lpn: individual error probability
    :param N: number of parties
    :param z: collusion parameter
    :param q: field size
    :param lam: security parameter
    """
    assert N >= z+2, "invalid collusion parameter"
    p_lpn = RR(p_lpn)

    # security of KAHE
    p_hint = get_p_hint(N, z, q, p_lpn, solver)
    k_lpn = prange_k(p_hint, lam)

    # correctness of KAHE
    p_total = get_piling(N, q, p_lpn)
    R = capacity(q, p_total)

    if R <= 0:
        return float('inf')
    else:
        comm_KAHE = k_info / R * math.log(q,2)

    # communication to decryptor assuming swiftagg+
    comm_decrpytor = swiftagg_plus_comm(k_lpn, q, min(q, N), z)

    return comm_decrpytor + comm_KAHE

def optimize_user_comm(k_info, N, z, q, lam,solver):
    """
    optimize per-user communication cost over the noise rate of the KAHE

    :param k_info: information dimension
    :param N: number of parties
    :param z: collusion parameter
    :param q: field size
    :param lam: security parameter
    """

    from scipy.optimize import minimize_scalar
    p_min = 1e-10
    p_max = 0.3
    res = minimize_scalar(lambda p_lpn: user_comm(k_info, p_lpn, N, z, q, lam, solver), bounds=(p_min, p_max), method='bounded')

    p_lpn = res.x

    p_hint = get_p_hint(N, z, q, p_lpn, solver)
    k_lpn = prange_k(p_hint, lam)

    return res.fun, p_lpn, k_lpn

def optimize_user_comm_crt(k_info, N, z, total_range, lam, solver, debug = 0):
    """
    optimize per-user communication cost using the Chinese Remainder Theorem (CRT)-technique
    
    :param k_info: information dimension
    :param N: number of parties
    :param z: collusion parameter
    :param total_range: required modulus to avoid overflow
    :param lam: security parameter
    """

    qs = list(primes(z+1,next_prime(200)+1))

    comms = []
    ps = []
    ks = []

    L = []

    for q_ in qs:
        if debug:
            print('-------------------------')
            print(f'optimize over field q={q_}')
            print('-------------------------')
        comm_q_, p_lpn, k_lpn = optimize_user_comm(k_info, N, z, q_, lam,solver)

        
        comms.append(comm_q_)
        ps.append(p_lpn)
        ps.append(p_lpn)
        ks.append(k_lpn)

        L.append( (comm_q_, math.log(q_,2), [q_]) )

    L.sort(key=lambda x: x[1]/x[0], reverse=True)

    if debug:
        for x in L:
            print(x[0], x[1], x[1]/x[0], x[2])


        print('optimizing composition')

    dp = {0: (0, [])}

    factor = 0.002 #lower: slower but more precise optimization

    for cost, log2q, q_set in L:

        new_dp = dp.copy()

        for q_scaled, (cost_sum, q_used) in dp.items():

            new_q_scaled = q_scaled + int(round(log2q/factor))
            new_cost_sum = cost_sum + cost
            new_q_set = q_set + q_used

            # keep best per bucket only cost only
            if new_q_scaled not in new_dp or new_cost_sum < new_dp[new_q_scaled][0]:
                new_dp[new_q_scaled] = (new_cost_sum, new_q_set)

        dp = new_dp

    # select allowed combinations only
    feasible = [(comm, q_set) for _, (comm, q_set) in dp.items() if prod(q_set) >= total_range]
    
    opt = min(feasible)
    opt_q = sorted(opt[1])
    I =  [qs.index(q_) for q_ in opt_q]
    opt_p = [ps[i] for i in I]
    opt_k = [ks[i] for i in I]

    total_comm = ceil(sum(comms[i] for i in I))

    return total_comm, opt_q, opt_p, opt_k


def swiftagg_plus_comm(k_info, q, N, z):
    """
    Per-user communication cost of pure swiftagg+

    :param k_info: information dimension
    :param N: number of parties
    :param z: collusion parameter
    """
    assert z < N <= q, "Invalid parameters for SwiftAgg+"
    return (1 + z/(N-z))*k_info* math.log(q,2)