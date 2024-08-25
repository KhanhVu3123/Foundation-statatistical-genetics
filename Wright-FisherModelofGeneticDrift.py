from scipy.stats import binom

# N is the number of individuals, m is the number of dominant alleles
# g is the number of generation in our problem and k is the number of recessive alleles that the population after g generations will have at least.

# We assume that an individual can mate with itself, the species is diploid.

N,m,g,k = (input().split())

N,m,g,k = int(N), int(m), int(g), int(k)

lst = [0]*(2*N+1)
# lst[i] is Probability that the population has i recessive alleles.

lst[2*N-m] = 1

# Loop through all generations
for i in range(g):
    new_lst = [0]*(2*N+1)
    # Loop through each prob of no of recessive alleles.
    for j in range(len(lst)):
        # The list probability of getting recessive alleles from j recessive alleles.
        prob_exactly = [float(binom.pmf(noRecessive,2*N,j/(2*N)) )* lst[j] for noRecessive in range(2*N+1)]
        for index in range(len(prob_exactly)):
            new_lst[index] = new_lst[index] + prob_exactly[index]
    lst = new_lst

ans = 0
for i in range(k,len(lst)):
    ans = ans + lst[i]

print(round(ans,3))    

