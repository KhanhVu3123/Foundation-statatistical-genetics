# Given the population number and number of recessive alleles, we aim to calculate the probability that 
# the recessive alleles will be wiped out completely from 1 to m generations.

from scipy.stats import binom
import math

# Constructing the Wright Fisher model with n diploid individuals and k recessive alleles.
# Return the probabilities of 0 to 2*n recessive alleles after 1 generation.
def WrightFisherModel_1gen(n, k):
    lst = [0] * (2*n+1)
    if(k ==0):
        lst[0] = 1
        return lst
    
    for i in range(len(lst)):
        lst[i] = binom.pmf(i, 2*n, k/(2*n))

    return lst

# Return all generations prob.
def WrightFisherModel(n,k,g):
    return_this = []
    lst = WrightFisherModel_1gen(n,k)
    return_this.append(lst)
    for j in range(g-1):
        one_gen_list = [0] * (2*n+1)
        for i in range(len(lst)):
            new_lst = WrightFisherModel_1gen(n,i)
            new_lst = [x * lst[i] for x in new_lst]
            for index in range(len(one_gen_list)):
                one_gen_list[index] += new_lst[index]
        lst = one_gen_list
        return_this.append(lst)
    return return_this

N, m = input("Enter population number and number of populations: ").split()

N, m = int(N), int(m)

no_recessive_lst = list(input("Enter the number of recessive alleles: ").split())

for i in range(len(no_recessive_lst)):
    no_recessive_lst[i] = int(no_recessive_lst[i])

ans_list = [[0] * len(no_recessive_lst) for _ in range(m)] 
# Loop through each number of recessive alleles for the population

for no_recessive in no_recessive_lst:
    no_recessive_prob = WrightFisherModel(N,no_recessive,m)
    for index in range(m):
        ans_list[index][no_recessive_lst.index(no_recessive)] = math.log10(no_recessive_prob[index][0])

for i in range(m):
    for j in range(len(no_recessive_lst)):
        print(ans_list[i][j], end=" ")
    print()

print(ans_list)
