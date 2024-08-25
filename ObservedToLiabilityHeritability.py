# As disease is observed from 0 to 1, we assume there is a liability model behind the disease that
# if the additive genetic value is larger than the threshold then the individual will have disease
# and the lower one does not. Observed heritablity therefore should be converted to liability scale
# to constuct the genetic architecture of the disease.

# To do this, we will use equation 23 from Sang Hong Lee et al, 2011.

from scipy import stats

h_obs = float(input("Enter observed heritability value: "))
K = float(input("Enter disease prevalence or number of cases worldwide: "))
P = float(input("Enter case/total ratio in the GWAS study: "))

if K > 1:
    K = K/(8*(10**9))
z = stats.norm.pdf(stats.norm.ppf(K,0,1),0,1)
# The height of the threshold in the liability threshold model.

h_lia = h_obs * (K*(1-K))**2 / (P*(1-P))/ (z**2)

print("Heritablity of the liability scale is: ",  str(h_lia))