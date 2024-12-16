N=27895
M=516348
B=126
C=2
P=6
import numpy as np
from scipy.sparse import load_npz

X_temp2=[]
for j in range(B):
    matrices=[]
    for p in range(P):
        X_temp = load_npz(
            '../test_site/N{}_M{}_C{}_P{}_B{}/Party_{}/X_block_{}.npz'.format(N, M, C,
                                                                                                              P, B,
                                                                                                              p+1,
                                                                                                              j +1))
        X_temp=X_temp.toarray()
        matrices.append(X_temp)
    X_temp2.append(np.vstack(matrices))
X=np.hstack(X_temp2)



matrices=[]
for j in range(P):
    y_t=np.load('../test_site/N{}_M{}_C{}_P{}_B{}/Party_{}/y.npy'.format(N, M, C, P, B,j+1))
    print(y_t.shape)
    matrices.append(y_t)
y=np.vstack(matrices)
print(y.shape)

matrices=[]
for j in range(P):
    Z_t=np.load('../test_site/N{}_M{}_C{}_P{}_B{}/Party_{}/Z.npy'.format(N, M, C, P, B, j+1))
    matrices.append(Z_t)
Z=np.vstack(matrices)

import os
individuals_per_family = N // 6

directory_path = '../test_site/N{}_M{}_C{}_P{}_B{}/REGENIE/'.format(N, M, C, P, B)
if not os.path.exists(directory_path):
    os.makedirs(directory_path, exist_ok=True)

file_path = os.path.join(directory_path, 'phenotype.txt')

with open(file_path, "w") as f:
    f.write("FID IID Y1\n")

    y_index = 0

    for family in range(1, 7):
        num_individuals = individuals_per_family if family < 6 else len(y) - 5 * individuals_per_family
        print(y.shape)
        print(y[50:1000])
        for individual in range(1, num_individuals + 1):
            FID = f"FAM{family}ID{individual}"
            IID = f"IND{family}ID{individual}"
            Y1 = y[y_index,0]  
            y_index += 1

            f.write(f"{FID} {IID} {Y1:.14f}\n")


print(f"Output file written to {file_path}")

with open('../test_site/N{}_M{}_C{}_P{}_B{}/REGENIE/covariates.txt'.format(N, M, C, P, B), "w") as f:
    f.write("FID IID V1 V2\n")
    y_index=0
    for family in range(1, 7):
        num_individuals = individuals_per_family if family < 6 else len(y) - 5 * individuals_per_family
        for individual in range(1, num_individuals + 1):
            FID = f"FAM{family}ID{individual}"
            IID = f"IND{family}ID{individual}"
            V1= np.squeeze(Z[y_index,0])
            V2= np.squeeze(Z[y_index,1])
            y_index += 1
            f.write(f"{FID} {IID} {V1} {V2} \n")

print("Output file written to covariates.txt")

import numpy as np
from scipy.sparse import load_npz
import gc 


def convert_to_genotype(snp):
    return 'A A' if snp == 0 else 'A C' if snp == 1 else 'C C'

with open(f'../test_site/N{N}_M{M}_C{C}_P{P}_B{B}/REGENIE/output.map', 'w') as f_map:
    for i in range(M):
        f_map.write(f'1 snp{i} 0 {i * 100}\n')

with open(f'../test_site/N{N}_M{M}_C{C}_P{P}_B{B}/REGENIE/output.ped', 'w') as f_ped:
    for p in range(P):
        X_party = []
        for j in range(B):
            path = f'../test_site/N{N}_M{M}_C{C}_P{P}_B{B}/Party_{p + 1}/X_block_{j + 1}.npz'
            X_block = load_npz(path).toarray()
            X_party.append(X_block)

        X_party_combined = np.hstack(X_party)

        for i in range(X_party_combined.shape[0]):
            y_formatted = '{:.14f}'.format(y[i][0])
            f_ped.write(f'FAM{p+1}ID{i+1} IND{p+1}ID{i+1} 0 0 1 {str(y_formatted)}')
            for k in range(X_party_combined.shape[1]):
                snp_val = X_party_combined[i, k]
                genotype = convert_to_genotype(snp_val)
                f_ped.write(f' {genotype}')
            f_ped.write('\n')

            sample_progress = (i + 1) / X_party_combined.shape[0] * 100
            print(f"Party {p+1}, Sample {i+1}/{X_party_combined.shape[0]}: {sample_progress:.2f}% completed")

        print(f"Party {p+1} processed")

        # Free memory
        del X_party
        del X_party_combined
        gc.collect()  
      
print("All parties processed.")
