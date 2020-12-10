import numpy as np

# Thermal energy in kBT 2.494353 kJ/mol at 300 K
kBT=2.494353
beta=1.0/kBT

# Total number of fes files in folder
total_files=201
# Min and max initial guesses
state_A_min = 0.200
state_A_max = 0.359
state_B_min = 0.361
state_B_max = 0.800

for i in range(total_files):
        file_name="fes_" + str(i) + ".dat"
        matrix=np.genfromtxt(file_name)
        cv=matrix[:,0]
        fes=matrix[:,1]
        fes -= np.amin(fes)
        sum_A = 0.0
        sum_B = 0.0
        for k in range(len(cv)):
            if ( cv[k] >= state_A_min and cv[k] <= state_A_max ):
                sum_A += np.exp(-beta*fes[k])
            elif ( cv[k] >= state_B_min and cv[k] <= state_B_max ):
                sum_B += np.exp(-beta*fes[k])
        FE_Diff_AB = 0.0
        if (sum_B > 0.0 and sum_A > 0.0): FE_Diff_AB = -kBT*np.log(sum_A/sum_B)
        print(str(i) + " " + str(FE_Diff_AB))
