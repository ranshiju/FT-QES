# FT-QES
* Codes for calculating the physical-bath Hamiltonians of 1D and 2D (homeycomb) models  
* Ref: "Quantum simulation for thermodynamics of infinite-size many-body systems by O(10) sites" ([arXiv:1810.01612](https://arxiv.org/abs/1810.01612))  
* The Hamiltonains are given by Hermitian matrices in the computational basis (the eigenstates of the z-component of spin-1/2 operator）  

## For 1D chain, see 'QES_1D_bath_H'  
Step 1. prepare 'Parameter.m'; one may edit the physical and computational parameters as your wish;  
Step 2. run 'MainBH.m';  
Step 3. Load the mat file in the current folder; the physical-bath Hamiltonians are saved inside as a cell.  

## For 2D honeycom lattice, see 'QES_honeycomb_bath_H'  
Step 1. prepare Parameter.m. One may edit the physical and computational parameters as your wish (first 22 rows);  
Step 2. run 'MainBH.m';  
Step 3. load the mat file in the current folder; the physical-bath Hamiltonians are saved inside as a cell.
