# Shi2014exp_BC
This model can have 2 back-stresses and is rate-independent. It only supports Armstrong-Frederick type of kinematic hardening. For a detailed description of the model formulation, please see the paper in ref.bib

# parameters
1. G (Shear modulus)
2. K (Bulk modulus)
3. Y0 (Initial yield limit) 
4. Hiso (Isotropic hardening modulus)
5. Rinf (Isotropic saturation stress)
6. delta (Fraction of Armstrong-Frederick (delta=1) and Burlet-Cailletaud (delta=0) type of kinematic hardening)
7. Hkin1 (Kinematic hardening modulus nr. 1)
8. Binf1 (Kinematic saturation stress nr 1)
9. Hkin2 (Kinematic hardening modulus nr. 2)
10. Binf2 (Kinematic saturation stress nr 2)
11. bD (Amount of dynamic hardening)
12. bL (Amount of latent hardening)
13. cD (Evolution parameter for CD)
14. cL (Evolution parameter for CL)


# state variables
1-9) Fp-I2 (Plastic deformation gradient (minus 2nd order identity))
10) lambda (Accumulated plastic deformation (time integral of plastic multiplier))
11-19) Fk1-I2 (Kinematic deformation gradient nr1 (minus 2nd order identity))
20-28) Fk2-I2 (Kinematic deformation gradient nr2 (minus 2nd order identity))
29-49) CD - I4 (Dynamic hardening 4th order tensor (minus 4th order symmetric identity tensor). Saved in special format using that it is minor and major symmetric)
50-70) CL - I4 (Latent hardening 4th order tensor (minus 4th order symmetric identity tensor). Saved in special format using that it is minor and major symmetric)