# gfs_bc
This model can have 2-4 back-stresses and is rate-independent. For a detailed description of the model formulation, please see the paper in ref.bib

# parameters
1. G (Shear modulus)
2. K (Bulk modulus)
3. Y0 (Initial yield limit) 
4. Hiso (Isotropic hardening modulus)
5. invYiso (Inverse of isotropic saturation stress)
6. delta (Scaling of armstrong frederick (delta=1) and burlet-cailletaud (delta=0))
7. Hk1 (Kinematic hardening modulus nr. 1)
8. invYk1 (Inverse of kinematic saturation stress nr 1)
9. mk1 (Not used in this model)
10. Hk2 (Kinematic hardening modulus nr. 2)
11. invYk2 (Inverse of kinematic saturation stress nr 2)
12. mk2 (Not used in this model)

...

# state variables
1) lambda (Total plastic multiplier)
2-10) Fp-I2 (Plastic deformation gradient (minus 2nd order identity))
11-19) Fk1-I2 (Kinematic deformation gradient nr1 (minus 2nd order identity))
20-28) Fk2-I2 (Kinematic deformation gradient nr2 (minus 2nd order identity))

...