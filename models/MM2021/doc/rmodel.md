# rmodel_ISO2_BC2

This model is generated in ``acegen_rmodel``. It can have 2 back-stresses and is rate-independent. The kinematic hardening is of Delobelle type (Combined AF and BC). For a detailed description of the model formulation, please see the paper in ref.bib, where this model is denoted by `H_r\neq 0` 

# parameters

1. G (Shear modulus)
2. K (Bulk modulus)
3. Y0 (Initial yield limit) 
4. k1 (Isotropic hardening constant 1)
5. Rinf1 (Isotropic saturation stress 1)
6. k2 (Isotropic hardening constant 2)
7. Rinf2 (Isotropic saturation stress 2)
8. delta (Amount of AF vs BC)
9. Hkin1 (Kinematic hardening modulus nr. 1)
10. Binf1 (Kinematic saturation stress nr 1)
11. Hkin2 (Kinematic hardening modulus nr. 2)
12. Binf2 (Kinematic saturation stress nr 2)
13. cc (Evolution parameter for C_c)
14. bc (Amount of cross hardening)
15. Hr (hardening modulus for r)
16. rinf (saturation value of r)


# state variables

1-9) Fp-I2 (Plastic deformation gradient (minus 2nd order identity))

10) lambda (Accumulated plastic deformation (time integral of plastic multiplier))

11-19) Fk1-I2 (Kinematic deformation gradient nr1 (minus 2nd order identity))

20-28) Fk2-I2 (Kinematic deformation gradient nr2 (minus 2nd order identity))

29-49) Cc (4th order cross hardening tensor. It is major and minor symmetric, and thus only 21 components [Note, also deviatoric and could theoretically be saved with only 15 components])

50-58) Fr-I2 (r deformation gradient minus 2nd order identity)