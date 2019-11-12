# ohnowang
This model can have up to 4 back-stresses and is rate-independent

# parameters
1. E (Elastic modulus)
2. nu (Poissons ratio)
3. Y0 (Initial yield limit) 
4. Hiso (Isotropic hardening modulus)
5. invYiso (Inverse of isotropic saturation stress)
6. delta (Scaling of armstrong frederick (delta=1) and burlet-cailletaud (delta=0))
7. Hk1 (Kinematic hardening modulus nr. 1)
8. invYk1 (Inverse of kinematic saturation stress nr 1)
9. mk1 (Exponent in kinematic evolution law nr 1)
10. Hk2 (Kinematic hardening modulus nr. 2)
11. invYk2 (Inverse of kinematic saturation stress nr 2)
12. mk2 (Exponent in kinematic evolution law nr 2)
...

# state variables
1) kappa (Isotropic hardening stress)

2) lambda (plastic multiplier)

3-8) bs1 (back-stress 1)

9-14) bs2 (back-stress 2)

...