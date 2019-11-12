! Module containing functions and subroutines for working with tensors in 9-voigt notation
! Convention used in these routines are the following
! Voigt 9 element notation [11, 22, 33, 12, 23, 13, 31, 21, 32]
! Voigt 6 element symmetric notation [11, 22, 33, 12, 23, 13]
! Thus                       |1 4 7|
! [1 2 3 4 5 6 7 8 9]   =>   |8 2 5|
!                            |6 9 3|
!
!==============================================================================
!        Summary of functions and subroutines contained in the module
!==============================================================================
!
! CONVERSIONS
! V6_2_V9           Transformation of a 2nd order tensor from 6 to 9 components in Voigt format
! V9_2_V6           Transformation of a 2nd order tensor from 9 to 6 components in Voigt format
! V9_2_M            Transformation of a 2nd order tensor with 9 components to matrix form of size(3,3)
! V6_2_M            Transformation of a 2nd order tensor with 6 components to matrix form of size(3,3)
! M_2_V9            Transformation of a matrix of size(3,3) to a second order tensor form with 9 components
! M_2_V6            Transformation of a matrix of size(3,3) to a second order tensor form with 9 components
! V9x9_2_V6x6       Transformation of a 4th order tensor of size(9,9) to size(6,6)
! V6x6_2_Sym_V9x9   Transformation of a 4th order tensor of size(6,6) to size(9,9)
! V9x9_2_T4         Transformation of a 4th order tensor of size(9,9) to pure tensor form of size(3,3,3,3)
! T4_2_V9x9         Transformation of a pure tensor form of size(3,3,3,3) to 4th order tensor of size(9,9) 
! minsym_V9x9       Minor symmetric part of 9x9 4th order tensor, not implemented. Required for V9x9_2_V6x6, thus in this file
! majsym_V9x9       Major symmetric part of 9x9 4th order tensor, not implemented. Put in this file to correspond to minsym_V9x9
! ----------------------------------------------------------------------------------------------------------
!
! OPEN_PRODUCTS
! op_ve             Open product between two 1st order tensors (3 components) [i x j] giving matrix 3x3 out
! op_ve_9           Open product between two 1st order tensors (3 components) [i x j] giving 9 voigt out
! op_V9             Open product between two 2nd order tensors (9 components) [ij x kl]
! op_a_V9           Non-standard open product "above" between two 2nd order tensors (9 components) [ik x jl]
! op_b_V9           Non-standard open product "below" between two 2nd order tensors (9 components) [il x jk]
! Pop_V9            c_ijkl = 0.25*( a_ik b_jl+a_il b_jk+b_ik a_jl+b_il a_jk )
!
! Open products between 4th and - 2nd order tensors
! v9_op_v9x9        a_ij B_klmn - 2nd order in first position
! v9x9_op_v9        A_ijkl b_mn - 2nd order in last position
! v9_opc_v9x9       c_kl D_ijmn - 2nd order in center position
! ilmn_jk           c_ijklmn = a_ilmn * b_jk
! ikmn_jl           c_ijklmn = a_ikmn * b_jl
!
! 6-component functions
! op_V6             Open product between two 2nd order tensors (6 components)
! ----------------------------------------------------------------------------------------------------------
!
! CONTRACTIONS
! V9_d_v3           Contraction between 2nd order tensor and 1st order tensor
! V9_d_V9           Contraction between two 2nd order tensors with 9 components
! V9_dd_V9          Double contraction between two 2th order tensor of size (9)
! V9x9_dd_V9x9      Double contraction between two 4th order tensor of size (9,9)
! V9x9_dd_V9        Double contraction between a 4th order and a 2nd order tensor of size (9,9) and (9) respectively
! V9_dd_V9x9        Double contraction between a 2nd order and a 4th order tensor of size (9) and (9,9) respectively
!
! Special cases
! V9dV9dV9          Contraction between three 2nd order tensors with 9 components
! V9_d_V9x9         Single contraction between a 2nd and a 4th order tensor
! V9x9_d_V9         Single contraction between a 4th and a 2nd order tensor
!
! Operations on ordinary matrices
! M_x_ve            Multiplication between a matrix and vector
! ve_x_M            Multiplication between a vector and matrix
! ----------------------------------------------------------------------------------------------------------
!
! STANDARD_TENSORS (PROJECTION TENSORS)
! eye_V9            2nd order identity tensor (9 components)
! idev_v9x9         4th order deviatoric identity tensor
! ----------------------------------------------------------------------------------------------------------
!
! OPERATIONS        
! trans_V9          Transpose of a 2nd order tensor with 9 components           
! dev_V9            Deviatoric part of a 2nd order tensor with 9 components     
! tr_V9             Trace of a 9 component voigt vector                         
! det_V9            Determinant of 2nd order 9-component tensor                 
! inv_V9            Inverse of 2nd order 9-component tensor                     
!
! Operations on full matrices
! det_M             Determinant of 1x1, 2x2 or 3x3 matrix                       
! inv_M3            Subroutine to calculate inverse of 3x3 matrix               
! inv_M             Subroutine to calculate inverse of general matrix           
! trans_M           Transpose of square matrix                                  
! tr_M              Trace of a general square matrix                            
! norm_V3           Norm of a 1st order tensor with 3 components
! norm_V9           Norm of a general vector
! ----------------------------------------------------------------------------------------------------------
!
! STIFFNESS_CONVERSIONS
! dSdE_2_dPdF       Convert dSdE to dPdF
! dSdE_2_dtaudF     Convert dSdE to d(tau)dF
! dPdF_2_dtaudF     Convert dPdF to d(tau)dF
! ----------------------------------------------------------------------------------------------------------
!
! ABAQUS CONVERSIONS
!   Functions               Input           Description
! sigma_v9_2_sigma_abaqus   sigma_v9        Convert sigma in 9-voigt to abaqus sigma 6-voigt and abaqus index order
! sigma_abaqus_2_sigma_v9   sigma_abaqus    Convert sigma in abaqus sigma 6-voigt and abaqus index order to 9-voigt
! dtaudF_2_DDSDDE           dtaudF, F       Convert dtaudF in 9x9voigt to Abaqus DDSDDE in 6x6 voigt and abaqus index order
! ----------------------------------------------------------------------------------------------------------
! 
! MISCELLANEOUS
! jacobi_mekh       ADD DESCRIPTION!
! eigsrt            ADD DESCRIPTION!
! M_pow             Calculate C^m where C is 3x3 matrix
! eig_prob          ADD DESCRIPTION!
! us_solve          ADD DESCRIPTION!
! us_gausjor        ADD DESCRIPTION!
!==============================================================================

include 'Tensors_module/Conversions.f90'       !Must always be included, required for multiple other modules
include 'Tensors_module/Open_products.f90'
include 'Tensors_module/Contractions.f90'
include 'Tensors_module/Standard_tensors.f90'
include 'Tensors_module/Operations.f90'
include 'Tensors_module/Stiffness_conversions.f90'
include 'Tensors_module/Abaqus_conversions.f90'
include 'Tensors_module/Miscellaneous.f90'


Module Tensors_module
use Conversions
use Open_products
use Contractions
use Standard_tensors
use Operations
use Stiffness_Conversions
use ABAQUS_Conversions
use Miscellaneous_tensors

Implicit none
    
END MODULE