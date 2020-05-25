
module crystal_orientations_mod
use Tensors_module
use orientation_mod
    implicit none
    private
    public  :: rotate_input
    public  :: rotate_output
    public  :: get_rotation_matrices
    
    integer, parameter, dimension(9)    :: aba6_2_v9_pos = (/1, 2, 3, 4, 6, 5, 5, 4, 6/) ! Index conversion
    integer, parameter, dimension(6)    :: v9_2_aba6_pos = (/1, 2, 3, 4, 6, 5/)          ! Index conversion

    contains
    
subroutine rotate_input(dfgrad, dfgrad_rot, qvec)
implicit none
    double precision    :: dfgrad(3,3), dfgrad_rot(3,3), qvec(9), qmat(3,3)
    
    qmat = v9_2_m(qvec)
    
    dfgrad_rot = matmul(transpose(qmat), matmul(dfgrad, qmat))
    
end subroutine

subroutine rotate_output(stress, ddsdde, qvec)
implicit none
    double precision    :: stress(6), ddsdde(6,6), qvec(9)
    double precision    :: stressmat(3,3), ddsdde9(9,9), qmat(3,3), rot9x9(9,9)
    
    ! Rotate stress
    qmat = v9_2_m(qvec)
    stressmat = v9_2_m(sigma_abaqus_2_sigma_v9(stress))             ! Convert stress abaqus voigt to matrix
    stressmat = matmul(qmat, matmul(stressmat, transpose(qmat)))    ! Rotate stress matrix
    stress = sigma_v9_2_sigma_abaqus(m_2_v9(stressmat))             ! Convert stress matrix to abaqus voigt
    
    
    rot9x9 = op_a_V9(qvec, qvec)                                    ! Create rotation matrix in voigt for 4th order tensor
    ddsdde9 = ddsdde(aba6_2_v9_pos, aba6_2_v9_pos)                  ! Convert ddsdde in abaqus indicies to 9-voigt indicies from tensor module
    ddsdde9 = matmul(rot9x9, matmul(ddsdde9, transpose(rot9x9)))    ! Rotate ddsdde9
    ddsdde = ddsdde9(v9_2_aba6_pos,v9_2_aba6_pos)                   ! Convert ddsdde9 in 9-voigt indicies to abaqus indices (6-voigt)
    
end subroutine
    
end module    