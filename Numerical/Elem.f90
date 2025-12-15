! -------------------------------------------------------------------
! File: Elem.f90
! COPYRIGHT (C) Erik Lund, Department of Materials and Production,
! Aalborg University, Denmark.
!
! List of modules, subroutines and functions:
!   module Elem
!     subroutine ElemStiffnessMatrix
!     subroutine ElemStresses
!     subroutine GetElemData
!     subroutine GetElemDisp
!     subroutine ElemLoadVector
!     subroutine StrainDisplacementJacobian
!     subroutine ConstitutiveMatrix
!
! Document:   "DMS7FEA: A Static FE Analysis Program Using Isoparametric 2D Elements" by Erik Lund
! Revision    date    engineer   review
!   0A:       4/10/99 EL         Initial version
!   0B:      22/10/03 EL         Revised version
!   0C:       3/12/04 EL         A few comments added in ElemLoadVector
!   0D:       4/11/25 DMS1-1     Added last two subroutines and finished ElemStiffnessMatrix and ElemStresses
! -------------------------------------------------------------------


module Elem
  use DataTypes
  implicit none

  public  :: ElemStiffnessMatrix, ElemStresses, &
             GetElemData, GetElemDisp, ElemLoadVector
             
!  private :: 

  contains

! -------------------------------------------------------------------
! Subroutine to calculate StrainDisplacementJacobian(Xi, Eta)
!
! This subroutine calculates the strain-displacement matrix and the Jacobian at the Gauss points
!
! Author: DMS1-1 2025
! -------------------------------------------------------------------

subroutine StrainDisplacementJacobian(E, Xi, Eta, B, detJ)

  implicit none  

  ! Defining global variables
  type(ElemType), intent(in)                          :: E
  real(double), intent(in)                            :: Xi, Eta
  real(double), dimension(3,2*E%nNodes), intent(out)  :: B
  real(double), intent(out)                           :: detJ
  
  ! Defining local variables
  real(double), dimension(2,E%nNodes)                 :: dN
  real(double), dimension(2,2)                        :: J, invJ
  real(double), dimension(4,2*E%nNodes)               :: B3
  real(double), dimension(4,4)                        :: B2
  real(double), dimension(3,4)                        :: B1
  
  ! Fill in the derivative matrix dN (derivatives of shape functions with respect to Xi and Eta)
  select case (E%nNodes) ! Defines dN based on number of nodes in element
    case (4)
      ! Derivatives w.r.t Xi
      dN(1,1) = -0.25d0 * (1.0d0 - Eta)
      dN(1,2) =  0.25d0 * (1.0d0 - Eta)
      dN(1,3) =  0.25d0 * (1.0d0 + Eta)
      dN(1,4) = -0.25d0 * (1.0d0 + Eta)
      ! Derivatives w.r.t Eta
      dN(2,1) = -0.25d0 * (1.0d0 - Xi)
      dN(2,2) = -0.25d0 * (1.0d0 + Xi)
      dN(2,3) =  0.25d0 * (1.0d0 + Xi)
      dN(2,4) =  0.25d0 * (1.0d0 - Xi)
    
    case(8)  
      ! Derivatives w.r.t Xi
      dN(1,8) = -0.5d0 * (1.0d0 - Eta**2)
      dN(1,7) = -Xi * (1.0d0 + Eta)
      dN(1,6) =  0.5d0 * (1.0d0 - Eta**2)
      dN(1,5) = -Xi * (1.0d0 - Eta)
      dN(1,4) = -0.25d0 * (1.0d0 + Eta) - 0.5d0 * dN(1,7) - 0.5d0 * dN(1,8)
      dN(1,3) =  0.25d0 * (1.0d0 + Eta) - 0.5d0 * dN(1,6) - 0.5d0 * dN(1,7)
      dN(1,2) =  0.25d0 * (1.0d0 - Eta) - 0.5d0 * dN(1,5) - 0.5d0 * dN(1,6)
      dN(1,1) = -0.25d0 * (1.0d0 - Eta) - 0.5d0 * dN(1,5) - 0.5d0 * dN(1,8)
      ! Derivatives w.r.t Eta
      dN(2,8) = -Eta * (1.0d0 - Xi)
      dN(2,7) =  0.5d0 * (1.0d0 - Xi**2)
      dN(2,6) = -Eta * (1.0d0 + Xi)
      dN(2,5) = -0.5d0 * (1.0d0 - Xi**2)
      dN(2,4) =  0.25d0 * (1.0d0 - Xi) - 0.5d0 * dN(2,7) - 0.5d0 * dN(2,8)
      dN(2,3) =  0.25d0 * (1.0d0 + Xi) - 0.5d0 * dN(2,6) - 0.5d0 * dN(2,7)
      dN(2,2) = -0.25d0 * (1.0d0 + Xi) - 0.5d0 * dN(2,5) - 0.5d0 * dN(2,6)
      dN(2,1) = -0.25d0 * (1.0d0 - Xi) - 0.5d0 * dN(2,5) - 0.5d0 * dN(2,8)
  end select
  
  ! Calculating the Jacobian Matrix
  J = matmul(dN,E%coor(1:E%nNodes,:)) ! E%coor contains node numbers and their xy coordinates
  
  ! Calculating the determinent of the Jacobian Matrix
  detJ = J(1,1) * J(2,2) - J(2,1) * J(1,2)
  
  ! Calculating the inverse of the Jacobian Matrix
  invJ = (1.0d0 / detJ) * reshape((/J(2,2), -J(2,1), -J(1,2), J(1,1)/), (/2, 2/))
  
  ! Fill the B1 matrix, which relates displacement derivatives to strains
  B1 = 0.0d0
  B1(1,1) = 1.0d0
  B1(2,4) = 1.0d0
  B1(3,2) = 1.0d0
  B1(3,3) = 1.0d0
  
  ! Fill the inverse Jacobian into B2
  B2 = 0.0d0
  B2(1:2,1:2) = invJ
  B2(3:4,3:4) = invJ
  
  ! Fill derivatives of shape functions w.r.t. Xi and Eta into B3
  B3 = 0.0d0
  B3(1,1:(2*E%nNodes)-1:2) = dN(1,:) ! Row 1, odd columns
  B3(2,1:(2*E%nNodes)-1:2) = dN(2,:) ! Row 2, odd columns
  B3(3,2:2*E%nNodes:2) = dN(1,:) ! Row 3, even columns
  B3(4,2:2*E%nNodes:2) = dN(2,:) ! Row 4, even columns
  
  ! Calculating the B matrix by collecting the matrices of equation 3.12, 3.13 and 3.14 from the project report
  B = matmul(B1,matmul(B2,B3))
  
  return
    
end subroutine StrainDisplacementJacobian
  
! -------------------------------------------------------------------
! Subroutine to calculate ConstitutiveMatrix(E, C)
!
! This subroutine calculates the Constitutive matrix
!
! Author: DMS1-1 2025
! -------------------------------------------------------------------

subroutine ConstitutiveMatrix(E, C)

  implicit none  
  
  ! Defining global variables
  type(ElemType), intent(in)                :: E
  real(double), dimension(3,3), intent(out) :: C
  
  ! Defining local variables
  real(double)                              :: Young, Nu, factor
  
  ! Assign material properties to local variables
  Young = E%EX
  Nu = E%NUXY
  
  ! Calculate the factor used in the constitutive matrix
  factor = Young / (1.0d0 - Nu**2)
  
  ! Fill in the constitutive matrix C for plane stress condition
  C = 0.0d0
  C(1,1) = factor
  C(1,2) = factor * Nu
  C(2,1) = factor * Nu
  C(2,2) = factor
  C(3,3) = factor * (1.0d0 - Nu) / 2.0d0
  
  return
  
end subroutine ConstitutiveMatrix
  
! -------------------------------------------------------------------
! Subroutine ElemStiffnessMatrix(E, K)
! 
! This subroutine computes the element stiffness matrix for a 2D solid isoparametric element,
! for example using 4-node 2D solid isoparametric elements.
!
! Input to the subroutine:
!   E: the element data
!
! Output:
!   K : the element stiffness matrix K of dimension 2*E%nNodes x 2*E%nNodes
! 
! Called subroutines:
!   ConstitutiveMatrix
!   StrainDisplacementJacobian
!
! Document:   "DMS7FEA: A Static FE Analysis Program Using Isoparametric 2D Elements" by Erik Lund
! Programmer: EL, 4/10/99
! Revision    date    engineer   review
!   0A:       4/10/99 EL         Initial version
!   0B:       4/11/25 DMS1-1     Added all called subroutines and all calculations
!  -------------------------------------------------------------------

subroutine ElemStiffnessMatrix(E, K)

  implicit none

  ! Global variables
  type(ElemType), intent(in)                    :: E
  real(double), dimension(:,:), intent(in out)  :: K

  ! Local variables
  real(double)                                  :: Xi, Eta, detJ
  real(double), dimension(3,2*E%nNodes)         :: B
  real(double), dimension(3,3)                  :: C
  real(double), allocatable, dimension(:)       :: W, XiGauss, EtaGauss
  integer                                       :: i, j
  
  ! Get the constitutive matrix C
  call ConstitutiveMatrix(E, C)
  
  select case (E%nNodes)
    case (4)
      ! Define the size of XiGauss, EtaGauss and W for 2x2 Gauss quadrature
      allocate(XiGauss(2), EtaGauss(2), W(2))
      
      ! Define the weight factor for 2x2 Gauss quadrature
      W = (/ 1.0d0, 1.0d0 /)
  
      ! Define the 2x2 Gauss Points
      XiGauss = (/ -1.0d0/sqrt(3.0d0), 1.0d0/sqrt(3.0d0) /)
      EtaGauss = (/ -1.0d0/sqrt(3.0d0), 1.0d0/sqrt(3.0d0) /)
      
    case (8)
      ! Define the size of XiGauss, EtaGauss and W for 3x3 Gauss quadrature
      allocate(XiGauss(3), EtaGauss(3), W(3))
      
      ! Define the weight factor for 3x3 Gauss quadrature
      W = (/ 5.0d0/9.0d0, 8.0d0/9.0d0, 5.0d0/9.0d0 /)
      
      ! Define the 3x3 Gauss Points
      XiGauss = (/ -sqrt(3.0d0/5.0d0), 0.0d0, sqrt(3.0d0/5.0d0) /)
      EtaGauss = (/ -sqrt(3.0d0/5.0d0), 0.0d0, sqrt(3.0d0/5.0d0) /)      
  end select
  
  ! Initialize K
  K = 0.0d0  
  
  ! Compute the element stiffness matrix K
  do i = 1,size(XiGauss)
    ! Get the Xi coordinate of the Gauss point
    Xi = XiGauss(i)
    
    do j = 1,size(EtaGauss)
      ! Get the Eta coordinate of the Gauss point
      Eta = EtaGauss(j)
      
      ! Compute B and detJ for the Gauss point
      call StrainDisplacementJacobian(E, Xi, Eta, B, detJ)
      
      ! Integrate stiffness matrix (E%t is the element thickness)
      K = K + matmul(transpose(B), matmul(C, B)) * W(i) * W(j) * detJ * E%t
    end do
  end do
  
  deallocate(XiGauss, EtaGauss, W)
  
  return
  
end subroutine ElemStiffnessMatrix

! -------------------------------------------------------------------
! Subroutine ElemStresses(E, ElemD, ElemStress)
! 
! This subroutine computes the element stresses for 2D solid isoparametric elements,
! for example by using 4-node 2D solid isoparametric elements
!
! Input to the subroutine:
!   E         : the element data
!   ElemD     : the element displacement vector of dimension 2*E%nNodes
!
! Output:
!   NodeStress: the node stresses of dimension 1:E%nNodes x 10 (SX, SY, TXY, SRef, S1, S2, ...)
!               Please note that NodeStress is defined statically to be of dimension 9 x 10
!   ElemStress: the element stresses of dimension 10 (SX, SY, TXY, SRef, ...)
!               The node/element stresses are stored in the order: SX, SY, TXY, SRef, S1, S2, .. (max 10 values)
!               For ElemStress, these other default positions are assumed in IO.f90:
!               - ElemStress(7): the principal stress direction
!               - ElemStress(8): element specific strain energy density
!               If you decide to store values at other positions, then remember to add I/O routines 
!               in WriteOutputToFEPlot in IO.f90 for these extra data.
! 
! Called subroutines:
!   ConstitutiveMatrix
!   StrainDisplacementJacobian
!
! Document:   "DMS7FEA: A Static FE Analysis Program Using Isoparametric 2D Elements" by Erik Lund
! Programmer: EL, 4/10/99
! Revision    date    engineer   review
!   0A:       4/10/99 EL         Initial version
!   0B:       1/12/23 EL         Updated the header
!   0C:       4/11/25 DMS1-1     Added all called subroutines and all calculations 
!  -------------------------------------------------------------------

subroutine ElemStresses(E, ElemD, NodeStress, ElemStress)

  implicit none

  ! Global variables
  type(ElemType), intent(in)                    :: E
  real(double), dimension(:), intent(in)        :: ElemD
  real(double), dimension(:,:), intent(in out)  :: NodeStress
  real(double), dimension(:), intent(in out)    :: ElemStress
  
  ! Local variables
  real(double), dimension(3,3)                  :: C
  real(double), dimension(3,2*E%nNodes)         :: B
  real(double), dimension(E%nNodes)             :: rCoor, sCoor
  real(double), dimension(4)                    :: N
  real(double), dimension(4,2)                  :: GaussPoints
  real(double), dimension(4,3)                  :: GaussPointStress
  real(double)                                  :: Xi, Eta, detJ, r, s
  integer                                       :: i
  
  ! Get the constitutive matrix C
  call ConstitutiveMatrix(E, C)
      
  ! Defining r and s coordinates of nodes
  select case (E%nNodes)
    case (4)
      rCoor = (/-sqrt(3.0d0),  sqrt(3.0d0), sqrt(3.0d0), -sqrt(3.0d0)/)
      sCoor = (/-sqrt(3.0d0), -sqrt(3.0d0), sqrt(3.0d0),  sqrt(3.0d0)/)
    
    case (8)
      rCoor = (/-sqrt(3.0d0),  sqrt(3.0d0), sqrt(3.0d0), -sqrt(3.0d0), &
                 0.0d0,        sqrt(3.0d0), 0.0d0,       -sqrt(3.0d0)/)
      sCoor = (/-sqrt(3.0d0), -sqrt(3.0d0), sqrt(3.0d0),  sqrt(3.0d0), &
                -sqrt(3.0d0),  0.0d0,       sqrt(3.0d0),  0.0d0      /)
  end select
  
  ! Define the 2x2 Gauss Points
  GaussPoints = reshape((/ &
    -1.0d0/sqrt(3.0d0), -1.0d0/sqrt(3.0d0), &
     1.0d0/sqrt(3.0d0), -1.0d0/sqrt(3.0d0), &
     1.0d0/sqrt(3.0d0),  1.0d0/sqrt(3.0d0), &
    -1.0d0/sqrt(3.0d0),  1.0d0/sqrt(3.0d0)  & 
      /), (/4,2/), order=(/2,1/) )
  
  ! Calculate stresses at Gauss points
  do i = 1, 4
    ! Get the Xi and Eta coordinates of the Gauss point
    Xi = GaussPoints(i,1)
    Eta = GaussPoints(i,2)
    
    ! Compute B for the Gauss point
    call StrainDisplacementJacobian(E, Xi, Eta, B, detJ)
    
    ! Compute the stress at the Gauss point
    GaussPointStress(i,1:3) = matmul(C, matmul(B, ElemD(1:2*E%nNodes)))
  end do
  
  ! Initialize NodeStress to zero
  NodeStress = 0.0d0
  
  ! Calculate node stresses
  do i = 1,E%nNodes
    ! Initialize r and s coordinates of the node
    r = rCoor(i)
    s = sCoor(i)

    ! Define the shape functions
    N(1) = 0.25d0 * (1.0d0 - r) * (1.0d0 - s)
    N(2) = 0.25d0 * (1.0d0 + r) * (1.0d0 - s)
    N(3) = 0.25d0 * (1.0d0 + r) * (1.0d0 + s)
    N(4) = 0.25d0 * (1.0d0 - r) * (1.0d0 + s)

    ! Calculate node stresses by extrapolation from Gauss point stresses
    NodeStress(i,1:3) = matmul(transpose(GaussPointStress), N)

    ! Calculate node reference stress
    NodeStress(i,4) = sqrt(NodeStress(i,1)**2 - NodeStress(i,1) * NodeStress(i,2) + NodeStress(i,2)**2 + 3 * NodeStress(i,3)**2)

    ! Calculate node principal stresses and principal stress direction
    NodeStress(i,5) = 0.5d0 * (NodeStress(i,1) + NodeStress(i,2)) + &
	    sqrt( ( (NodeStress(i,1) - NodeStress(i,2)) /2 )**2 + NodeStress(i,3)**2 )
    NodeStress(i,6) = 0.5d0 * (NodeStress(i,1) + NodeStress(i,2)) - &
	    sqrt( ( (NodeStress(i,1) - NodeStress(i,2)) /2 )**2 + NodeStress(i,3)**2 )
    NodeStress(i,7) = atan2(2.0d0 * NodeStress(i,3), NodeStress(i,1) - NodeStress(i,2)) / 2.0d0
  end do
  
  ! Calculate element stress at origin of natural coordinate system
  Xi = 0.0d0
  Eta = 0.0d0
  call StrainDisplacementJacobian(E, Xi, Eta, B, detJ)
  ElemStress(1:3) = matmul(C, matmul(B, ElemD(1:2*E%nNodes)))
  
  ! Calculate element reference stress
  ElemStress(4) = sqrt(ElemStress(1)**2 - ElemStress(1) * ElemStress(2) + ElemStress(2)**2 + 3 * ElemStress(3)**2)
  
  ! Calculate element principal stresses and principal stress direction
  ElemStress(5) = 0.5d0 * (ElemStress(1) + ElemStress(2)) + &
    sqrt( ( (ElemStress(1) - ElemStress(2)) /2 )**2 + ElemStress(3)**2 )
  ElemStress(6) = 0.5d0 * (ElemStress(1) + ElemStress(2)) - &
    sqrt( ( (ElemStress(1) - ElemStress(2)) /2 )**2 + ElemStress(3)**2 )
  ElemStress(7) = atan2(2.0d0 * ElemStress(3), ElemStress(1) - ElemStress(2)) / 2.0d0
  
  return

end subroutine ElemStresses

! -------------------------------------------------------------------
! Subroutine GetElemData(G, E)
! 
! This subroutine copies FE info from G to the element.
!
! Input to the subroutine:
!   G     : all information about the model
!   ElemNo: the element number
!
! Output:
!   E     : all data for the element
! 
! Called subroutines:
!   none
!
! Document:   "DMS7FEA: A Static FE Analysis Program Using Isoparametric 2D Elements" by Erik Lund
! Programmer: EL, 4/10/99
! Revision    date    engineer   review
!   0A:       4/10/99 EL         Initial version
!   0B:
!  -------------------------------------------------------------------

subroutine GetElemData(G, ElemNo, E)

  implicit none

  type(GlobalFEData), intent(in)   :: G
  integer, intent(in)              :: ElemNo
  type(ElemType), intent(in out)   :: E

  integer :: NodeNo, DimNo, i

  ! Copy the element number
  E%No = ElemNo

  ! Set the number of nodes for the element
  E%nNodes = G%nNodesPerElem(ElemNo)

  ! The element type can be 0: plane stress or 2: plane strain
  E%Type = G%EGroup( G%ElemInfo(ElemNo)%EGroupNo )%Type
  
  ! Copy material properties EX, NUXY and thickness t
  E%EX   = G%MProp( G%ElemInfo(ElemNo)%MPropNo )%EX
  E%NUXY = G%MProp( G%ElemInfo(ElemNo)%MPropNo )%NUXY
  E%t    = G%RConst( G%ElemInfo(ElemNo)%RConstNo )%t

  ! Copy list of nodes defining the element to Elem%Nodes
  NodesLoop: do NodeNo = 1,E%nNodes
    E%Nodes(NodeNo) = G%ElemConnect(ElemNo,NodeNo)
  end do NodesLoop

  ! Copy coordinates of the nodes to E%Coor(NodeNo,DimNo)
  CoorLoop: do NodeNo = 1,E%nNodes
    do DimNo = 1,2
      E%Coor(NodeNo,DimNo) = G%Coor(E%Nodes(NodeNo),DimNo)
    end do
  end do CoorLoop

  ! Copy DOF numbers to E%DOF
  i = 0
  DOFLoop: do NodeNo = 1,E%nNodes
    do DimNo = 1,2
      i = i+1
      E%DOF(i) = G%DOF(E%Nodes(NodeNo),DimNo)
    end do
  end do DOFLoop

  return

end subroutine GetElemData


! -------------------------------------------------------------------
! Subroutine GetElemDisp(A, E, ElemD)
! 
! This subroutine displacements stored in A%Uf and 
! A%Up to element displacement vector ElemD
!
! Input to the subroutine:
!   A    : analysis info
!   E    : all data for the element
!
! Output:
!   ElemD: the element displacement vector of dimension 2*E%nNodes
! 
! Called subroutines:
!   none
!
! Document:   "DMS7FEA: A Static FE Analysis Program Using Isoparametric 2D Elements" by Erik Lund
! Programmer: EL, 4/10/99
! Revision    date    engineer   review
!   0A:       4/10/99 EL         Initial version
!   0B:
!  -------------------------------------------------------------------

subroutine GetElemDisp(A, E, ElemD)

  implicit none

  type(AnalysisData), intent(in)             :: A
  type(ElemType), intent(in)                 :: E
  real(double), dimension(:), intent(in out) :: ElemD

  integer :: DOFNo

  ! Copy displacements stored in A%Uf and A%Up
  ! to element displacement vector ElemD
  DOFNoLoop: do DOFNo = 1, 2*E%nNodes
    if (E%DOF(DOFNo) > 0) then
      ElemD(DOFNo) = A%Uf(E%DOF(DOFNo))
    else
      ElemD(DOFNo) = A%Up(-E%DOF(DOFNo))
    end if
  end do DOFNoLoop


  return

end subroutine GetElemDisp



! -------------------------------------------------------------------
! Subroutine ElemLoadVector(G, E)
! 
! This subroutine computes the consistent element load vector ElemR
! according to pressure specifications in G%PEl
! The nodal load calculation is only implemented for linear elements,
! so if you are going to implement quadratic elements, then you need
! to add code in this subroutine - it's not so difficult :-)
!
! Input to the subroutine:
!   PEl  : type containg element pressure specifications
!   E    : the element data
!
! Output:
!   ElemR: the consistent element load vector of length 2*E%nNodes
! 
! Called subroutines:
!   none
!
! Document:   "DMS7FEA: A Static FE Analysis Program Using Isoparametric 2D Elements" by Erik Lund
! Programmer: EL, 4/10/99
! Revision    date    engineer   review
!   0A:       4/10/99 EL         Initial version
!   0B:       17/9/02 EL         Updated with linearly varying loads
!   0C:       3/12/04 EL         A few comments added (together with FaceNodes(3))
!  -------------------------------------------------------------------

subroutine ElemLoadVector(PEl, E, ElemR)

  implicit none

  type(PElType), intent(in)                  :: PEl
  type(ElemType), intent(in)                 :: E
  real(double), dimension(:), intent(in out) :: ElemR

  integer      :: FaceNodes(3), NodeNo, DimNo, DOFNo
  real(double) :: Length, NodalForce(2), FaceVec(2), FaceNormalVec(2)

  ! In case of the bilinear Q4 element with constant pressure loading,
  ! the loads are distributed evenly between the two face nodes

  ! Initialize ElemR
  ElemR(:) = 0.0

  ! First, we find the two nodes for which the pressure apply
  ! PEl%FaceNo is the face number (1-4 for a quadrilateral element)
  if (PEl%FaceNo <= 3) then
    FaceNodes(1) = PEl%FaceNo
    FaceNodes(2) = PEl%FaceNo+1
  else
    FaceNodes(1) = PEl%FaceNo
    FaceNodes(2) = 1
  end if
  ! Take care of quadratic elements
  if (E%nNodes == 8) FaceNodes(3) = PEl%FaceNo+4
  if (E%nNodes == 6) FaceNodes(3) = PEl%FaceNo+3

  ! The vector going from node 1 to node 2 on this face (boundary) is called FaceVec
  SetFaceVecLoop: do DimNo = 1,2
    FaceVec(DimNo) = E%Coor(FaceNodes(2),DimNo) - E%Coor(FaceNodes(1),DimNo)
  end do SetFaceVecLoop

  ! The lenght of the face is called Length
  Length = sqrt(FaceVec(1)**2 + FaceVec(2)**2)

  ! First deal with directions X: 1 or Y: 2
  if (PEl%Dir == 1 .or. PEl%Dir == 2) then
    ! Distribute the loads in a consistent way
    select case (E%nNodes)

      case (3,4)
        ! The pressure is multiplied by face area and distributed between the two nodes
        NodalForce(1) = Length*(2.0*E%t*PEl%Value(1) + E%t*PEl%Value(2))/6.0
        NodalForce(2) = Length*(2.0*E%t*PEl%Value(2) + E%t*PEl%Value(1))/6.0
        ! Add NodalForce to the face nodes in the direction PEl%Dir
        do NodeNo = 1,2
          DOFNo = (FaceNodes(NodeNo)-1)*2 + PEl%Dir
          ElemR(DOFNo) = ElemR(DOFNo) + NodalForce(NodeNo)
        end do

      ! Add code for quadratic elements
      case (6,8,9)

    end select
  end if


  ! Finally, we deal with pressure normal to the face
  if (PEl%Dir == 4) then
  
    ! First, compute a normal vector to the face
    FaceNormalVec(1) = -FaceVec(2)/Length
    FaceNormalVec(2) = FaceVec(1)/Length

    ! Distribute the loads in a consistent way
    select case (E%nNodes)

      case (3,4)
        ! The pressure is multiplied by face area and distributed between the 2 nodes
        NodalForce(1) = Length*(2.0*E%t*PEl%Value(1) + E%t*PEl%Value(2))/6.0
        NodalForce(2) = Length*(2.0*E%t*PEl%Value(2) + E%t*PEl%Value(1))/6.0
        ! Loop over the 2 face nodes and add projected nodal forces to ElemR
        do NodeNo = 1,2
          do DimNo = 1,2
            DOFNo = (FaceNodes(NodeNo)-1)*2 + DimNo
            ElemR(DOFNo) = ElemR(DOFNo) + NodalForce(NodeNo)*FaceNormalVec(DimNo)
          end do
        end do

      ! Add code for quadratic elements
      case (6,8,9)

    end select

  end if

  return

end subroutine ElemLoadVector



end module Elem
