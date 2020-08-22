program main
    implicit none
    real*4 :: beta=1.72
    integer*4 :: windows=5
    integer*4 :: frames=10000
    character(len=50) :: filename
    character(len=80) :: line
    integer*4 :: i,j,k
    real*4,allocatable :: energy(:,:,:)
    allocate(energy(windows,frames,windows))
    write(*,*)"*****************************************"
    write(*,*)"Read the data from outfile"
    open(101,file='filename.dat')
    do i=1,windows
        j=0
        read(101,*),filename
        open(102,file=filename)
        do
        read(102,'(a80)')line
        if (line(7:15).eq.'A V E R A') exit
        if (line(1:4).eq.'MBAR')then
            j=j+1    
            do k=1,windows
                read(102,'(a80)')line
                if (line(20:20).eq.'*')then
                    energy(i,j,k)=10000000
                else
                    read(line(21:31),*)energy(i,j,k)
                end if
            enddo
        end if
        enddo
        write(*,*)j
        close(102)
    enddo
    close(101)
call FEPcalculate(beta,energy,frames,windows)
call Barcalculate(beta,energy,frames,windows)
end program main
subroutine Barcalculate(beta,energy,frames,windows)
    implicit none
    integer*4 :: i,j,k
    integer*4 :: frames,windows
    real*4 :: beta
    real*4 :: energy(windows,frames,windows)
    real*4 :: A0_1(frames)
    real*4 :: A1_0(frames)
    real*4 :: Ave_A0_1,Ave_A1_0
    real*4 :: C
    real*4 :: Delta_A
    real*4 :: F0_1(frames)
    real*4 :: F1_0(frames)
    real*4 :: Sum_F0_1,Sum_F1_0,Ave_F0_1,Ave_F1_0
    real*4 :: SumDelta,Ave_error10,Ave_error01
    real*4 :: F0_1_error,F1_0_error,Error,SumError
    SumDelta=0.0
    SumError=0.0
    write(*,*)"***********************************************"
    write(*,*)"The module of Bar:"
    do i=1,windows-1
        Ave_F0_1=0.0
        Ave_F1_0=0.0
        F0_1_error=0.0
        F1_0_error=0.0
        C=0.0
        do j=1,frames
            A0_1(j)=energy(i+1,j,i)-energy(i+1,j,i+1)
            A1_0(j)=energy(i,j,i+1)-energy(i,j,i)
        enddo
        do
        do j=1,frames
            F0_1(j)=1.0/(1.0+exp(beta*(A0_1(j)+C)))
            F1_0(j)=1.0/(1.0+exp(beta*(A1_0(j)-C)))
            F0_1_error=F0_1_error+(1.0/(1.0+exp(beta*(A0_1(j)+C))))**2
            F1_0_error=F1_0_error+(1.0/(1.0+exp(beta*(A1_0(j)-C))))**2
            Ave_F0_1=Ave_F0_1+F0_1(j)
            Ave_F1_0=Ave_F1_0+F1_0(j)
        end do
        Ave_F0_1=Ave_F0_1/real(frames)
        Ave_F1_0=Ave_F1_0/real(frames)
        Ave_error10=F1_0_error/real(frames)
        Ave_error01=F0_1_error/real(frames)
        if (abs(Ave_F0_1-Ave_F1_0).lt.1E-6)then
            Delta_A=C
            Error=(Ave_error10/(Ave_F1_0**2)+Ave_error01/(Ave_F0_1**2)-2.0)/real(frames)
            exit
        else
            Delta_A=log(Ave_F0_1/Ave_F1_0)/beta+C
            C=Delta_A
        end if
        enddo
        SumDelta=SumDelta+Delta_A
        SumError=SumError+Error
        write(*,'(f9.4,5x,f9.4)')Delta_A,Error
    enddo
    write(*,*)"Bar_value with deviation (kcal/mol):"
    write(*,'(f9.4,5x,f9.4)')SumDelta,SumError
end subroutine Barcalculate
subroutine FEPcalculate(beta,energy,frames,windows)
    implicit none
    integer*4 :: i,j,k
    integer*4 :: frames
    integer*4 :: windows
    real*4 :: beta
    real*4 :: energy(windows,frames,windows)
    real*4 :: delta_E(frames)
    real*4 :: List_E(windows)
    real*4 :: List_Std(windows)
    real*4 :: Stdsum_E
    real*4 :: Std_E
    real*4 :: Exp_E
    real*4 :: Sum_E
    real*4 :: Ave_delta_E
    real*4 :: Log_E
    real*8 :: Fep_value
    real*8 :: Fep_std
    write(*,*)"*****************************************"
    write(*,*)"The module of FEP:"
    do i=1,windows-1
        Sum_E=0.0
        Stdsum_E=0.0
        do j=1,frames
            delta_E(j)=energy(i,j,i+1)-energy(i,j,i)
            Exp_E=exp(-beta*delta_E(j))
            Sum_E=Sum_E+Exp_E
            Stdsum_E=Stdsum_E+(Exp_E-Ave_delta_E)**2
        enddo
        Stdsum_E=Stdsum_E/REAL(frames)
        Ave_delta_E=Sum_E/real(frames)
        Std_E=sqrt(Stdsum_E/real(frames))/(beta*Ave_delta_E)
        Log_E=-log(Ave_delta_E)/beta
        List_E(i)=Log_E
        List_Std(i)=Std_E
        write(*,'(f9.4,f9.4)')List_E(i),List_Std(i)
    enddo
    Fep_value=0
    do i=1,windows-1
        Fep_value=Fep_value+List_E(i)
        Fep_std=Fep_std+List_Std(i)
    end do
    write(*,*)"Fep_value with deviation (kcal/mol):"
    write(*,'(f9.4,5x,f6.4)')Fep_value,Fep_std
end subroutine FEPcalculate
