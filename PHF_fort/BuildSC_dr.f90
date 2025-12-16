    Module BuildSC_dr
	Use Precision
	Use Constants
	Implicit None

	Contains

	Subroutine BuildSC_Vec(c0,c1,c2,t0,t1,t2,Y,NSO)
	Implicit None
	Integer,            Intent(in)  :: NSO
	Complex(kind=pr),   Intent(in)  :: c0, c1(NSO,NSO)
	Complex(kind=pr),   Intent(in)  :: c2(NSO,NSO,NSO,NSO)
	Complex(kind=pr),   Intent(in)  :: Y(NSO,NSO)
	Complex(kind=pr),   Intent(out) :: t0, t1(NSO,NSO)
	Complex(kind=pr),   Intent(out) :: t2(NSO,NSO,NSO,NSO)
	complex(kind=pr), dimension(:, :, :, :), allocatable :: tau0

    complex(kind=pr), dimension(:, :), allocatable :: tau1

    complex(kind=pr), dimension(:, :), allocatable :: tau2

    complex(kind=pr) :: tau3

    complex(kind=pr), dimension(:, :, :, :), allocatable :: tau4

    complex(kind=pr), dimension(:, :), allocatable :: tau5

    complex(kind=pr), dimension(:, :), allocatable :: tau6

    complex(kind=pr), dimension(:, :), allocatable :: tau7

    complex(kind=pr) :: tau8

    complex(kind=pr), dimension(:, :, :, :), allocatable :: tau9

    complex(kind=pr), dimension(:, :, :, :), allocatable :: tau10

    complex(kind=pr), dimension(:, :, :, :), allocatable :: tau11

    complex(kind=pr), dimension(:, :, :, :), allocatable :: tau12

    complex(kind=pr), dimension(:, :), allocatable :: tau13

    complex(kind=pr), dimension(:, :), allocatable :: tau14

    complex(kind=pr), dimension(:, :), allocatable :: tau15

    complex(kind=pr), dimension(:, :, :, :), allocatable :: tau16

    complex(kind=pr), dimension(:, :), allocatable :: tau17

    complex(kind=pr), dimension(:, :), allocatable :: tau18

    complex(kind=pr), dimension(:, :, :, :), allocatable :: tau19

    complex(kind=pr), dimension(:, :, :, :), allocatable :: tau20

    complex(kind=pr), dimension(:, :), allocatable :: tau21

    complex(kind=pr), dimension(:, :, :, :), allocatable :: tau22

    complex(kind=pr) :: tau23

    complex(kind=pr), dimension(:, :), allocatable :: tau24

    complex(kind=pr), dimension(:, :), allocatable :: tau25

    complex(kind=pr), dimension(:, :), allocatable :: tau26

    complex(kind=pr), dimension(:, :), allocatable :: tau27
	Integer :: a, b, c, d, e, f

	allocate(tau0(nso, nso, nso, nso))
    allocate(tau1(nso, nso))
    allocate(tau2(nso, nso))
    allocate(tau6(nso, nso))
    allocate(tau13(nso, nso))
    allocate(tau5(nso, nso))
    allocate(tau14(nso, nso))
    allocate(tau15(nso, nso))
    allocate(tau16(nso, nso, nso, nso))
    allocate(tau21(nso, nso))
    allocate(tau22(nso, nso, nso, nso))
    allocate(tau4(nso, nso, nso, nso))
    allocate(tau7(nso, nso))
    allocate(tau9(nso, nso, nso, nso))
    allocate(tau10(nso, nso, nso, nso))
    allocate(tau11(nso, nso, nso, nso))
    allocate(tau12(nso, nso, nso, nso))
    allocate(tau19(nso, nso, nso, nso))
    allocate(tau20(nso, nso, nso, nso))
    allocate(tau17(nso, nso))
    allocate(tau18(nso, nso))
    allocate(tau27(nso, nso))
    allocate(tau24(nso, nso))
    allocate(tau25(nso, nso))
    allocate(tau26(nso, nso))

    !$omp parallel default(shared)

    !$omp single
    tau0 = 0.0
    !$omp end single

    !$omp do schedule(static)
    do d=1, nso
        do c=1, nso
            do b=1, nso
                do a=1, nso
    
                    tau0(a, b, c, d) = tau0(a, b, c, d) + ( &
                        c2(a, b, c, d)&
                    )
    
                end do
            end do
        end do
    end do
    !$omp end do

    !$omp do schedule(static)
    do d=1, nso
        do c=1, nso
            do b=1, nso
                do a=1, nso
    
                    tau0(a, b, c, d) = tau0(a, b, c, d) - ( &
                        c2(b, a, c, d)&
                    )
    
                end do
            end do
        end do
    end do
    !$omp end do

    !$omp single
    tau1 = 0.0
    !$omp end single

    !$omp do schedule(static)
    do b=1, nso
        do a=1, nso
            do d=1, nso
                do c=1, nso
                    tau1(a, b) = tau1(a, b) + ( &
                        Y(c, d) * tau0(d, a, b, c)&
                    )
                end do
            end do
        end do
    end do
    !$omp end do

    !$omp single
    tau2 = 0.0
    !$omp end single

    !$omp do schedule(static)
    do b=1, nso
        do a=1, nso
    
            tau2(a, b) = tau2(a, b) - ( &
                tau1(a, b)&
            )
    
        end do
    end do
    !$omp end do

    !$omp single
    tau6 = 0.0
    !$omp end single

    !$omp do schedule(static)
    do b=1, nso
        do a=1, nso
    
            tau6(a, b) = tau6(a, b) - ( &
                tau1(a, b)&
            )
    
        end do
    end do
    !$omp end do

    !$omp single
    tau8 = 0.0
    !$omp end single

    !$omp do schedule(static) reduction(+:tau8)
    
    do b=1, nso
        do a=1, nso
            tau8 = tau8 - ( &
                Y(a, b) * tau1(b, a)&
            )
        end do
    end do
    
    !$omp end do

    !$omp single
    t2 = 0.0
    !$omp end single

    !$omp do schedule(static)
    do d=1, nso
        do c=1, nso
            do b=1, nso
                do a=1, nso
    
                    t2(a, b, c, d) = t2(a, b, c, d) + ( &
                        tau8 * Y(a, c) * Y(b, d)&
                    )
    
                end do
            end do
        end do
    end do
    !$omp end do

    !$omp do schedule(static)
    do d=1, nso
        do c=1, nso
            do b=1, nso
                do a=1, nso
    
                    t2(a, b, c, d) = t2(a, b, c, d) - ( &
                        tau8 * Y(a, d) * Y(b, c)&
                    )
    
                end do
            end do
        end do
    end do
    !$omp end do

    !$omp single
    tau13 = 0.0
    !$omp end single

    !$omp do schedule(static)
    do b=1, nso
        do a=1, nso
    
            tau13(a, b) = tau13(a, b) - ( &
                tau1(a, b)&
            )
    
        end do
    end do
    !$omp end do


    !$omp single
    tau5 = 0.0
    !$omp end single

    !$omp do schedule(static)
    do b=1, nso
        do a=1, nso
            do d=1, nso
                do c=1, nso
                    tau5(a, b) = tau5(a, b) + ( &
                        Y(c, d) * tau0(a, d, c, b)&
                    )
                end do
            end do
        end do
    end do
    !$omp end do


    !$omp do schedule(static)
    do b=1, nso
        do a=1, nso
    
            tau6(a, b) = tau6(a, b) - ( &
                tau5(a, b)&
            )
    
        end do
    end do
    !$omp end do

    !$omp do schedule(static)
    do b=1, nso
        do a=1, nso
    
            tau13(a, b) = tau13(a, b) - ( &
                tau5(a, b)&
            )
    
        end do
    end do
    !$omp end do


    !$omp single
    tau14 = 0.0
    !$omp end single

    !$omp do schedule(static)
    do b=1, nso
        do a=1, nso
            do c=1, nso
                tau14(a, b) = tau14(a, b) + ( &
                    Y(a, c) * tau13(c, b)&
                )
            end do
        end do
    end do
    !$omp end do

    !$omp single
    tau15 = 0.0
    !$omp end single

    !$omp do schedule(static)
    do b=1, nso
        do a=1, nso
            do c=1, nso
                tau15(a, b) = tau15(a, b) + ( &
                    Y(c, a) * tau14(b, c)&
                )
            end do
        end do
    end do
    !$omp end do


    !$omp single
    tau16 = 0.0
    !$omp end single

    !$omp do schedule(static)
    do d=1, nso
        do c=1, nso
            do b=1, nso
                do a=1, nso
    
                    tau16(a, b, c, d) = tau16(a, b, c, d) - ( &
                        Y(a, c) * tau15(d, b)&
                    )
    
                end do
            end do
        end do
    end do
    !$omp end do


    !$omp single
    tau21 = 0.0
    !$omp end single

    !$omp do schedule(static)
    do b=1, nso
        do a=1, nso
            do c=1, nso
                tau21(a, b) = tau21(a, b) + ( &
                    Y(c, a) * tau13(b, c)&
                )
            end do
        end do
    end do
    !$omp end do


    !$omp single
    tau22 = 0.0
    !$omp end single

    !$omp do schedule(static)
    do d=1, nso
        do c=1, nso
            do b=1, nso
                do a=1, nso
    
                    tau22(a, b, c, d) = tau22(a, b, c, d) + ( &
                        Y(a, b) * tau21(c, d)&
                    )
    
                end do
            end do
        end do
    end do
    !$omp end do


    !$omp do schedule(static)
    do b=1, nso
        do a=1, nso
    
            tau2(a, b) = tau2(a, b) + ( &
                c1(a, b)&
            )
    
        end do
    end do
    !$omp end do

    !$omp single
    tau3 = 0.0
    !$omp end single

    !$omp do schedule(static) reduction(+:tau3)
    
    do b=1, nso
        do a=1, nso
            tau3 = tau3 + ( &
                Y(a, b) * tau2(b, a)&
            )
        end do
    end do
    
    !$omp end do


    !$omp single
    t0 = 0.0
    !$omp end single

    !$omp single
    
    
    t0 = t0 + ( &
        tau3&
    )
    
    
    !$omp end single

    !$omp single
    t1 = 0.0
    !$omp end single

    !$omp do schedule(static)
    do b=1, nso
        do a=1, nso
    
            t1(a, b) = t1(a, b) + ( &
                tau3 * Y(a, b)&
            )
    
        end do
    end do
    !$omp end do

    !$omp single
    tau4 = 0.0
    !$omp end single

    !$omp do schedule(static)
    do d=1, nso
        do c=1, nso
            do b=1, nso
                do a=1, nso
    
                    tau4(a, b, c, d) = tau4(a, b, c, d) - ( &
                        c2(a, b, c, d)&
                    )
    
                end do
            end do
        end do
    end do
    !$omp end do

    !$omp do schedule(static)
    do d=1, nso
        do c=1, nso
            do b=1, nso
                do a=1, nso
    
                    tau4(a, b, c, d) = tau4(a, b, c, d) + ( &
                        c2(b, a, c, d)&
                    )
    
                end do
            end do
        end do
    end do
    !$omp end do

    !$omp single
    tau7 = 0.0
    !$omp end single

    !$omp do schedule(static)
    do b=1, nso
        do a=1, nso
            do d=1, nso
                do c=1, nso
                    tau7(a, b) = tau7(a, b) - ( &
                        Y(c, d) * tau4(d, a, b, c)&
                    )
                end do
            end do
        end do
    end do
    !$omp end do

    !$omp do schedule(static)
    do b=1, nso
        do a=1, nso
            do d=1, nso
                do c=1, nso
                    tau7(a, b) = tau7(a, b) - ( &
                        Y(c, d) * tau4(a, d, c, b)&
                    )
                end do
            end do
        end do
    end do
    !$omp end do


    !$omp do schedule(static)
    do b=1, nso
        do a=1, nso
    
            tau6(a, b) = tau6(a, b) + ( &
                c1(a, b)&
            )
    
        end do
    end do
    !$omp end do

    !$omp do schedule(static)
    do b=1, nso
        do a=1, nso
            do c=1, nso
                tau7(a, b) = tau7(a, b) + ( &
                    Y(a, c) * tau6(c, b)&
                )
            end do
        end do
    end do
    !$omp end do


    !$omp do schedule(static)
    do b=1, nso
        do a=1, nso
    
            tau7(a, b) = tau7(a, b) - ( &
                c1(a, b)&
            )
    
        end do
    end do
    !$omp end do

    !$omp do schedule(static)
    do b=1, nso
        do a=1, nso
            do c=1, nso
                t1(a, b) = t1(a, b) - ( &
                    Y(c, b) * tau7(a, c)&
                )
            end do
        end do
    end do
    !$omp end do


    !$omp single
    tau9 = 0.0
    !$omp end single

    !$omp do schedule(static)
    do d=1, nso
        do c=1, nso
            do b=1, nso
                do a=1, nso
                    do e=1, nso
                        tau9(a, b, c, d) = tau9(a, b, c, d) + ( &
                            Y(e, a) * c2(b, c, d, e)&
                        )
                    end do
                end do
            end do
        end do
    end do
    !$omp end do

    !$omp single
    tau10 = 0.0
    !$omp end single

    !$omp do schedule(static)
    do d=1, nso
        do c=1, nso
            do b=1, nso
                do a=1, nso
                    do e=1, nso
                        tau10(a, b, c, d) = tau10(a, b, c, d) + ( &
                            Y(e, a) * tau9(b, c, d, e)&
                        )
                    end do
                end do
            end do
        end do
    end do
    !$omp end do


    !$omp single
    tau11 = 0.0
    !$omp end single

    !$omp do schedule(static)
    do d=1, nso
        do c=1, nso
            do b=1, nso
                do a=1, nso
                    do e=1, nso
                        tau11(a, b, c, d) = tau11(a, b, c, d) + ( &
                            Y(a, e) * tau10(b, c, d, e)&
                        )
                    end do
                end do
            end do
        end do
    end do
    !$omp end do

    !$omp single
    tau12 = 0.0
    !$omp end single

    !$omp do schedule(static)
    do d=1, nso
        do c=1, nso
            do b=1, nso
                do a=1, nso
                    do e=1, nso
                        tau12(a, b, c, d) = tau12(a, b, c, d) + ( &
                            Y(a, e) * tau11(b, c, d, e)&
                        )
                    end do
                end do
            end do
        end do
    end do
    !$omp end do


    !$omp do schedule(static)
    do d=1, nso
        do c=1, nso
            do b=1, nso
                do a=1, nso
    
                    tau16(a, b, c, d) = tau16(a, b, c, d) + ( &
                        tau12(a, b, c, d)&
                    )
    
                end do
            end do
        end do
    end do
    !$omp end do


    !$omp do schedule(static)
    do d=1, nso
        do c=1, nso
            do b=1, nso
                do a=1, nso
    
                    t2(a, b, c, d) = t2(a, b, c, d) + ( &
                        tau16(a, b, c, d)&
                    )
    
                end do
            end do
        end do
    end do
    !$omp end do

    !$omp do schedule(static)
    do d=1, nso
        do c=1, nso
            do b=1, nso
                do a=1, nso
    
                    t2(a, b, c, d) = t2(a, b, c, d) - ( &
                        tau16(a, b, d, c)&
                    )
    
                end do
            end do
        end do
    end do
    !$omp end do

    !$omp do schedule(static)
    do d=1, nso
        do c=1, nso
            do b=1, nso
                do a=1, nso
    
                    t2(a, b, c, d) = t2(a, b, c, d) - ( &
                        tau16(b, a, c, d)&
                    )
    
                end do
            end do
        end do
    end do
    !$omp end do

    !$omp do schedule(static)
    do d=1, nso
        do c=1, nso
            do b=1, nso
                do a=1, nso
    
                    t2(a, b, c, d) = t2(a, b, c, d) + ( &
                        tau16(b, a, d, c)&
                    )
    
                end do
            end do
        end do
    end do
    !$omp end do


    !$omp single
    tau19 = 0.0
    !$omp end single

    !$omp do schedule(static)
    do d=1, nso
        do c=1, nso
            do b=1, nso
                do a=1, nso
    
                    tau19(a, b, c, d) = tau19(a, b, c, d) - ( &
                        tau10(a, b, c, d)&
                    )
    
                end do
            end do
        end do
    end do
    !$omp end do

    !$omp do schedule(static)
    do d=1, nso
        do c=1, nso
            do b=1, nso
                do a=1, nso
    
                    tau19(a, b, c, d) = tau19(a, b, c, d) + ( &
                        tau10(a, b, d, c)&
                    )
    
                end do
            end do
        end do
    end do
    !$omp end do

    !$omp single
    tau20 = 0.0
    !$omp end single

    !$omp do schedule(static)
    do d=1, nso
        do c=1, nso
            do b=1, nso
                do a=1, nso
                    do e=1, nso
                        tau20(a, b, c, d) = tau20(a, b, c, d) + ( &
                            Y(a, e) * tau19(b, c, d, e)&
                        )
                    end do
                end do
            end do
        end do
    end do
    !$omp end do


    !$omp do schedule(static)
    do d=1, nso
        do c=1, nso
            do b=1, nso
                do a=1, nso
    
                    tau22(a, b, c, d) = tau22(a, b, c, d) - ( &
                        tau20(a, b, c, d)&
                    )
    
                end do
            end do
        end do
    end do
    !$omp end do


    !$omp do schedule(static)
    do d=1, nso
        do c=1, nso
            do b=1, nso
                do a=1, nso
    
                    t2(a, b, c, d) = t2(a, b, c, d) + ( &
                        tau22(a, c, d, b)&
                    )
    
                end do
            end do
        end do
    end do
    !$omp end do

    !$omp do schedule(static)
    do d=1, nso
        do c=1, nso
            do b=1, nso
                do a=1, nso
    
                    t2(a, b, c, d) = t2(a, b, c, d) - ( &
                        tau22(a, d, c, b)&
                    )
    
                end do
            end do
        end do
    end do
    !$omp end do

    !$omp do schedule(static)
    do d=1, nso
        do c=1, nso
            do b=1, nso
                do a=1, nso
    
                    t2(a, b, c, d) = t2(a, b, c, d) - ( &
                        tau22(b, c, d, a)&
                    )
    
                end do
            end do
        end do
    end do
    !$omp end do

    !$omp do schedule(static)
    do d=1, nso
        do c=1, nso
            do b=1, nso
                do a=1, nso
    
                    t2(a, b, c, d) = t2(a, b, c, d) + ( &
                        tau22(b, d, c, a)&
                    )
    
                end do
            end do
        end do
    end do
    !$omp end do


    !$omp do schedule(static)
    do d=1, nso
        do c=1, nso
            do b=1, nso
                do a=1, nso
    
                    t2(a, b, c, d) = t2(a, b, c, d) + ( &
                        tau10(c, d, a, b)&
                    )
    
                end do
            end do
        end do
    end do
    !$omp end do

    !$omp do schedule(static)
    do d=1, nso
        do c=1, nso
            do b=1, nso
                do a=1, nso
    
                    t2(a, b, c, d) = t2(a, b, c, d) - ( &
                        tau10(c, d, b, a)&
                    )
    
                end do
            end do
        end do
    end do
    !$omp end do

    !$omp do schedule(static)
    do d=1, nso
        do c=1, nso
            do b=1, nso
                do a=1, nso
    
                    t2(a, b, c, d) = t2(a, b, c, d) - ( &
                        tau10(d, c, a, b)&
                    )
    
                end do
            end do
        end do
    end do
    !$omp end do

    !$omp do schedule(static)
    do d=1, nso
        do c=1, nso
            do b=1, nso
                do a=1, nso
    
                    t2(a, b, c, d) = t2(a, b, c, d) + ( &
                        tau10(d, c, b, a)&
                    )
    
                end do
            end do
        end do
    end do
    !$omp end do


    !$omp single
    tau17 = 0.0
    !$omp end single

    !$omp do schedule(static)
    do b=1, nso
        do a=1, nso
            do c=1, nso
                tau17(a, b) = tau17(a, b) + ( &
                    Y(c, a) * c1(b, c)&
                )
            end do
        end do
    end do
    !$omp end do

    !$omp single
    tau18 = 0.0
    !$omp end single

    !$omp do schedule(static)
    do b=1, nso
        do a=1, nso
            do c=1, nso
                tau18(a, b) = tau18(a, b) + ( &
                    Y(a, c) * tau17(b, c)&
                )
            end do
        end do
    end do
    !$omp end do

    !$omp do schedule(static)
    do d=1, nso
        do c=1, nso
            do b=1, nso
                do a=1, nso
    
                    t2(a, b, c, d) = t2(a, b, c, d) - ( &
                        Y(a, c) * tau18(b, d)&
                    )
    
                end do
            end do
        end do
    end do
    !$omp end do

    !$omp do schedule(static)
    do d=1, nso
        do c=1, nso
            do b=1, nso
                do a=1, nso
    
                    t2(a, b, c, d) = t2(a, b, c, d) + ( &
                        Y(b, c) * tau18(a, d)&
                    )
    
                end do
            end do
        end do
    end do
    !$omp end do


    !$omp do schedule(static)
    do d=1, nso
        do c=1, nso
            do b=1, nso
                do a=1, nso
    
                    t2(a, b, c, d) = t2(a, b, c, d) + ( &
                        Y(a, c) * tau17(d, b)&
                    )
    
                end do
            end do
        end do
    end do
    !$omp end do

    !$omp do schedule(static)
    do d=1, nso
        do c=1, nso
            do b=1, nso
                do a=1, nso
    
                    t2(a, b, c, d) = t2(a, b, c, d) - ( &
                        Y(b, c) * tau17(d, a)&
                    )
    
                end do
            end do
        end do
    end do
    !$omp end do


    !$omp single
    tau23 = 0.0
    !$omp end single

    !$omp do schedule(static) reduction(+:tau23)
    
    do b=1, nso
        do a=1, nso
            tau23 = tau23 + ( &
                Y(b, a) * c1(a, b)&
            )
        end do
    end do
    
    !$omp end do

    !$omp single
    tau27 = 0.0
    !$omp end single

    !$omp do schedule(static)
    do b=1, nso
        do a=1, nso
    
            tau27(a, b) = tau27(a, b) + ( &
                tau23 * Y(b, a)&
            )
    
        end do
    end do
    !$omp end do

    !$omp single
    tau24 = 0.0
    !$omp end single

    !$omp do schedule(static)
    do b=1, nso
        do a=1, nso
            do c=1, nso
                tau24(a, b) = tau24(a, b) + ( &
                    Y(a, c) * c1(c, b)&
                )
            end do
        end do
    end do
    !$omp end do

    !$omp single
    tau25 = 0.0
    !$omp end single

    !$omp do schedule(static)
    do b=1, nso
        do a=1, nso
    
            tau25(a, b) = tau25(a, b) + ( &
                tau24(a, b)&
            )
    
        end do
    end do
    !$omp end do


    !$omp do schedule(static)
    do b=1, nso
        do a=1, nso
    
            tau25(a, b) = tau25(a, b) - ( &
                c1(a, b)&
            )
    
        end do
    end do
    !$omp end do

    !$omp single
    tau26 = 0.0
    !$omp end single

    !$omp do schedule(static)
    do b=1, nso
        do a=1, nso
            do c=1, nso
                tau26(a, b) = tau26(a, b) + ( &
                    Y(c, a) * tau25(b, c)&
                )
            end do
        end do
    end do
    !$omp end do


    !$omp do schedule(static)
    do b=1, nso
        do a=1, nso
    
            tau27(a, b) = tau27(a, b) - ( &
                tau26(a, b)&
            )
    
        end do
    end do
    !$omp end do


    !$omp do schedule(static)
    do b=1, nso
        do a=1, nso
    
            tau27(a, b) = tau27(a, b) + ( &
                ! Not supported in Fortran:
    ! IndexedBase
    c0 * Y(b, a)&
            )
    
        end do
    end do
    !$omp end do

    !$omp do schedule(static)
    do d=1, nso
        do c=1, nso
            do b=1, nso
                do a=1, nso
    
                    t2(a, b, c, d) = t2(a, b, c, d) + ( &
                        Y(b, d) * tau27(c, a)&
                    )
    
                end do
            end do
        end do
    end do
    !$omp end do

    !$omp do schedule(static)
    do d=1, nso
        do c=1, nso
            do b=1, nso
                do a=1, nso
    
                    t2(a, b, c, d) = t2(a, b, c, d) - ( &
                        Y(a, d) * tau27(c, b)&
                    )
    
                end do
            end do
        end do
    end do
    !$omp end do


    !$omp single
    
    
    t0 = t0 + ( &
        ! Not supported in Fortran:
    ! IndexedBase
    c0&
    )
    
    
    !$omp end single

    !$omp do schedule(static)
    do b=1, nso
        do a=1, nso
    
            t1(a, b) = t1(a, b) + ( &
                ! Not supported in Fortran:
    ! IndexedBase
    c0 * Y(a, b)&
            )
    
        end do
    end do
    !$omp end do

    !$omp end parallel

	deallocate(tau1)
    deallocate(tau0)
    deallocate(tau5)
    deallocate(tau14)
    deallocate(tau15)
    deallocate(tau13)
    deallocate(tau21)
    deallocate(tau2)
    deallocate(tau4)
    deallocate(tau6)
    deallocate(tau7)
    deallocate(tau9)
    deallocate(tau11)
    deallocate(tau12)
    deallocate(tau16)
    deallocate(tau19)
    deallocate(tau20)
    deallocate(tau22)
    deallocate(tau10)
    deallocate(tau18)
    deallocate(tau17)
    deallocate(tau24)
    deallocate(tau25)
    deallocate(tau26)
    deallocate(tau27)

	End Subroutine BuildSC_Vec

	End Module BuildSC_dr
