! The present program simulates a simple outstar network. For simplicity we have
! set the propagation delay to zero, resulting a system of differential (as
! opposed to functional-differential) equations.

program outstar
    implicit none

    ! Constants
    integer, parameter :: n = 2        ! number of input units
    integer, parameter :: nsteps = 100 ! number of steps to take
    real, parameter    :: dt = 0.1     ! (fixed) step size

    ! Model parameters
    real               :: A(0:n)       ! STM decay rates
    real               :: B(n)         ! LTM decay rate
    real               :: I(0:n)       ! Input activities
    real               :: u(0:2*n)     ! STM traces x0 and xi (in u(0:n)) and
                                       ! LTM traces z0i (u(n+1:))

    ! Simulation state
    real               :: t            ! time step

    ! Initialize model parameters
    call random_number(u)
    call random_number(I)
    A = 1.
    B = 1. ! More generally B can be a function of u(0) (e.g., the identity).
    t = 0.

    ! Learning
    call simulation

    ! We expect the normalized STM- and LTM traces to become equal to the
    ! normalized activities in the limit (Grossberg's outstar theorem).
    print *, 'Original input activities: ', I(1:) / sum(I(1:))

contains

    subroutine simulation
        integer :: i ! loop variable

        ! Initialization
        i = 0

        ! Main loop
        do
            ! Print diagnostics (relative STM- and LTM traces)
            print *, t, u(0), u(1:n) / sum(u(1:n)), u(n+1:) / sum(u(n+1:))

            ! Check no. of iterations
            if (i == nsteps) then
                exit
            end if

            ! Take next Runge-Kutta step
            call rk

            ! Update state
            i = i + 1
            t = t + dt
        end do
    end subroutine simulation

    ! Single Runge-Kutta step
    subroutine rk
        real :: k(0:2*n, 4) ! slopes

        ! Beginning of the interval. I is assumed to already contain the inputs
        ! at this time.
        call dudt(k(:, 1), u)

        ! Interval midpoint. I is recalculated.
        call update_inputs
        call dudt(k(:, 2), u + 0.5 * dt * k(:, 1))
        call dudt(k(:, 3), u + 0.5 * dt * k(:, 2))

        ! End of the interval. I is again recalculated, to be reused in the next
        ! RK step for the beginning of the next interval.
        call update_inputs
        call dudt(k(:, 4), u + dt * k(:, 3))

        ! Update the activities.
        u = u + (dt / 6) * (k(:, 1) + 2 * k(:, 2) + 2 * k(:, 3) + k(:, 4))
    end subroutine rk

    ! Evaluates the STM- and LTM equations at the current time step
    subroutine dudt(k, u)
        real, intent(out) :: k(0:2*n) ! slopes
        real, intent(in)  :: u(0:2*n) ! STM- and LTM traces
        real              :: S1       ! signal used in STM update
        real              :: S2       ! signal used in LTM update

        ! Calculate signals for the STM- and LTM equations. We have here picked
        ! two different (sigmoid) functions to illustrate that the signals need
        ! not be the same.
        S1 = 1 / (1 + exp(-u(0)))
        S2 = tanh(u(0))

        ! Calculate dx0dt
        k(0) = -A(0) * u(0) + I(0)

        ! Calculate dxidt, 1 <= i <= n
        k(1:n) = -A(1:n) * u(1:n) + S1 * u(n+1:) + I(1:n)

        ! Calculate dzdt (LTM traces)
        k(n+1:) = -B * u(n+1:) + S2 * u(1:n)
    end subroutine dudt

    ! (Re)calculates the input activities for a given time step
    subroutine update_inputs
        ! Increase the inputs at a fixed rate, maintaining their relative
        ! intensities.
        I = 1.01 * I
    end subroutine update_inputs

end program outstar
