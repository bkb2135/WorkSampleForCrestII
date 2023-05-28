module  test
    contains
        function compute_light_curve(time, p, tau) result(numerical_flux)
        !Compute the light curve for a given time array, planet radius, and transit duration
            real, dimension(2157), intent(in) :: time
            real, intent(in) :: p, tau
            real, dimension(2157) :: numerical_flux, z_arr
            z_arr = (time)/tau
            do i=1, 2157
                flux_obscured = integrator(numerator, p, z_arr(i))
                flux_unobscured = integrator(denominator, p, z_arr(i))
                numerical_flux(i) = flux_obscured/flux_unobscured
            end do
        !compute_light_curve = numerical_flux
        end function compute_light_curve

        real function delta(r, p, z)
        !Delta from Mandel & Agol 2003
            real, intent(in) :: p, r, z
            real, parameter :: PI=4.0*ATAN(1.0)

            IF (r >= z + p .OR. r <= z-p)THEN
                delta = 0.0
            ELSE IF (r + z <= p) THEN
                delta = 1.0
            ELSE 
                delta = (1/PI)*ACOS((z*z - p*p + r*r) / (2.0*z*r))
            END IF
        end function delta

        real function limb_darkening(r) result(ld)
        !Limb darkening profile from Mandel & Agol 2003
            real, intent(in) :: r
            real, parameter :: u1=0.3, u2=0.2
            real :: mu
            mu = sqrt(1.0 - r*r)
            ld = 1.0 - u1*(1.0-mu) - u2*(1.0-mu)*(1.0-mu)
            !ld = 1.0
        end function limb_darkening

        real function denominator(r)
        !Denominator from Mandel & Agol 2003
            real, intent(in) :: r
            denominator = limb_darkening(r)*2*r
        end function denominator

        real function numerator(r, p, z) 
        !Numerator from Mandel & Agol 2003
            real, intent(in) :: p, r, z 
            numerator = limb_darkening(r)*(1-delta(r, p, abs(z)))*2*r
        end function numerator

        real function integrator(f, p, z)
        !Integrate a function f from 0 to 1 using Simpson's rule, taking p and z as parameters
            real, intent(in) :: p, z
            real :: b, a, sum
            integer :: i, n
            b = 1.0
            a = 0.0
            n = 1000
            h = (b-a)/n
            sum = f(a, p, z) + f(b, p, z)   
            do i=1, n-1, 2
                x = a + i*h
                sum = sum + 4.0*f(x, p, z)
            end do
            do i=2, n-2, 2
                x = a + i*h
                sum = sum + 2.0*f(x, p, z)
            end do
            integrator = sum*h/3.0
        end function integrator

        real function compute_chi_squared(y, pred, y_err) result(chi_squared)
        !Currently computes the non-reduced chi-squared
        real, dimension(2157) :: y, pred, y_err
        real :: sum = 0
        integer :: i
        do i = 1, 2157
            sum  = sum+ ((y(i)-pred(i))/pred(i))**2
        end do
        chi_squared = sum
        end function compute_chi_squared

        function mcmc(starting_point, nsteps, step_size, time, flux, flux_err) result(chain)
            real, dimension(2), intent(in) :: starting_point
            integer, intent(in) :: nsteps
            real, intent(in):: step_size
            real, dimension(2157), intent(in) :: time, flux, flux_err
            real, dimension(2, nsteps) :: chain
            real, dimension(2, nsteps) :: prob_chain
            real :: log_prob, log_prob_new
            real, dimension(2) :: x_new, x_old
            integer :: i
            call srand(0)
            x_old = starting_point
            log_prob = lnprob_log(time, x_old, flux, flux_err)
            do i=1, nsteps

                x_new(1) = x_old(1) + step_size*rand() - step_size*0.5
                x_new(2) = x_old(2) + step_size*rand() - step_size*0.5
                log_prob_new = lnprob_log(time, x_new, flux, flux_err)
                if (log_prob_new < log_prob) then
                    chain(:,i) = x_new
                    x_old = x_new
                    log_prob = log_prob_new
                else if (log_prob_new - log_prob < rand()) then
                    chain(:,i) = x_new
                    x_old = x_new
                    log_prob = log_prob_new
                else
                    chain(:,i) = x_old

                end if
                prob_chain(:,i) = [x_old]
            end do
            !chain = prob_chain
        end function mcmc

        real function lnprob_log(time, x, y, y_err) result(log_probability)
        !Actually does not calculate the log probability, but the chi-squared, future versions will correct this
            real,dimension(2157), intent(in) :: y, y_err, time
            real, dimension(2157) :: flux_prediction
            real, dimension(2), intent(in) :: x
            real :: p, tau

            p = x(1)
            tau = x(2)
            flux_prediction = compute_light_curve(time, p, tau)
            log_probability = compute_chi_squared(y, flux_prediction, y_err)

        end function lnprob_log

        function gradient_descent(pin, tauin, time, flux, flux_err) result(minimum)
        !Gradient descent algorithm to find the minimum of the chi-squared function
            real, intent(in) :: pin, tauin
            real, dimension(2157), intent(in) :: time, flux, flux_err
            real :: p, tau
            real, dimension(2) :: minimum
            real :: step_size, chi_squared, chi_squared_new, p_new, tau_new
            !real :: p_gradient, tau_gradient
            integer :: i
            p = pin
            tau = tauin
            step_size = 0.0001
            chi_squared = compute_chi_squared(flux, compute_light_curve(time, p, tau), flux_err)
            do i=1, 10000
                p_new = p - step_size*gradient_p(p, tau, time, flux, flux_err)
                tau_new = tau - step_size*gradient_tau(p, tau, time, flux, flux_err)
                chi_squared_new = compute_chi_squared(flux, compute_light_curve(time, p_new, tau_new), flux_err)
                if (chi_squared_new < chi_squared) then
                    p = p_new
                    tau = tau_new
                    chi_squared = chi_squared_new
                else if (chi_squared_new - chi_squared < rand()) then
                    p = p_new
                    tau = tau_new
                    chi_squared = chi_squared_new
                end if
            end do
            minimum = [p, tau]
        end function gradient_descent


        real function gradient_p(p, tau, time, flux, flux_err) result(p_gradient)
        !Computes the gradient of the chi-squared function with respect to p
            real, intent(in) :: p, tau
            real, dimension(2157), intent(in) :: time, flux, flux_err
            real :: chi_squared, chi_squared_new, p_new
            real :: step_size
            step_size = 0.0001
            chi_squared = compute_chi_squared(flux, compute_light_curve(time, p, tau), flux_err)
            p_new = p + step_size
            chi_squared_new = compute_chi_squared(flux, compute_light_curve(time, p_new, tau), flux_err)
            p_gradient = (chi_squared_new - chi_squared)/step_size
        end function gradient_p

        real function gradient_tau(p, tau, time, flux, flux_err) result(tau_gradient)
        !Computes the gradient of the chi-squared function with respect to tau
            real, intent(in) :: p, tau
            real, dimension(2157), intent(in) :: time, flux, flux_err
            real :: chi_squared, chi_squared_new, tau_new
            real :: step_size
            step_size = 0.0001
            chi_squared = compute_chi_squared(flux, compute_light_curve(time, p, tau), flux_err)
            tau_new = tau + step_size
            chi_squared_new = compute_chi_squared(flux, compute_light_curve(time, p, tau_new), flux_err)
            tau_gradient = (chi_squared_new - chi_squared)/step_size
        end function gradient_tau

        character(len=20) function int2str(k)
        !Convert an integer to string
            integer, intent(in) :: k
            write (int2str, *) k
            int2str = adjustl(int2str)
        end function int2str

end module test

program transit
    use test
    implicit none
    integer, parameter :: arr_size = 2157
    real, dimension(arr_size) :: time, flux, flux_err
    integer :: i
    !Read in both text files and save to arrays
    open(1, file='time.txt', status='old')
    open(2, file='flux.txt', status='old')
    open(3, file='flux_err.txt', status='old')
    do i=1, arr_size
        read(1,*) time(i)
        read(2,*) flux(i)
        read(3,*) flux_err(i)
    end do
    close(1)
    close(2)
    close(3)
    !run the mcmc function 
    call srand(1)
    !Set to one to conserve run time, 
    do i = 1, 1
        open(i, file='chain'//trim(adjustl(int2str(i)))//'.txt', status='unknown')
        write(i,*) mcmc([0.02+0.01*(i-1)*(rand()-0.5), 0.12+(i-1)*0.01*(rand()-0.5)], 1000, 0.001, time, flux, flux_err)
        close(i)
    end do
    !Output can be redirected to a text file and then plotted in python
    open(11, file='grad.txt', status='unknown')
    write(11,*) gradient_descent(0.02, 0.12, time, flux, flux_err)
    close(11)
    
end program