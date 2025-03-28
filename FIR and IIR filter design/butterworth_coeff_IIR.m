function p = butterworth_coeff_IIR(N, theta_c)
    wc = tan(theta_c/2);
    beta = pi/N;

    par = ~mod(N, 2);
    if par
        disp("par")
        i = -N/2:(N/2)-1;
        p = wc.*exp(1j*(pi + (i.*beta) + (beta/2)));
    else
        disp("impar")
        i = (-1*(N-1)/2):((N-1)/2);
        p = wc.*exp(1j*(pi + (i.*beta)));
    end
end