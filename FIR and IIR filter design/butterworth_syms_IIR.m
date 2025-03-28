function hz = butterworth_syms_IIR(p, tp, wp1, wp2)
    syms s z
    digits(4)
    sympref('FloatingPointOutput',true);

    % Cambio el signo si es impar
    par = ~mod(length(p), 2);
    sign = 1;
    if ~par
        sign = -1;
    end

    % Obtengo numerador y denominador
    num = sign*prod(p);
    den = 1;
    for i=p
        den = den*(s-i);
    end

    % Obtengo H(s)
    hs = eval(num/den);

    % Transofrmacion bilineal a H(z)
    s_to_z = (1-z^-1)/(1+z^-1);
    hz = subs(hs, s, s_to_z);
    [n, d] = numden(hz);
    n = expand(n);
    d = expand(d);
    [~, npowers] = coeffs(n);
    [dcoeffs , dpowers] = coeffs(d);

    % Hallo el exponente mas alto entre el numerador y denominador
    nexponents = mapSymType(npowers, 'power', @(Z) children(Z,2));
    if npowers(end) == 1; nexponents(end) = 0; end

    dexponents = mapSymType(dpowers, 'power', @(Z) children(Z,2));
    if dpowers(end) == 1; dexponents(end) = 0; end
    maxexp = max(double(nexponents(1)), double(dexponents(1)));

    % Divido entre z^maxexp (el exponente maximo de z en la fracción)
    n = expand(n/z^maxexp);
    d = expand(d/z^maxexp);

    % Divido entre el coeficiente libre del denominador
    n = expand(n/dcoeffs(1));
    d = expand(d/dcoeffs(1));

    % Obtengo hz
    hz = n/d;

    %% Rechazabanda
    % % Calculo el Z^-1 rechazabanda
    % alpha = (cos((wp2+wp1)/2))/(cos((wp2-wp1)/2));
    % k = tan((wp2-wp1)/2)*tan(tp/2);
    % z_bandstop_num = (z^-2)-(((2*alpha)/(1+k))*z^-1)+((1-k)/(1+k));
    % z_bandstop_den = (((1-k)/(1+k))*z^-2)-((2*alpha/(1+k))*z^-1)+1;
    % z_bandstop = z_bandstop_num/z_bandstop_den;
    % 
    % % Reemplazo en hz
    % hz = subs(hz, z, z_bandstop);
    % 
    % % Trato de hallar coeficientes nuevamente para hz rechazabanda
    % [n, d] = numden(hz);
    % n = expand(n);
    % d = expand(d);
    % [~, npowers] = coeffs(n);
    % [dcoeffs , dpowers] = coeffs(d);
    % 
    % % Hallo el exponente mas alto entre el numerador y denominador
    % nexponents = mapSymType(npowers, 'power', @(Z) children(Z,2));
    % if npowers(end) == 1; nexponents(end) = 0; end
    % 
    % dexponents = mapSymType(dpowers, 'power', @(Z) children(Z,2));
    % if dpowers(end) == 1; dexponents(end) = 0; end
    % maxexp = max(double(nexponents(1)), double(dexponents(1)));
    % 
    % % Divido entre z^maxexp (el exponente maximo de z en la fracción)
    % n = expand(n/z^maxexp);
    % d = expand(d/z^maxexp);
    % 
    % % Divido entre el coeficiente libre del denominador
    % n = expand(n/dcoeffs(1));
    % d = expand(d/dcoeffs(1));
    % 
    % % Obtengo hz rechaza banda
    % hz = n/d;
end