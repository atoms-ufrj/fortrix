function [lngm] = unifac(T,x,nu,R,G,A)
    r = nu*R;
    q = nu*Q;
    m = 1 - 5*q;
    Omega = nu * diag(Q);
    G = Omega*exp(-A/T);
    epsilon = q + (Omega .* log(G))*ones(Q);
    phi = r/(r'*x);
    theta = q/(q'*x);
    L = G*diag(1 ./ (G'*x));
    lngm = m.*log(phi) + 4*q.*log(theta) - phi*m'*x + m + epsilon - Omega*log(G'*x) - L*Omega'*x;
endfunction
