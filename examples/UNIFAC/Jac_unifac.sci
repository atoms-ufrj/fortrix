function [J] = Jac_unifac(T,x,nu,R,G,A)
    r = nu*R;
    q = nu*Q;
    m = 1 - 5*q;
    Omega = nu * diag(Q);
    G = Omega*exp(-A/T);
    phi = r/(r'*x);
    theta = q/(q'*x);
    L = G*diag(1 ./ (G'*x));
    J = -(m*phi' + phi*m' + L*Omega' + Omega*L') - 4*q*theta' + phi*m'*x*phi' + L*diag(Omega'*x)*L';
endfunction
