function [Coefs] = latdirStabDerivs(A,B,m, u0, theta0, Ixx, Izz, Izx, rho, S, b)
% Calculate augmented moments of inertia
Ixp = (Ixx*Izz-Izx^2)/Izz;
Izp = (Ixx*Izz-Izx^2)/Ixx;
Izxp = Izx/(Ixx*Izz-Izx^2);

% Turn A and B matrices into Dimensional Derivatives
% Y
    Yv = A(1,1)*m*u0;
    Yp = A(1,2)*m*u0;
    Yr = (A(1,3)-1)*m*u0;
    Yda = B(1,1)*m*u0;
    Ydr = B(1,2)*m*u0;
% Solve for L and N variables using systems of equations
% Lv Nv
    a(1,1) = 1/Ixp;
    a(1,2) = Izxp;
    a(2,1) = Izxp;
    a(2,2) = 1/Izp;
    c(1,1) = A(2,1);
    c(2,1) = A(3,1);
    X = a\c;
    Lv = X(1);
    Nv = X(2);
 % Lp Np
    a(1,1) = 1/Ixp;
    a(1,2) = Izxp;
    a(2,1) = Izxp;
    a(2,2) = 1/Izp;
    c(1,1) = A(2,2);
    c(2,1) = A(3,2);
    X = a\c;
    Lp = X(1);
    Np = X(2);
 % Lp Np
    a(1,1) = 1/Ixp;
    a(1,2) = Izxp;
    a(2,1) = Izxp;
    a(2,2) = 1/Izp;
    c(1,1) = A(2,3);
    c(2,1) = A(3,3);
    X = a\c;
    Lr = X(1);
    Nr = X(2);
  % Lda Nda
    a(1,1) = 1/Ixp;
    a(1,2) = Izxp;
    a(2,1) = Izxp;
    a(2,2) = 1/Izp;
    c(1,1) = B(2,1);
    c(2,1) = B(3,1);
    X = a\c;
    Lda = X(1);
    Nda = X(2);
  % Lda Nda
    a(1,1) = 1/Ixp;
    a(1,2) = Izxp;
    a(2,1) = Izxp;
    a(2,2) = 1/Izp;
    c(1,1) = B(2,2);
    c(2,1) = B(3,2);
    X = a\c;
    Ldr = X(1);
    Ndr = X(2);
    
  % Non Dimensionalize
    Cyb = Yv/(0.5*rho*u0^2*S);
    Cyp = Yp/(0.25*rho*u0*b*S);
    Cyr = Yr/(0.25*rho*u0*b*S);
    Cyda = Yda/(0.5*rho*u0^2*S);
    Cydr = Ydr/(0.5*rho*u0^2*S);
    
    Clb = Lv/(0.5*rho*u0^2*b*S);
    Clp = Lp/(0.25*rho*u0*b^2*S);
    Clr = Lr/(0.25*rho*u0*b^2*S);
    Clda = Lda/(0.5*rho*u0^2*b*S);
    Cldr = Ldr/(0.5*rho*u0^2*b*S);
    
    Cnb = Nv/(0.5*rho*u0^2*b*S);
    Cnp = Np/(0.25*rho*u0*b^2*S);
    Cnr = Nr/(0.25*rho*u0*b^2*S);
    Cnda = Nda/(0.5*rho*u0^2*b*S);
    Cndr = Ndr/(0.5*rho*u0^2*b*S);
    
    Coefsarray = [Cyb;Cyp;Cyr;Cyda;Cydr;Clb;Clp;Clr;Clda;Cldr;Cnb;Cnp;Cnr;Cnda;Cndr];
    Coefs = array2table(Coefsarray,'RowNames',{'Cyb','Cyp','Cyr','Cyda','Cydr','Clb','Clp','Clr','Clda','Cldr','Cnb','Cnp','Cnr','Cnda','Cndr'});
end
    