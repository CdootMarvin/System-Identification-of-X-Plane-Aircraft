function [Coefs] = lonStabDerivs(A, B, m, u0, Iyy, rho, S, CL0, theta0, cbar)
% Turn the A and B matrices to Dimensional Derivatives
% X
    Xu = -A(1,1)*m;
    Xalp = -A(1,2)*m;
    Xq = -A(1,3)*m;
    Xde = -B(1,1)*m;
% Z
    Zu = -A(2,1)*m*u0;
    Zalp = -A(2,2)*m*u0;
    Zq = -(A(2,3)-1)*m*u0;
    Zde = -B(2,1)*m*u0;
    
 % M
    Mu = A(3,1)*Iyy;
    Malp = A(3,2)*Iyy;
    Mq = A(3,3)*Iyy;
    Mde = B(3,1)*Iyy;
 % Non-Dimensionalize
 % D
    Cdu = (Xu-rho*u0*S*CL0*sin(theta0))/(0.5*rho*u0*S);
    Cdalp = Xalp/(0.5*rho*u0^2*S);
    Cdq = Xq/(0.5*rho*u0*cbar*S);
    Cdde = Xde/(0.5*rho*u0^2*S);
 % L
    CLu = (Zu+rho*u0*S*CL0*cos(theta0))/(0.5*rho*u0*S);
    CLalp = Zalp/(0.5*rho*u0^2*S);
    CLq = Zq/(0.25*rho*u0*cbar*S);
    CLde = Zde/(0.5*rho*u0^2*S);
 % M
    Cmu = Mu/(0.5*rho*u0*cbar*S);
    Cmalp = Malp/(0.5*rho*u0^2*cbar*S);
    Cmq = Mq/(0.25*rho*u0*cbar^2*S);
    Cmde = Mde/(0.5*rho*u0^2*cbar*S);
    
    coefsarray = [Cdu;Cdalp;Cdq;CLu;CLalp;CLq;Cmu;Cmalp;Cmq;Cdde;CLde;Cmde];
    Coefs = array2table(coefsarray,'RowNames',{'CDu','CDalpha','CDq','CLu','CLalpha','CLq','CMu','CMalpha','CMq','CDde','CLde','CMde'});
end
    
    