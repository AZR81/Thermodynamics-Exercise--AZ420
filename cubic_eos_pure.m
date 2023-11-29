%  Solving Cubic Equations of State
%  Javascript program written by Patrick Barrie
%  with online "help" on how the calculation is performed
%  Program is Copyright Patrick J. Barrie 30 October 2003
%  Please acknowledge me if you use it!

% Peng-Robinson is 4

function cubic_eos_pure(component, temperature, pressure)
    % main routine to do the calculation
    global form
    form.eos = 4;
    form.Tc = 0;  % K
    form.Pc = 0;  % bar
    form.Wc = 0;
    
    form.T = 0;  % K
    form.P = 0;  % bar
    
    form.Z0 = 0;
    form.Z1 = 0;
    form.Z2 = 0;
    form.state0 = 0;
    form.state1 = 0;
    form.state2 = 0;
    form.vol0 = 0;
    form.vol2 = 0;
    form.phi0 = 0;
    form.phi2 = 0;
    form.Hdep0 = 0;
    form.Hdep2 = 0;
    form.Sdep0 = 0;
    form.Sdep2 = 0;
    
    global solver_data
    if component == 1
        % THF
        form.Tc = 540.1;  % K
        form.Pc = 51.9;  % bar
        form.Wc = 0.217;
    else
        % Water
        form.Tc = 647.3;  % K
        form.Pc = 221.2;  % bar
        form.Wc = 0.344;
    end
    form.T = temperature;  % K
    form.P = pressure;  % bar

    calculation();
    
    solver_data.phi = form.phi0;
    
    function rb = calculation()
    
	    R = 8.3145;
	    Tc = (form.Tc);
	    Pc = (form.Pc);
	    Wc = (form.Wc);
	    T = (form.T);	
	    P = form.P;
    
	    paramsab = zeros(4, 1);
	    coeffs = zeros(6, 1);
	    Z = zeros(5, 1);
	    depfns = zeros(4, 1);
	    depfns2 = zeros(4, 1);
	    state = zeros(4, 1);
    
    % check entries on form are sensible:
    
	    temp1 = checkcritparams(Tc, Pc, Wc);
	    temp2 = checkparams(T, P);
	    if ((~(temp1))||(~(temp2)))
	        rb = false;
	    end
    
	    Tr = T/Tc;
	    
    % calculate the coefficients of the cubic equation 
    
	    paramsab = calcab(R, Tc, Pc, Wc, Tr, T, P);
	    coeffs = calccoeffs(paramsab(2),paramsab(3), R, T, P);
    
    % solve the cubic equation for Z
    
	    Z = solvecubic(coeffs(1),coeffs(2),coeffs(3));
    
    % determine fugacity coefficient of 1st root and departure functions 
    
	    depfns = calcdepfns(coeffs(4),coeffs(5),paramsab(1),Z(1), R, T, Tr);
    
    % determine fugacity coefficient of 2nd root (if applicable)
    
	    if (Z(3)>0)
	    
	        depfns2 = calcdepfns(coeffs(4),coeffs(5),paramsab(1),Z(3), R, T, Tr);
	    else
	    
	        depfns2(1)=""; depfns2(2)=""; depfns2(3)="";
	    end
    
    % now determine the state(s)
    
	    state = getstate(Z(4),Z(3),depfns(1),depfns2(1));
    
    % now get molar volumes and round answers to 4 or 5 signficant figures
    
        vol0=Z(1)*R*T/(P*1E5);
        vol1="";
        vol2="";
    
        if (Z(3)>0)
            vol2=Z(3)*R*T/(P*1E5);
        end
    
    % now return values to form
    
	    form.Z0 = Z(1);
	    form.Z1 = Z(2);
	    form.Z2 = Z(3);
	    form.state0 = state(1);
	    form.state1 = state(2);
	    form.state2 = state(3);
	    form.vol0 = vol0;
	    form.vol2 = vol2;
	    form.phi0 = depfns(1);
	    form.phi2 = depfns2(1);
	    form.Hdep0 = depfns(2);
	    form.Hdep2 = depfns2(2);
	    form.Sdep0 = depfns(3);
	    form.Sdep2 = depfns2(3);
    
	    rb = true;
    end
    
    % end of main calculation script
    
    % below are the functions used in the main calculation script
    
    % get the critical parameters of the substance
    
    function filling()
    
	    params = zeros(3, 1);
    
        substance=form.substance;
    
        if (substance~=1)
            params = lookupcritical(substance);
            form.Tc=params(1);
            form.Pc=params(2);
            form.Wc=params(3);
        else
            form.Tc="";
            form.Pc="";
            form.Wc="";
        end
    end
    
    % look up Tc, Pc and Wc if a particular substance is selected
    
    function ra = lookupcritical(substance)
    
        t = zeros(20, 1);
        p = zeros(20, 1);
        w = zeros(20, 1);
	    t(2)=0; p(2)=0; w(2)=0;
	    t(3)=126.1; p(3)=33.94; w(3)=0.040;
	    t(4)=154.6; p(4)=50.43; w(4)=0.022;
	    t(5)=33.3; p(5)=12.97; w(5)=-0.215;
	    t(6)=304.2; p(6)=73.82; w(6)=0.228;
	    t(7)=190.6; p(7)=46.04; w(7)=0.011;
	    t(8)=305.4; p(8)=48.80; w(8)=0.099;
	    t(9)=369.8; p(9)=42.49; w(9)=0.152;
	    t(10)=425.2; p(10)=37.97; w(10)=0.193;
	    t(11)=469.7; p(11)=33.69; w(11)=0.249;
	    t(12)=507.4; p(12)=30.12; w(12)=0.305;
	    t(13)=282.4; p(13)=50.32; w(13)=0.085;
	    t(14)=364.8; p(14)=46.13; w(14)=0.142;
	    t(15)=308.3; p(15)=61.39; w(15)=0.187;
	    t(16)=562.2; p(16)=48.98; w(16)=0.211;
	    t(17)=591.8; p(17)=41.09; w(17)=0.264;
	    t(18)=647.3; p(18)=221.2; w(18)=0.344;
    
	    ra = [t(substance),p(substance),w(substance)];
    end
    
    % check that the critical parameters entered are sensible
    
    function ba = checkcritparams(Tc, Pc, Wc)
        if ((Tc==0)||(Pc==0))
            disp('Please enter critical parameters');
            ba = false;
            return
        end
        if (Tc<0)
            disp('Critical temperature cannot be negative!');
            ba = false;
            return
        end
        if (Pc<0)
            disp('Critical pressure cannot be negative!');
            ba = false;
            return
        end
        
        if ((isnan(Tc))||(isnan(Pc))||(isnan(Wc)))
            disp('Critical parameters must be numbers');
            ba = false;
            return
        end
        
        ba = true;
    end
    
    % check that T and P are entered are sensible
    
    function ba = checkparams(T, P)
        if ((T==0)||(P==0))
            disp('Please enter T and P of interest');
            ba = false;
            return
        end
        if (T<0)
            disp('Absolute temperature cannot be negative!');
            ba = false;
            return
        end
        if (P<0)
            disp('Pressure cannot be negative!');
            ba = false;
            return
        end
    
        if (isnan(T))
            disp('Temperature must be a number');
            ba = false;
            return
        end
        if (isnan(P))
            disp('Pressure must be a number');
            ba = false;
            return
        end
    
        ba = true;
    end
    
    % calculate a and b parameters (depending on form.eos)
    
    function ra = calcab(R, Tc, Pc, Wc, Tr, T, P)
        alpha = 0;
        kappa = 0;
        a = 0;
        b = 0;
        switch (form.eos)
            case 1
                a = 27*R*R*Tc*Tc/(64*Pc*1E5);
                b = R*Tc/(8*Pc*1E5);
	            kappa = 0;
        
            case 2
                a = 0.42748*R*R*power(Tc,2.5)/(Pc*1E5);
                b = 0.08664*R*Tc/(Pc*1E5);
                kappa = 0;
            case 3
                kappa = 0.480 + 1.574*Wc - 0.176*Wc*Wc;
                alpha = power(1+kappa*(1-sqrt(T/Tc)),2);
                a = 0.42748*R*R*Tc*Tc*alpha/(Pc*1E5);
                b = 0.08664*R*Tc/(Pc*1E5);
        
            case 4
                kappa = 0.37464 + 1.54226*Wc - 0.26992*Wc*Wc;
                alpha = power(1+kappa*(1-sqrt(Tr)),2);
                a = 0.45724*R*R*Tc*Tc*alpha/(Pc*1E5);
                b = 0.07780*R*Tc/(Pc*1E5);
        
            case 5
	            kappa = 0.134 + 0.508*Wc - 0.0467*Wc*Wc;
	            alpha = exp((2.00+0.836*Tr)*(1-power(Tr,kappa)));
                a = 0.45724*R*R*Tc*Tc*alpha/(Pc*1E5);
                b = 0.07780*R*Tc/(Pc*1E5);
        
            otherwise
                kappa = 0; a = 0; b = 0;
        end
    
	    ra = [kappa,a,b];
    end
    
    % calculate coefficients in the cubic equation of state
    
    function ra = calccoeffs(a,b, R, T, P)
        A = 0;
        B = 0;
        C0 = 0;
        C1 = 0;
        C2 = 0;
        switch (form.eos)
            case 1
                A = a*P*1E5/power(R*T,2);
                B = b*P*1E5/(R*T);
                C2 = -1 - B;
                C1 = A;
                C0 = -A*B;
        
            case 2
                A = a*P*1E5/(sqrt(T)*power(R*T,2));
                B = b*P*1E5/(R*T);
                C2 = -1;
                C1 = A - B - B*B;
                C0 = -A*B;
        
            case 3
                A = a*P*1E5/power(R*T,2);
                B = b*P*1E5/(R*T);
                C2 = -1;
                C1 = A - B - B*B;
                C0 = -A*B;
        
            case 4
                A = a*P*1E5/power(R*T,2);
                B = b*P*1E5/(R*T);
                C2 = -1+B;
                C1 = A - 3*B*B - 2*B;
                C0 = -A*B + B*B + B*B*B;
        
            case 5
                A = a*P*1E5/power(R*T,2);
                B = b*P*1E5/(R*T);
                C2 = -1+B;
                C1 = A - 3*B*B - 2*B;
                C0 = -A*B + B*B + B*B*B;
        
            otherwise
                C2 = 0; C1 = 0; C0 = 0;
        end
    
        ra = [C0,C1,C2,A,B];
    end
    
    % function for solving cubic equations
    
    function ra = solvecubic(C0,C1,C2)
        Q1 = C2*C1/6 - C0/2 - C2*C2*C2/27;
        P1 = C2*C2/9 - C1/3;
        D = Q1*Q1 - P1*P1*P1;
        temp = 0;
        temp1 = 0;
        temp2 = 0;
        Z0 = 0;
        Z1 = 0;
        Z2 = 0;
        if (D>=0)
            temp1 = power(abs(Q1+sqrt(D)),0.3333333333);
            temp1 = temp1 * (Q1 + sqrt(D)) / abs(Q1 + sqrt(D));
            temp2 = power(abs(Q1-sqrt(D)),0.3333333333);
            temp2 = temp2 * (Q1 - sqrt(D))/abs(Q1 - sqrt(D));
            Z0 = temp1 + temp2 - C2/3;
            % Z1 = "";
            % Z2 = "";
            Z1 = 0;
            Z2 = 0;
        else
            temp1 = Q1*Q1/(P1*P1*P1);
            temp2 = sqrt(1-temp1)/sqrt(temp1);
            temp2 = temp2 * Q1/abs(Q1);
            Phi = atan(temp2);
            if (Phi<0) 
                Phi=Phi+pi; 
            end
            Z0 = 2*sqrt(P1)*cos(Phi/3) - C2/3;
            Z1 = 2*sqrt(P1)*cos((Phi+2*pi)/3) - C2/3;
            Z2 = 2*sqrt(P1)*cos((Phi+4*pi)/3) - C2/3;
            if (Z0<Z1) 
                temp = Z0; Z0 = Z1; Z1 = temp; 
            end
            if (Z1<Z2) 
                temp = Z1; Z1 = Z2; Z2 = temp; 
            end
            if (Z0<Z1) 
                temp = Z0; Z0 = Z1; Z1 = temp; 
            end
        end
        
        ra = [Z0,Z1,Z2,D];
    end
    
    % calculate fugacity coefficient and departure functions
    
    function ra = calcdepfns(A,B,kappa,Z, R, T, Tr)
        phi = 0;
        Hdep = 0;
        Sdep = 0;
        switch (form.eos)
            case 1
                phi = exp(Z - 1 - log(Z-B) - A/Z);
                Hdep = R*T*(Z - 1 - A/Z);
                Sdep = R*log(Z-B);
        
            case 2
                phi = exp(Z - 1 - log(Z-B) - A*log(1+B/Z)/B);
                Hdep = R*T*(Z - 1 - 1.5*A*log(1+B/Z)/B);
                Sdep = R*(log(Z-B)-0.5*A*log(1+B/Z)/B);
        
            case 3
                phi = exp(Z - 1 - log(Z-B) - A*log(1+B/Z)/B);
                temp1 = sqrt(Tr);
                temp2 = kappa*temp1/(1+kappa*(1-temp1));
                Hdep = R*T*(Z - 1 - A*(1+temp2)*log(1+B/Z)/B);
                Sdep = R*(log(Z-B) - A*temp2*log(1+B/Z)/B);
        
            case 4
                temp = log((Z+2.414213562*B)/(Z-0.414213562*B));
                phi = exp(Z - 1 - log(Z-B) - A*temp/(2.828427125*B));
                temp1 = sqrt(Tr);
                temp2 = -kappa*temp1/(1+kappa*(1-temp1));
                Hdep = R*T*(Z - 1 + A*(temp2-1)*temp/(2.828427125*B));
                Sdep = R*(log(Z-B) + A*temp2*temp/(2.828427125*B));
        
            case 5
                temp = log((Z+2.414213562*B)/(Z-0.414213562*B));
                phi = exp(Z - 1 - log(Z-B) - A*temp/(2.828427125*B));
                temp1 = power(Tr,kappa);
                temp2 = 0.836*Tr*(1-(kappa+1)*temp1) - 2.00*kappa*temp1;
                Hdep = R*T*(Z-1 + A*(temp2-1)*temp/(2.828427125*B));
                Sdep = R*(log(Z-B) + A*temp2*temp/(2.828427125*B));
        
            otherwise
                phi = 0; 
                Hdep = ""; 
                Sdep = "";
        end
        ra = [phi,Hdep,Sdep];
    end
    
    % work out which state(s) are present
    
    function ra = getstate(D,Z2,phi0,phi2)
        state0="";
        state1="";
        state2="";
        if (D>=0)
            state0 = "only real root";
        else
            state1 = "no meaning";
            if ((Z2>0)&&(phi0<phi2))
                state0 = "vapour (stable)";
                state2 = "liquid (metastable)";
            end
            if ((Z2>0)&&(phi0>phi2))
                state0 = "vapour (metastable)";
                state2 = "liquid (stable)";
            end
            if ((Z2>0)&&(phi0==phi2))
                state0 = "vapour (stable)";
                state2 = "liquid (stable)";
            end
            if (Z2<0)
                state2 = "no meaning";
                state0 = "only sensible root";
            end
        end
        ra = [state0,state1,state2];
    end
end

