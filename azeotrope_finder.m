function fobj = azeotrope_finder(input)
    global solver_data
    
    mode = solver_data.finder_mode;
    T_or_P = solver_data.finder_T_or_P;  % Could be T or P depending on the mode
    R = 8.3145;
    fobj = [0, 0];
    l12_11 = solver_data.finder_lambdas(1);
    l21_22 = solver_data.finder_lambdas(2);
    
    if mode == "T"
        T = input(2);
        P = T_or_P;
        V1 = 8.62074E-05;
        V2 = 1.83779E-05;
    else
        T = T_or_P + 273.15;
        P = input(2);
        if T_or_P == 25
            V1 = 8.16512E-05;
            V2 = 1.80542E-05;
        else
            V1 = 8.43256E-05;
            V2 = 1.82186E-05;
        end
    end

    x1 = input(1);
    x2 = 1 - x1;

    L12 = (V2/V1) * exp(-l12_11/(R*T));
    L21 = (V1/V2) * exp(-l21_22/(R*T));

    x1_x2L12 = x1 + x2 * L12;
    x2_x1L21 = x2 + x1 * L21;
    bracket = (L12 / x1_x2L12) - (L21 / x2_x1L21);
    y1 = exp(-log(x1_x2L12) + x2 * bracket);
    y2 = exp(-log(x2_x1L21) - x1 * bracket);

    if solver_data.use_equation_1 == true
        fobj(1) = (y1 - (P/solve_antoine([T, 1]))) * 100;
        fobj(2) = (y2 - (P/solve_antoine([T, 2]))) * 100;
    else
        cubic_eos_pure(1, T, P);
        phi_1 = solver_data.phi;
        cubic_eos_pure(1, T, solve_antoine([T, 1]));
        phi_1_sat = solver_data.phi;
        eq2_correction_1 = (phi_1_sat/phi_1) * exp((V1 * 100000 * (P - solve_antoine([T, 1]))) / (R*T));
        
        cubic_eos_pure(2, T, P);
        phi_2 = solver_data.phi;
        cubic_eos_pure(2, T, solve_antoine([T, 2]));
        phi_2_sat = solver_data.phi;
        eq2_correction_2 = (phi_2_sat/phi_2) * exp((V2 * 100000 * (P - solve_antoine([T, 2]))) / (R*T));
        solver_data.Results.poynting = [exp((V1 * 100000 * (P - solve_antoine([T, 1]))) / (R*T)), exp((V2 * 100000 * (P - solve_antoine([T, 2]))) / (R*T))];
        solver_data.Results.phi = [phi_1, phi_1_sat; phi_2, phi_2_sat];
        fobj(1) = (y1 - (P/(eq2_correction_1 * solve_antoine([T, 1])))) * 100;
        fobj(2) = (y2 - (P/(eq2_correction_2 * solve_antoine([T, 2])))) * 100;
        solver_data.Results.gamma = [y1, y2];
    end

    if mode == "T"
        solver_data.Results.finder_T_or_P = T;
        solver_data.Results.finder_y1 = x1;
    else
        solver_data.Results.finder_T_or_P = P;
        solver_data.Results.finder_y1 = x1;
    end
end