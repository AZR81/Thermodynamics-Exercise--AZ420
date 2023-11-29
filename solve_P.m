function fobj = solve_P(delta_lambdas)
    global data
    global solver_data

    l12_11 = delta_lambdas(1);
    l21_22 = delta_lambdas(2);
    R = 8.3145;
    T = 273.15 + solver_data.Temp;
    fobj = 0;

    if solver_data.Temp == 25
        x_array = data.x1_25c;
        P_array = data.P_25c;
        V1 = 8.16512E-05;
        V2 = 1.80542E-05;
    else
        x_array = data.x1_50c;
        P_array = data.P_50c;
        V1 = 8.43256E-05;
        V2 = 1.82186E-05;
    end

    for i=1:length(x_array)
        x1 = x_array(i);
        x2 = 1 - x1;
        L12 = (V2/V1) * exp(-l12_11/(R*T));
        L21 = (V1/V2) * exp(-l21_22/(R*T));
        x1_x2L12 = x1 + x2 * L12;
        x2_x1L21 = x2 + x1 * L21;
        bracket = (L12 / x1_x2L12) - (L21 / x2_x1L21);
        y1 = exp(-log(x1_x2L12) + x2 * bracket);
        y2 = exp(-log(x2_x1L21) - x1 * bracket);
        
        P_solved = (y1 * x1 * solve_antoine([T, 1]) + y2 * x2 * solve_antoine([T, 2])) * 100000;
        fobj = fobj + ((P_array(i) - P_solved) ^ 2);
    end
    solver_data.Results.dls = delta_lambdas;
    solver_data.Results.rmse = sqrt(fobj / length(x_array));
end