function fobj = wilson_temp_fitter(delta_lambdas)
    global data
    global solver_data

    l12_11 = delta_lambdas(1);
    l21_22 = delta_lambdas(2);
    R = 8.3145;
    fobj = 0;
    P = 1.013;

    T_array = data.T_1b;
    x_array = data.x1_1b;
    y_array = data.y1_1b;

    V1 = 8.62074E-05;
    V2 = 1.83779E-05;

    for i=1:length(T_array)
        T = T_array(i);

        x1 = x_array(i);
        x2 = 1 - x1;

        L12 = (V2/V1) * exp(-l12_11/(R*T));
        L21 = (V1/V2) * exp(-l21_22/(R*T));

        x1_x2L12 = x1 + x2 * L12;
        x2_x1L21 = x2 + x1 * L21;
        bracket = (L12 / x1_x2L12) - (L21 / x2_x1L21);
        y1 = exp(-log(x1_x2L12) + x2 * bracket);

        y_solved = (y1 * x1 * solve_antoine([T, 1])) / P;
        fobj = fobj + (((y_array(i) - y_solved) * 100) ^ 2);
    end
    solver_data.Results.dls = delta_lambdas;
    solver_data.Results.rmse = sqrt(fobj / (length(x_array)));
end