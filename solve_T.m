function fobj = solve_T(T_array)
    global solver_data
    
    l12_11 = solver_data.delta_lambdas(1);
    l21_22 = solver_data.delta_lambdas(2);
    if solver_data.triple_precision == true
        x_array1 = linspace(0, solver_data.lower_step, solver_data.x_count);
        x_array2 = linspace(solver_data.lower_step, 0.1, solver_data.x_count);
        x_array3 = linspace(0.1, 1, solver_data.x_count);
        x_array = [x_array1, x_array2, x_array3];

        iter_limit = 3*solver_data.x_count;
        fobj = zeros(3 * solver_data.x_count, 1);
        v1 = zeros(3 * solver_data.x_count, 1);
    else
        x_array = linspace(0, 1, solver_data.x_count);

        iter_limit = solver_data.x_count;
        fobj = zeros(solver_data.x_count, 1);
        v1 = zeros(solver_data.x_count, 1);
    end
    R = 8.3145;
    P = solver_data.pressure;

    fobj_val = 0;

    V1 = 8.16512E-05;
    V2 = 1.80542E-05;
    if P == 7
        V1 = 9.97656E-05;
        V2 = 1.96891E-05;
    end
    
    if solver_data.use_equation_1 == true
        for i=1:iter_limit
            T = T_array(i);
    
            x1 = x_array(i);
            x2 = 1 - x1;
    
            L12 = (V2/V1) * exp(-l12_11/(R*T));
            L21 = (V1/V2) * exp(-l21_22/(R*T));
    
            x1_x2L12 = x1 + x2 * L12;
            x2_x1L21 = x2 + x1 * L21;
            bracket = (L12 / x1_x2L12) - (L21 / x2_x1L21);
            y1 = exp(-log(x1_x2L12) + x2 * bracket);
            y2 = exp(-log(x2_x1L21) - x1 * bracket);
            
            P_solved = (y1 * x1 * solve_antoine([T, 1]) + y2 * x2 * solve_antoine([T, 2]));
            v1(i) = (y1 * x1 * solve_antoine([T, 1])) / P;
            fobj(i) = (((P - P_solved) * 100000) ^ 2);
            fobj_val = fobj_val + fobj(i);
        end
    else
        for i=1:iter_limit
            T = T_array(i);
    
            x1 = x_array(i);
            x2 = 1 - x1;
    
            L12 = (V2/V1) * exp(-l12_11/(R*T));
            L21 = (V1/V2) * exp(-l21_22/(R*T));
    
            x1_x2L12 = x1 + x2 * L12;
            x2_x1L21 = x2 + x1 * L21;
            bracket = (L12 / x1_x2L12) - (L21 / x2_x1L21);
            y1 = exp(-log(x1_x2L12) + x2 * bracket);
            y2 = exp(-log(x2_x1L21) - x1 * bracket);
            
            cubic_eos_pure(1, T, P);
            phi_1 = solver_data.phi;
            cubic_eos_pure(1, T, solve_antoine([T, 1]));
            phi_1_sat = solver_data.phi;
            eq2_correction_1 = (phi_1_sat/phi_1) * exp((V1 * (P - solve_antoine([T, 1]))) / (R*T));
            
            cubic_eos_pure(2, T, P);
            phi_2 = solver_data.phi;
            cubic_eos_pure(2, T, solve_antoine([T, 2]));
            phi_2_sat = solver_data.phi;
            eq2_correction_2 = (phi_2_sat/phi_2) * exp((V2 * (P - solve_antoine([T, 2]))) / (R*T));
            
            P_solved = (eq2_correction_1 * y1 * x1 * solve_antoine([T, 1]) + eq2_correction_2 * y2 * x2 * solve_antoine([T, 2]));
            v1(i) = (eq2_correction_1 * y1 * x1 * solve_antoine([T, 1])) / P;
            fobj(i) = (((P - P_solved) * 100000) ^ 2);
            fobj_val = fobj_val + fobj(i);
        end
    end

    solver_data.Results.rmse = sqrt(fobj_val / length(x_array));
    solver_data.Results.abs_error = sqrt(fobj);
    solver_data.Results.T = T_array;
    solver_data.Results.y1 = v1;
end