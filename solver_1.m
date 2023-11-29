global solver_data
%% User-specified data
solver_data.l25 = [2293.2, 8098.2];  % Fitted 25 C
solver_data.l50 = [3572.1, 8312.7];  % Fitted 50 C
solver_data.l1013 = [5174.4, 8581.4];  % Extrapolated for 1.013 bar
solver_data.l1013f = [5410.6, 8381.7];  % Fitted 1.013 bar

solver_data.x_count = 50;  % Number of x1 points
solver_data.triple_precision = true;  % Increase the number of points between 0 <= x1 <= 0.1
solver_data.lower_step = 0.0002;  % Lower key value for triple precision, ...
% so that 0 to lower_step, lower_Step to 0.1 and 0.1 to 1 all have the same number of points
solver_data.use_original = false;  % Use the experimental values of x1

%% Output variables
solver_data.Results.dls = 0;  % Fitted Wilson parameters
solver_data.Results.rmse = 0;  % Root-mean square error
solver_data.Results.T = 300;  % Solved temperatures for a given x1 array
solver_data.Results.y1 = 0;  % Solved y1 for a given x1 array
solver_data.Results.finder_T_or_P = 0;  % Azeotrope temperature or pressure
solver_data.Results.finder_y1 = 0;  % Azeotrope composition
solver_data.Results.abs_error = 0;  % Absolute error at each x1
solver_data.Results.poynting = 0;  % Poynting correction terms
solver_data.Results.phi = 0;  % Fugacity coefficients (regular and saturated)
solver_data.Results.gamma = 0;  % Activity coefficients

%% Communication between functions
solver_data.use_equation_1 = true;  % used by solve_T
solver_data.delta_lambdas = 0;  %  Values of the Wilson parameters currently being used
solver_data.Temp = 0;  % Temperature in C
solver_data.pressure = 0;  % bar
solver_data.phi = 0;  % Vapour fugacity coefficient

% azeotrope_finder data
solver_data.finder_mode = 0;  % T or P
solver_data.finder_lambdas = 0;  % Values of Wilson parameters
solver_data.finder_T_or_P = 0;  % Specified temperature or pressure

%% Variables used in this file
dataholder;

display_option = 'iter';
% display_option = 'final-detailed';

% Options for fminunc
options1 = optimoptions('fminunc','Display', display_option,'TolFun',1e-4,'MaxFunEvals',1e20,'MaxIterations',10000);

% Options for fsolve
options2 = optimoptions('fsolve','Display', display_option,'TolFun',1e-12,'MaxFunEvals',1e20,'MaxIterations',10000, 'StepTolerance', 1e-12);
options2.Algorithm = "levenberg-marquardt";

%% Question 1
solver_data.Temp = 25;

% Question a (fit the Wilson intraction parameters)
fminunc(@solve_P, [1, 1], options1);
solver_data.Results.rmse
solver_data.Results.dls

% Question b and c (bubble point and azeotrope plots)
graph_type = "bubble_point";  % set to "bubble_point" for the bubble_point graph, set to "azeotrope" for the azeotrope y-x graph
plot_P(solver_data.Results.dls, graph_type)

% Question c (azeotrope composition and temperature)
solver_data.finder_lambdas = solver_data.l25;
solver_data.finder_mode = 'P';
solver_data.finder_T_or_P = 25;
fsolve(@azeotrope_finder, [0.5, 1], options2);
azeotrope_T = solver_data.Results.finder_T_or_P  % Azeotrope temperature
azeotrope_x = solver_data.Results.finder_y1  % Azeotrope composition

%% Question 2
solver_data.Temp = 50;

% Question a (fit the Wilson intraction parameters)
fminunc(@solve_P, [1, 1], options1);
solver_data.Results.rmse
solver_data.Results.dls

% Question b and c (bubble point and azeotrope plots)
graph_type = "bubble_point";  % set to "bubble_point" for the bubble_point graph, set to "azeotrope" for the azeotrope y-x graph
plot_P(solver_data.Results.dls, graph_type)

% Question c (azeotrope composition and temperature)
solver_data.finder_lambdas = solver_data.l50;
solver_data.finder_mode = 'P';
solver_data.finder_T_or_P = 50;
fsolve(@azeotrope_finder, [0.5, 1], options2);
azeotrope_T = solver_data.Results.finder_T_or_P  % Azeotrope temperature
azeotrope_x = solver_data.Results.finder_y1  % Azeotrope composition

%% Question 3
% Question a
% Extrapolate the parameters
p_25_mean = mean(data.P_25c)/1e5;
p_50_mean = mean(data.P_50c)/1e5;
pressures = [p_25_mean, p_50_mean];
l12 = [solver_data.l25(1), solver_data.l50(1)];
l21 = [solver_data.l25(2), solver_data.l50(2)];

l12_1013bar = interp1(pressures, l12, 1.013, "linear", 'extrap');
l21_1013bar = interp1(pressures, l21, 1.013, "linear", 'extrap');

% Compare the experimental and predicted values. Not used in the report.
plot_y([l12_1013bar, l21_1013bar])

% Questions a and b (plots)
% Find T for a range of x1 values
guess_T = zeros(solver_data.x_count, 1);
for i=1:3 * solver_data.x_count
    guess_T(i) = 300;
end
solver_data.delta_lambdas = [l12_1013bar, l21_1013bar];
solver_data.pressure = 1.013;
fsolve(@solve_T, guess_T, options2);
solver_data.Results.rmse
graph_type = "equilibrium";  % set to "equilibrium" for equilibrium y-x graph (a), set to "phase" for the phase T-x,y diagram (b)
plot_xyT(graph_type)

% Question c (azeotrope composition and pressure)
solver_data.finder_lambdas = solver_data.l1013;
solver_data.finder_mode = 'T';
solver_data.finder_T_or_P = 1.013;
fsolve(@azeotrope_finder, [0.5, 300], options2);
azeotrope_P = solver_data.Results.finder_T_or_P
azeotrope_x = solver_data.Results.finder_y1

%% Question 4
% fit the Wilson intraction parameters to the experimental data at 1.013 bar
fminunc(@wilson_temp_fitter, [1, 1], options1);
solver_data.Results.dls
solver_data.Results.rmse

p_25_mean = mean(data.P_25c)/1e5;
p_50_mean = mean(data.P_50c)/1e5;
pressures = [p_25_mean, p_50_mean, 1.013];
l12 = [solver_data.l25(1), solver_data.l50(1), solver_data.l1013f(1)];
l21 = [solver_data.l25(2), solver_data.l50(2), solver_data.l1013f(2)];

l12_7bar = interp1(pressures, l12, 7, "linear", 'extrap');
l21_7bar = interp1(pressures, l21, 7, "linear", 'extrap');

[l12_7bar, l21_7bar]  % Print the extrapolated Wilson parameters

solver_data.delta_lambdas = [l12_7bar, l21_7bar];  % Use the extrapolated values
solver_data.delta_lambdas = solver_data.l1013f;  % Use the values fitted to the 1.013 bar data

guess_T = zeros(solver_data.x_count, 1);
for i=1:3 * solver_data.x_count
    guess_T(i) = 450;
end
solver_data.pressure = 7;

% Questions a and b
solver_data.use_equation_1 = true;  % set to false to use equation 2 (not shown in report)

fsolve(@solve_T, guess_T, options2);
solver_data.Results.rmse
solver_data.Results.abs_error

graph_type = "equilibrium";  % set to "equilibrium" for equilibrium y-x graph (a), set to "phase" for the phase T-x,y diagram (b)
plot_xyT(graph_type)

% Question c and d
solver_data.finder_lambdas = solver_data.delta_lambdas;
solver_data.finder_mode = 'T';
solver_data.finder_T_or_P = 7;

% Using equation 1 (c)
fsolve(@azeotrope_finder, [0.5, 425], options2);
azeotrope_T_Equation_1 = solver_data.Results.finder_T_or_P
azeotrope_x_Equation_1 = solver_data.Results.finder_y1

% Using equation 2 (d)
solver_data.use_equation_1 = false;
fsolve(@azeotrope_finder, [0.5, 425], options2);
azeotrope_T_Equation_2 = solver_data.Results.finder_T_or_P
azeotrope_x_Equation_2 = solver_data.Results.finder_y1
azeotrope_poynting_values = solver_data.Results.poynting
azeotrope_fugacity_values = solver_data.Results.phi
azeotrope_activity_coefficients = solver_data.Results.gamma



%% Plotters
function plot_P(delta_lambdas, graph_type)
    global data
    global solver_data
    l12_11 = delta_lambdas(1);
    l21_22 = delta_lambdas(2);
    R = 8.3145;
    T = 273.15 + solver_data.Temp;

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
    if solver_data.triple_precision == true && solver_data.use_original == false
        original_x = x_array;
        x_array1 = linspace(0, 0.01, solver_data.x_count);
        x_array2 = linspace(0.01, 0.1, solver_data.x_count);
        x_array3 = linspace(0.1, 1, solver_data.x_count);
        x_array = [x_array1, x_array2, x_array3];
    elseif solver_data.triple_precision == false && solver_data.use_original == false
        original_x = x_array;
        x_array = linspace(0, 1, solver_data.x_count);
    end

    P_solved = zeros(length(x_array), 1);
    v1 = zeros(length(x_array), 1);

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
        P_solved(i) = (y1 * x1 * solve_antoine([T, 1]) + y2 * x2 * solve_antoine([T, 2])) * 100000;
        v1(i) = (y1 * x1 * solve_antoine([T, 1]) * 100000)/P_solved(i);
    end
    figure
    hold on
    
    if graph_type == "bubble_point"
        % Bubble point data (prediction vs experiments)
        plot(x_array, P_solved);
        plot(original_x, P_array, 'x')
        legend("Predicted", "Experimental Data", "Location", "southeast");
        xlabel("x_1")
        ylabel("Pressure (Pa)")
        set(gca,'fontname','times');
        set(gca,'fontsize',12);
        box on
    elseif graph_type == "azeotrope"
        % xy Azeotrope plot
        plot(x_array, v1);
        plot(x_array, x_array);
        xlabel("x_1")
        ylabel("y_1")
        xticks(0:0.1:1)
        set(gca,'fontname','times');
        set(gca,'fontsize',12);
        box on
    end

    hold off
end


function plot_y(delta_lambdas)
    global data
    global solver_data
    l12_11 = delta_lambdas(1);
    l21_22 = delta_lambdas(2);
    R = 8.3145;
    P = 1.013;
    
    T_array = data.T_1b;
    x_array = data.x1_1b;
    y_array = data.y1_1b;

    y_solved = zeros(length(x_array), 1);
    sqerror = 0;
    
    V1 = 8.62074E-05; % 1.013
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
        y_solved(i) = (y1 * x1 * solve_antoine([T, 1])) / P;
        sqerror = sqerror + (((y_array(i) - y_solved(i))*100) ^ 2);
    end
    solver_data.Results.rmse = sqrt(sqerror / length(x_array));
    figure
    hold on

    % Bubble point data (prediction vs experiments)
    plot(x_array, y_solved);
    plot(x_array, y_array, 'x')
    plot(x_array, x_array, "Color", "black");
    legend("Predicted", "Experimental Data", "Location", "southeast");
    xlabel("x_1")
    ylabel("y_1")
    set(gca,'fontname','times');
    set(gca,'fontsize',12);
    box on

    hold off
end


function plot_xyT(graph_type)
    global solver_data
    global data
    
    T_array = solver_data.Results.T;
    x_array1 = linspace(0, solver_data.lower_step, solver_data.x_count);
    x_array2 = linspace(solver_data.lower_step, 0.1, solver_data.x_count);
    x_array3 = linspace(0.1, 1, solver_data.x_count);
    x_array = [x_array1, x_array2, x_array3];
    y_array = solver_data.Results.y1;
    d_array = zeros(3*solver_data.x_count, 1);
    for i=1:2*solver_data.x_count
        d_array(i) = y_array(i) - x_array(i);
    end

    figure
    hold on
    
    if graph_type == "equilibrium"
        % Equilibrium calculations
        plot(x_array, y_array);
        plot(x_array, x_array, "Color", "black");
        if solver_data.pressure ~= 7
            plot(data.x1_1b, data.y1_1b, 'x');
            legend("Predicted", "Experimental Data", "Location", "southeast");
        end
        xlim([0, 1])
        ylim([0, 1])
        xlabel("x_1")
        ylabel("y_1")
        set(gca,'fontname','times');
        set(gca,'fontsize',12);
        xticks(0:0.1:1)
        box on
    elseif graph_type == "phase"
        % Phase diagram
        plot(x_array, T_array);
        plot(y_array, T_array);
        if solver_data.pressure ~= 7
            plot(data.x1_1b, data.T_1b, 'x');
            plot(data.y1_1b, data.T_1b, 'x');
            legend("Predicted Liquid", "Predicted Vapour", "Experimental Liquid", "Experimental Vapour", "Location", "northeast");
        else
            legend("Liquid", "Vapour");
        end
        xlim([0, 1])
        xlabel("x_1, y_1")
        ylabel("Temperature (K)")
        set(gca,'fontname','times');
        set(gca,'fontsize',12);
        xticks(0:0.1:1)
        box on
    end
    
    hold off
end