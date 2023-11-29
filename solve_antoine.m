function output=solve_antoine(inputs)
    global data
    component = inputs(2);
    reverse = false;
    if length(inputs) == 3
        reverse = true;
        P = inputs(1);
    else
        T = inputs(1);
    end
    if reverse == true  % Find saturation temperature
        if component == 1
            output = (data.antoine_b_c1_2 / (data.antoine_a_c1_2 - log10(P))) - data.antoine_c_c1_2;
            if output <= data.antoine_range_c1_1(2)
                output = (data.antoine_b_c1_1 / (data.antoine_a_c1_1 - log10(P))) - data.antoine_c_c1_1;
            end
        else
            output = (data.antoine_b_c2_2 / (data.antoine_a_c2_2 - log10(P))) - data.antoine_c_c2_2;
            if output <= data.antoine_range_c2_1(2)
                output = (data.antoine_b_c2_1 / (data.antoine_a_c2_1 - log10(P))) - data.antoine_c_c2_1;
            end
        end
    else  % Find saturation pressure
        if component == 1
            if T <= data.antoine_range_c1_1(2)
                output = power(10, data.antoine_a_c1_1 - (data.antoine_b_c1_1 / (data.antoine_c_c1_1 + T)));
            else
                output = power(10, data.antoine_a_c1_2 - (data.antoine_b_c1_2 / (data.antoine_c_c1_2 + T)));
            end
        else
            if T <= data.antoine_range_c2_1(2)
                output = power(10, data.antoine_a_c2_1 - (data.antoine_b_c2_1 / (data.antoine_c_c2_1 + T)));
            else
                output = power(10, data.antoine_a_c2_2 - (data.antoine_b_c2_2 / (data.antoine_c_c2_2 + T)));
            end
        end
    end
end