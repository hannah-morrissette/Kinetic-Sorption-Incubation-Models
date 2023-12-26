%% Saline Analytical Function

function y = fun_analytical_SfixS0(beta,x)
    global y0S
    ka = beta(1);
    kd = beta(2);
    LS0 = y0S(1,1);  %low S0
    HS0 = y0S(1,2);  %high S0
    LW0 = y0S(2,1);  %low w0
    HW0 = y0S(2,2);  %high w0
    z = exp(-(ka + kd)*x);
    y = zeros(30,1);
    for i = 1:15
        y(i) = (1 / (ka + kd))*((LS0 * kd * (1 - z(i))) + (LW0 * (ka + (kd * z(i)))));
    end
    for i = 16:30
        y(i) = (1 / (ka + kd))*((HS0 * kd * (1 - z(i))) + (HW0 * (ka + (kd * z(i)))));
    end
end