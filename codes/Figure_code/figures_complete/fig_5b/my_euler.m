function [ Tspan, f, rhs ] = my_euler( P_LF, S_LF, P_DG, S_DG, Tspan, step, IC )

f = zeros(1, length(Tspan));
rhs = zeros(1, length(Tspan)-1);
f(1) = IC;

for i = 1:(length(Tspan)-1)

rhs(i) = right_hand_side_ratio(Tspan(i), f(i), P_LF, S_LF, P_DG, S_DG);
    
f(i + 1) = f(i) + step * rhs(i);

end