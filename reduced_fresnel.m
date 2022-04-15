function rho = reduced_fresnel(theta1, theta2)
    if theta1 == 0
        rho = 0;
    else
    rho = 0.5 * ( ((tand(theta1 - theta2))^2)/((tand(theta1 + theta2))^2)...
    + ((sind(theta1 - theta2))^2)/((sind(theta1 + theta2))^2));
    end
   
end