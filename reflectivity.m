function rho = reflectivity(theta1,n1,n2,k2)

    term1 = ((n2^2) - (k2^2)-(n1^2) * (sind(theta1))^2)^2;
    term2 = 4 * (n2^2) * (k2^2);
    term3 = (n2^2) - (k2^2)-(n1^2) * (sind(theta1))^2;%sqrt(term1);

    psquared = 0.5 * ((sqrt(term1 + term2)) + term3);
    qsquared = 0.5 * ((sqrt(term1 + term2)) - term3);

    p = sqrt(psquared);
    q = sqrt(qsquared);

    % Calculate reflectivity

    numterm1 = (n1*cosd(theta1)-p)^2;
    denomterm1 = (n1*cosd(theta1)+p)^2;

    numterm2 = (p - n1*sind(theta1)*tand(theta1))^2;
    denomterm2 = (p + n1*sind(theta1)*tand(theta1))^2;

    leftPortion = 0.5 * ((numterm1 + q^2)/(denomterm1 + q^2));
    rightPortion = 1 + ((numterm2 + q^2)/(denomterm2 + q^2));
    rho = leftPortion * rightPortion;
end