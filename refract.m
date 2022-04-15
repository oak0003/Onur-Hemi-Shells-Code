function [s2,theta_1,theta_2] = refract(n1,n2,k2,s1,n)


s1 = s1./(norm(s1));
n = n./(norm(n));    

% Calculate incidence angle

theta_1 = acosd(dot(s1,-n));

% Find intermediate variable P

term1 = ((n2^2) - (k2^2)-(n1^2) * (sind(theta_1))^2)^2;
term2 = 4 * (n2^2) * (k2^2);
term3 = (n2^2) - (k2^2)-(n1^2) * (sind(theta_1))^2;%sqrt(term1);

psquared = 0.5 * ((sqrt(term1 + term2)) + term3);

p = sqrt(psquared);

% Compute refraction angle via generalized Snell's Law
theta_2 = atand((n1*sind(theta_1))/p);

A = sind(theta_2)/sind(180-theta_1); B = sind(theta_1-theta_2)/sind(180-theta_1);

% Compute direction vector of refracted angle
s2 = A.*s1 + B.*(-n);

s2 = s2./(norm(s2));







end