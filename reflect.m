function [s2,theta_1] = reflect(s1,n)

    s1 = s1./(norm(s1));
    
    n = n./(norm(n));
    
    theta_1 = acosd(dot(s1,-n));
    
    s2 = s1 - 2*dot(s1,n) .*  n;

end