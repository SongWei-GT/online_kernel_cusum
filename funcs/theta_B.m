function [my_theta_B] = theta_B(b,c3)

    term_1 = -1 + sqrt(1+2.*b.*c3);

    my_theta_B = term_1./c3;

end