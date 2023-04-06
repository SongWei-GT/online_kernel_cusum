function [ARL] = ARL_approx_skew_correct(b,omega_B,c3B)

    term_1 = sqrt(2*pi)./b;

    term_2 = 0;
    B_idx = 0;
    for B = omega_B
        B_idx = B_idx + 1;
        temp_term_2 = (2*B - 1)./B./(B-1);
        temp_theta_B = theta_B(b,c3B(B_idx));
        temp_term_3 = temp_theta_B.^2./2 + c3B(B_idx).*temp_theta_B.^3./6 - temp_theta_B.*b;
        
        term_2 = term_2 + exp(temp_term_3).*temp_term_2 .* func_nu(temp_theta_B.*sqrt(2.*temp_term_2));
    end

    ARL = term_1./term_2;

end