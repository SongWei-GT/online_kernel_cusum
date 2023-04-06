function [my_nu] = func_nu(x)

    term_1 = 2./x.*(normcdf(x./2) - 0.5);

    term_2 = x./2.*normcdf(x./2) + normpdf(x./2);

    my_nu = term_1./term_2;

end