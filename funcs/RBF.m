function [RBF1] = RBF(x1,x2,sgma)
    RBF1 = exp(-norm(x1-x2,2).^2./2./(sgma.^2));       
end


