function [MMD1] = h_RBF(x1,x2,y1,y2,sgma)
    MMD1 = RBF(x1,x2,sgma) + RBF(y1,y2,sgma) - RBF(x1,y2,sgma) - RBF(x2,y1,sgma);             
end