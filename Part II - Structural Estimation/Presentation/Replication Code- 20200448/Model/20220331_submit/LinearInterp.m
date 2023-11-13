function [out_value] = LinearInterp(Par, Par_expand, in_value)

% This function impletes linear interpolation. The inputs are 
% Par: parameters
% in_value: the values to tbe interpolated. It can be a value function or a distribution density
% Par_expand: parameters containing a finer set of grids

% value is stored as [nf*nA]*[nE*nL]
    
    value      = reshape(in_value, Par.nA,[]);    
    value_A    = interp1(Par.AGrid, value, Par_expand.AGrid,'linear','extrap');    
    
    value_A    = reshape(value_A', Par.nf,[]);
    value_Af   = interp1(Par.fGrid, value_A, Par_expand.fGrid,'linear','extrap');
    
    value_Af   = reshape(value_Af', Par.nL,[]);
    value_AfL  = interp1(Par.LGrid, value_Af, Par_expand.LGrid,'linear','extrap');    
    
    value_AfL  = reshape(value_AfL', Par.nE,[]);
    value_AfLE = interp1(Par.EGrid, value_AfL, Par_expand.EGrid,'linear','extrap');    
    
    value_AfLE = reshape(value_AfLE',Par_expand.nA*Par_expand.nf,[]);
    out_value  = value_AfLE; 
    
end

