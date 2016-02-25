function [ y ] = parseNaN( x )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    if isfinite(x)
        y = x;
    else
        y = 0;
    end

end

