function [t] = MarlaRandStream(s)
%MARLARANDSTREAM Summary of this function goes here
%   Detailed explanation goes here
    if isa(s, 'RandStream')
        t = s;
    else  % assume s is an integer seed
        seed = mod(s, 2^32);
        t = RandStream('mt19937ar', 'Seed', seed, ...
            'NormalTransform', 'Ziggurat');
    end
end

