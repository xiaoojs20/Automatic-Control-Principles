function syms_frac = ss2tf_fraction_from_NUM_DEN(NUM,DEN, s)
%ss2tf_fraction_from_NUM_DEN Creates a fraction from the NUM and DEN
%vectors and the symbolic variable s
syms numerator(s)
syms denominator(s)
syms func(s)
numerator(s) = 0;
denominator(s) = 0;
num_power = length(NUM);
den_power = length(DEN);
for i = 1:num_power
    numerator(s) = numerator(s)+NUM(i)*s^(num_power-i);
end
for i = 1:den_power
    denominator(s) = denominator(s)+DEN(i)*s^(den_power-i);
end

syms_frac = numerator/denominator;
end