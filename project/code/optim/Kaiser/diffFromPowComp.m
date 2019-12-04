function val = diffFromPowComp(alpha,M,m)
L = 2*m*M;
if (alpha >= 21) && (alpha <= 50)
    beta = 0.5842*(alpha-21)^0.4+0.07886*(alpha-21);
elseif alpha > 50
    beta = 0.1102*(alpha-8.7);
end
p0 = kaiser(L,beta);
c = p0(1:M/2);
s = flipud(p0(M/2+1:M));
sqsum = c.^2+s.^2;
sqsumdB = mag2db(sqsum);
diff = max(sqsumdB)-sqsumdB;
val = max(diff);