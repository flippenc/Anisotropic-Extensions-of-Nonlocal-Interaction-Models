function dydt = paramGeneralLcGrad(~, y, a, b, c, q, p, N)

dydt = zeros(2*N,1);

for i=1:N
    xTerm = 0;
    yTerm = 0;
    % xCoeff = (a^c) * (y(i).^(c-1));
    % yCoeff = (b^c) * (y(i+N).^(c-1));
    for j=1:N
        if j ~= i
            i1nonABS = y(i)   - y(j)  ;
            i2nonABS = y(i+N) - y(j+N);
            i1ABS    = abs(y(i)   - y(j));
            i2ABS    = abs(y(i+N) - y(j+N));
            
            xCoeff = ((a*i1ABS).^c)/i1nonABS;
            yCoeff = ((b*i2ABS).^c)/i2nonABS;
            
            A = (a*(i1ABS)).^c;
            B = (b*(i2ABS)).^c;
            X = ( A + B );
            Q = (q/c - 1);
            P = (p/c - 1);
            K = (X.^Q - X.^P);
            xTerm = xTerm + (xCoeff .* K);
            yTerm = yTerm + (yCoeff .* K);
        end
    end
    dydt(i) = (-1/N) .* xTerm;
    dydt(i+N) = (-1/N) .* yTerm;
end
