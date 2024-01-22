function dydt = parameterizedSystemGrad(~, y, a, b, c, q, p, N)

dydt = zeros(2*N,1);

for i=1:N
    xTerm = 0;
    yTerm = 0;
    % xCoeff = (a^c) * (y(i).^(c-1));
    % yCoeff = (b^c) * (y(i+N).^(c-1));
    for j=1:N
        if j ~= i
            xCoeff = (a^c) * (y(i) - y(j) ).^(c-1);
            yCoeff = (b^c) * (y(i+N) - y(j+N) ).^(c-1);
            A = ( a*( y(i) - y(j) ) ).^c;
            B = ( b*( y(i+N) - y(j+N) ) ).^c;
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
