function dydt = LinfODEsolverSystemGrad(~, y)
global q p N

dydt = zeros(2*N,1);

for i=1:N
    xTerm = 0;
    yTerm = 0;
    for j=1:N
        if j ~= i
            i1nonABS     = y(i)   - y(j)  ;
            i2nonABS     = y(i+N) - y(j+N);
            i1ABS        = abs(i1nonABS);
            i2ABS        = abs(i2nonABS);
            nonABSdiff   = i1ABS - i2ABS;
            i2nonABSdiff = i2ABS - i1ABS;
            ABSdiff      = abs( nonABSdiff );
            sum          = i1ABS + i2ABS;

            xCoeff = 1/2 .* (i1nonABS ./ i1ABS) .* ( ( nonABSdiff   ./ ABSdiff) + 1);
            yCoeff = 1/2 .* (i2nonABS ./ i2ABS) .* ( ( i2nonABSdiff ./ ABSdiff) + 1);
            
            A = ( ABSdiff ./ 2) + ( sum ./ 2);
            Q = q-1;
            P = p-1;
            K = (A.^Q - A.^P);
            
            xTerm = xTerm + (xCoeff .* K);
            yTerm = yTerm + (yCoeff .* K);
        end
    end
    dydt(i) = (-1/N) .* xTerm;
    dydt(i+N) = (-1/N) .* yTerm;
end