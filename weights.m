function w = weights(p,m)
    %{
    Function to generate weights, w for a derivative of order m,
    in our case m=1, central differnence stencil for order p.
    found at: https://en.wikipedia.org/wiki/Finite_difference_coefficient
    %}
    p = p/2; % order of central difference stencil
    %m = 1; % order of derivative
    
    % algoritmen är stulen direkt från wikipedialänken ovan
    A = ones(2*p+1);
    b = zeros(2*p+1,1);
    ROW = -p:p; % varje rad består av talen mellan -p till p upphöjt till radnumret med start på 0
    b(m+1) = factorial(m); % man sätter rad (m+1) i b vektorn till m!
    for i=1:2*p+1
        A(i,:) = ROW.^(i-1); % sätt rad nummer (i) till ROW.^(i-1), det här är ekvivalent med loopen nedanför
        % for j=1:2*p+1
        %     A(i,j) = A(i,j)*(-p+j-1)^(i-1);
        % end

    end
    w = (A\b)'; % vi får inte glömma att dividera med h^m efter/medan vi assemblar matrisen
end