v = [1;1;1;1;1;1]/sqrt(6);
A = [1 1 1 1 1 1; 1 2 3 4 5 6; 1 3 6 10 15 21; 1 4 10 20 35 56; 1 5 15 35 70 126; 1 6 21 56 126 252];
s = 1000;
e = 10.^(-5);
sz = size(A);
eigA = eig(A);

v = v/(norm(v,2));
mu = v'*A*v;

l = RQI(v,A,sz,e);

% RQI is a function that returns the result of Rayleigh Quotient Iteration
% The return is an eigenvalue that is closest to the mu from the initial
% vector that we picked.
function l = RQI(v, A, sz, e)
    v = v/(norm(v,2));
    l = v'*A*v;
    oldl = 0;
    while abs(oldl - l) > e
        newA = A - l*eye(sz);
        w = newA\v;
        newv = w/(norm(w,2));
        oldl = l;
        l = newv'*A*newv;
    end
end