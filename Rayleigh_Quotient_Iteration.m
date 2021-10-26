v = [1;1;1;1;1;1]/6;
A = [1 1 1 1 1 1; 1 2 3 4 5 6; 1 3 6 10 15 21; 1 4 10 20 35 56; 1 5 15 35 70 126; 1 6 21 56 126 252];
s = 1000;
e = 10.^(-10);
sz = size(A);
eigA = eig(A);

l = RQI(v,A,sz, e);

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