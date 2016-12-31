N=11; % number of vertecied
d=N-1; % dimentionarity
V0=simplex(N);

R=random_rotation(d);

V=R*V;

% simple projection t 2d:
V2=V(1:2,:);

% draw all edges, each vertex with each vertex:
for c1=1:N-1
    for c2=c1+1:N
        plot([V2(1,c1) V2(1,c2)],[V2(2,c1) V2(2,c2)],'k.-');
        hold on;
    end
end