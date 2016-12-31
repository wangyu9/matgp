function [t] = render_component(V,F,dV,varargin)

% V is the mean shape or first frame
% F is triangle list
% dV is the deformation component

dV = reshape(dV,size(V)); % if dV is organized in a n*dim by 1 vector
CV = sqrt(sum(dV.*dV,2));
% normalization
dV = dV / max(CV);
CV = CV / max(CV);

non_neg = false;
A = 1; % default magnitude
nf = 40;
iteration = 10;

nvar = length(varargin);
ii = 1;
while(ii<=nvar)
   if(strcmp(varargin{ii},'positive'))
      non_neg = true;
   elseif(strcmp(varargin{ii},'magnitude'))
      A = varargin{ii+1};
      ii = ii + 1;
   elseif(strcmp(varargin{ii},'frames'))
      nf = varargin{ii+1};
      ii = ii + 1;
   elseif(strcmp(varargin{ii},'iterations'))
      iteration = varargin{ii+1};
      ii = ii + 1;   
   else
      break; % stop for unknown para and bypass to next level.
   end
   ii = ii + 1;
end

A = A * 0.07 * mean(max(V)-min(V));

%% render input frames.
hold off;
axis manual;

if(ii<=nvar)
    % bypass remaining arguments if there is any
    [t,~,~] = render_mesh(V,F,'view',[90 90],'ScaleColor',CV,varargin(ii:end));
else
    [t,~,~] = render_mesh(V,F,'view',[90 90],'ScaleColor',CV);
end

%%
for iter = 1:iteration
    for i=0:nf
       factor = A * sin(i/nf*2*pi);
       if(non_neg)
            factor = abs(factor);
       end
       set(t,'Vertices',V+dV*factor);
       drawnow
    end
end