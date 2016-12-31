function [V,F,D,A] = read_splocs(fname)

%% For Debuggeing
%fname = './face_sploc_init.h5'
%%
verbose = false;
%%
if(verbose)
    h5disp(fname);
end
info = h5info(fname);
%% Easy part, read V,F
tris = h5read(fname,'/tris');
default = h5read(fname,'/default');
V = double(default)';
F = double(tris)'+1;
%% Read W
D = [];
for ii=1:size(info.Datasets)
    str = [info.Name,info.Datasets(ii).Name];
    is_comp = regexp(str,'\/comp[0-9][0-9][0-9]$');
    if(length(is_comp)==1&&is_comp==1) % this is equivalent to is_comp==1 as []==1 will be [], if([]) will not execute
        comp = permute(double(h5read(fname,str)),[2 1]);
        comp = comp - V; % note this!! 
        comp = comp(:);
        D = [D,comp];
    end
end
%% Read A
A = [];
for ii=1:size(info.Datasets)
    str = [info.Name,info.Datasets(ii).Name];
    is_actv = regexp(str,'\/actv[0-9][0-9][0-9]$');
    if(length(is_actv)==1&&is_actv==1)
        actv = double(h5read(fname,str));
        A = [A;actv'];
    end
end


%%
%tsurf(F,V(:,1),V(:,2),V(:,3), W(1:size(V,1,:),:));