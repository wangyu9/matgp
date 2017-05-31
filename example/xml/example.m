%%

%s = loadxml('curve0.xml');
s = xml2struct('curve0.xml');

%%
tree = struct2xml(s);
save(tree,'tmp.xml');