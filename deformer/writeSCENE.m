function writeSCENE(xmlFileName,V,F,UV,texture,isoline,argu)

    docRoot = com.mathworks.xml.XMLUtils.createDocument('scene');

    %docRootNode.setAttribute('attr_name','attr_value');

    %     thisElement = docNode.createElement('TriangleMesh'); 
    %     docNode2 = thisElement.createElement('positions'); 
    %     thisElement2.appendChild... 
    %         (docNode2.createTextNode(sprintf('%i',i)));
    %     thisElement.appendChild(thisElement2);

    write_camera(docRoot);
    write_light(docRoot);
    
    M = material_property('black');
    %M.shadeColor = [0 0 1];
    if ~isempty(isoline.source_point)
        write_points(docRoot,isoline.source_point,M);
    end
    
    M = material_property('default');
    write_textured_triangle_mesh(docRoot,V,F,UV,M,texture);
    
    if(argu.draw_plane)
        fr = min(V(:,2));
        Vf = [   05	fr -05;...
             -05	fr -05;...
             -05	fr  05;...
              05	fr  05];
        Ff = [0 1 2;0 2 3]+1; 
        M = material_property('blank');
        M = material_property('grey');
        write_triangle_mesh(docRoot,Vf,Ff,[],M);
    end
    
    if(isoline.draw)
        M = material_property('black');
        write_isolines(docRoot,isoline.V,isoline.F,isoline.u,20,M);
    end
    
    if(size(argu.quiver,1)>0)
        scale = mean(max(V)-min(V)) * 0.05;
        M = material_property('black');
        for i=1:size(argu.quiver)
            WN = argu.quiver(i,4:6);
            [QV,QF] = quiver_mesh(WN./norm(WN),scale*norm(WN));
            QV = QV + argu.quiver(i,1:3);
            if(size(argu.quiver,2)>6)
               M = material_property_quiver(argu.quiver(i,7:9));
            end
            write_triangle_mesh(docRoot,QV,QF,[],M);
        end
    end
    
    xmlwrite(xmlFileName,docRoot);
    %type(xmlFileName);

end

function [] = write_isolines(docRoot,V,F,u,k,M)

    add_floor = false;
    if add_floor
        iso = min(u):(max(u)-min(u))/k:max(u);
        %iso(1) = min(u)+(max(u)-min(u))/800;
        iso(1) = min(u)+(max(u)-min(u))/8000;
        [~,~,t] = IsoLine_helper({F,V},u,iso);
    else
        %[~,~,t] = IsoLine_helper({F,V},u,k);
        iso = min(u):(max(u)-min(u))/k:max(u);
        iso(end) = [];
        iso(1) = [];
        [~,~,t] = IsoLine_helper({F,V},u,iso);
    end
    
    for i=1:size(t.isolines,1)
       X = t.isolines{i,1};
       
%        X = X(~any(isnan(X),2),:);
%        x = size(X,1);
%        L = 1:x;
%        write_isolines(docRoot,[X,0.002*ones(x,1)],L-1,M);
       
       linewidth = 0.002;
       if add_floor          
           if i==1
               linewidth = 0.005;
           end
       end
       
       x = size(X,1);
       ss = 1;
       while(isnan(X(ss,1)))
           ss = ss + 1;
       end
       ee = ss;
       while(ee<x)
           if(isnan(X(ee,1)))
               assert(~isnan(X(ss,1)));
               VL = X(ss:ee-1,:);
               vl = size(VL,1);
               write_linesegments(docRoot,[VL,linewidth*ones(vl,1)],0:vl-1,M);
               ss = ee + 1;
               while(isnan(X(ss,1)))
                   ss = ss + 1;
               end
               ee = ss;
           end
           ee = ee + 1;
       end
       VL = X(ss:ee-1,:);
       vl = size(VL,1);
       write_linesegments(docRoot,[VL,linewidth*ones(vl,1)],0:vl-1,M);
    end
end

function write_points(docRoot,P,M)
    assert(size(P,2)==3);
    [Vp,Fp] = readOBJ([path_to_matgp,'/data/mesh/sphereSmall.obj']);
    Vp = Vp * 0.005/0.2;
    for i = 1:size(P,1)
        write_triangle_mesh(docRoot,bsxfun(@plus,Vp,P(i,:)),Fp,[],M);
    end
end

function [M] = material_property(material)
    switch(material)
        case 'default'
            M = struct();
            M.eta = 1.45;
            M.glitterColor = [0 0 0];
            M.glitterSpread = 0.1;
            M.shadeColor = [0.5 0 0];
        case 'blank'
            M = struct();
            M.eta = 1;
            M.glitterColor = [0 0 0];
            M.glitterSpread = 0.0;
            M.shadeColor = [1 1 1];
        case 'grey'
            M = struct();
            M.eta = 1.5;
            M.glitterColor = [0 0 0];
            M.glitterSpread = 0.0;
            M.shadeColor = [1 1 1] * 0.2;
        case 'black'
            M = struct();
            M.eta = 2;
            M.glitterColor = [0 0 0];
            M.glitterSpread = 0.0;
            M.shadeColor = [0 0 0];
        case 'reflect'
            M = struct();
            M.eta = 1.45;
            M.glitterColor = [0 0 0];
            M.glitterSpread = 0.1;
            M.shadeColor = [0.5 0.5 0.5];
    end
end

function [M] = material_property_quiver(rgb)
    M = struct();
    M.eta = 1.45;
    M.glitterColor = [0 0 0];
    M.glitterSpread = 0.1;
    M.shadeColor = rgb;
end

function append_float3(self,docRoot,name,value)
    float_data = docRoot.createElement('float3'); 
    self.appendChild(float_data);
    float_data.setAttribute('name',name);
    float_data.appendChild(docRoot.createTextNode(sprintf('%f %f %f',value(1),value(2),value(3))));
end

function append_float(self,docRoot,name,value)
    float_data = docRoot.createElement('float'); 
    self.appendChild(float_data);
    float_data.setAttribute('name',name);
    float_data.appendChild(docRoot.createTextNode(sprintf('%f',value)));
end

function append_texture3d(self,docRoot,name,src)
    float_data = docRoot.createElement('texture3d'); 
    self.appendChild(float_data);
    float_data.setAttribute('name',name);
    float_data.setAttribute('src',sprintf('%s',src));
end

function write_material(docRoot,TriangleMesh,M,texture)

    material = docRoot.createElement('material'); 
    TriangleMesh.appendChild(material);
    
        code = docRoot.createElement('code'); 
        material.appendChild(code);
        code.appendChild(docRoot.createTextNode('"OBJ"'));
       
        parameters = docRoot.createElement('parameters'); 
        material.appendChild(parameters);
        
        append_float3(parameters,docRoot,'Kd',[1,1,1]*0.8);
        append_float3(parameters,docRoot,'Ks',[1,1,1]*0.3);
        
        if ~isempty(texture.map_d)
            append_texture3d(parameters,docRoot,'map_d',texture.map_d);
        end
        if ~isempty(texture.map_Kd)
            append_texture3d(parameters,docRoot,'map_Kd',texture.map_Kd);
        end
        if ~isempty(texture.map_Ks)
            append_texture3d(parameters,docRoot,'map_Ks',texture.map_Ks);
        end
        if ~isempty(texture.map_Ns)
            append_texture3d(parameters,docRoot,'map_Ns',texture.map_Ns);
        end
        if ~isempty(texture.map_Bump)
            append_texture3d(parameters,docRoot,'map_Bump',texture.map_Bump);
        end      
end

function write_material_old(docRoot,TriangleMesh,M)

    material = docRoot.createElement('material'); 
    TriangleMesh.appendChild(material);
    
        code = docRoot.createElement('code'); 
        material.appendChild(code);
        code.appendChild(docRoot.createTextNode('"MetallicPaint"'));
    
        parameters = docRoot.createElement('parameters'); 
        material.appendChild(parameters);
        
        append_float(parameters,docRoot,'eta',M.eta);      
        append_float3(parameters,docRoot,'glitterColor',M.glitterColor);
        append_float(parameters,docRoot,'glitterSpread',M.glitterSpread);      
        append_float3(parameters,docRoot,'shadeColor',M.shadeColor);
        
%         append_float(parameters,docRoot,'eta',1.45);      
%         append_float3(parameters,docRoot,'glitterColor',[0 0 0]);
%         append_float(parameters,docRoot,'glitterSpread',0.1);      
%         append_float3(parameters,docRoot,'shadeColor',[0.5 0 0]);
    
end

function write_float_array(docRoot,TriangleMesh,name,V)
    V = V';
    V = V(:);
    positions = docRoot.createElement(name); 

    positions.appendChild(docRoot.createTextNode(sprintf('%f ',V)));

    TriangleMesh.appendChild(positions);
end

function write_int_array(docRoot,TriangleMesh,name,F)
    F = F';
    F = F(:);
    triangles = docRoot.createElement(name); 
    triangles.appendChild(docRoot.createTextNode(sprintf('%d ',F)));
    TriangleMesh.appendChild(triangles);
end

function write_positions(docRoot,TriangleMesh,V)
    V = V';
    V = V(:);
    positions = docRoot.createElement('positions'); 

    positions.appendChild(docRoot.createTextNode(sprintf('%f ',V)));

    TriangleMesh.appendChild(positions);
end

function write_triangles(docRoot,TriangleMesh,F)
    F = F';
    F = F(:);
    triangles = docRoot.createElement('triangles'); 
    triangles.appendChild(docRoot.createTextNode(sprintf('%d ',F)));
    TriangleMesh.appendChild(triangles);
end

function write_texcoords(docRoot,TriangleMesh,UV)
    UV = UV';
    UV = UV(:);
    texcoords = docRoot.createElement('texcoords'); 
    texcoords.appendChild(docRoot.createTextNode(sprintf('%d ',UV)));
    TriangleMesh.appendChild(texcoords);
end

function write_texture3d(docRoot,TriangleMesh,tfname)
    texture3d = docRoot.createElement('texture3d'); 
    texture3d.appendChild(docRoot.createTextNode(tfname));
    TriangleMesh.appendChild(texture3d);
end

function write_textured_triangle_mesh(docRoot,V,F,UV,M,texture)

    if isunix
        F = F-1;
    else
    F = F - 1;
    end

    scene = docRoot.getDocumentElement;

    TriangleMesh = docRoot.createElement('TriangleMesh'); 
    scene.appendChild(TriangleMesh);
    TriangleMesh.appendChild(docRoot.createTextNode(sprintf('\n')));
    
    write_positions(docRoot,TriangleMesh,V);
    write_triangles(docRoot,TriangleMesh,F);
    if(~isempty(UV))
        write_texcoords(docRoot,TriangleMesh,UV);
    end
    write_material(docRoot,TriangleMesh,M,texture);
    
    docRoot.appendChild(docRoot.createComment('this is a comment'));

end

function write_triangle_mesh(docRoot,V,F,UV,M)

    F = F - 1;

    scene = docRoot.getDocumentElement;

    TriangleMesh = docRoot.createElement('TriangleMesh'); 
    scene.appendChild(TriangleMesh);
    TriangleMesh.appendChild(docRoot.createTextNode(sprintf('\n')));
    
    write_positions(docRoot,TriangleMesh,V);
    write_triangles(docRoot,TriangleMesh,F);
    if(~isempty(UV))
        write_texcoords(docRoot,TriangleMesh,UV);
    end
    write_material_old(docRoot,TriangleMesh,M);
    
    docRoot.appendChild(docRoot.createComment('this is a comment'));

end


function write_linesegments(docRoot,V,L,M)

    scene = docRoot.getDocumentElement;

    LineSegments = docRoot.createElement('LineSegments'); 
    scene.appendChild(LineSegments);
    LineSegments.appendChild(docRoot.createTextNode(sprintf('\n')));
    
    write_float_array(docRoot,LineSegments,'positions',V);
    write_int_array(docRoot,LineSegments,'indices',L);
    
    write_material_old(docRoot,LineSegments,M);
    
    docRoot.appendChild(docRoot.createComment('this is a comment'));

end

function write_camera2(docRoot)

    scene = docRoot.getDocumentElement;

    Camera = docRoot.createElement('PerspectiveCamera'); 
    scene.appendChild(Camera);
    %Camera.appendChild(docRoot.createTextNode(sprintf('\n')));
    
    Camera.setAttribute('from','0 1 -5');
    Camera.setAttribute('to','0 0 0');
    Camera.setAttribute('up','0 1 -0.2');
    Camera.setAttribute('fov','37');
end

function write_camera(docRoot)
    C = {{'from','0 1.2 -6'},...
    {'to','0 0 0'},...
    {'up','0 1 -0.2'},...
    {'fov','11'}};
%     C = {{'from','0 5 0'},...
%     {'to','0 0 0'},...
%     {'up','1 0 0'},...
%     {'fov','11'}};
    C = {{'from','0 0 -5'},...
    {'to','0 0 0'},...
    {'up','0 1 0'},...
    {'fov','14'}};
% donut & dirac
%     C = {{'from','0 3 -6'},...
%     {'to','0 0 0'},...
%     {'up','0 1 -0.5'},...
%     {'fov','11'}};
    C = {{'from','0 2 -6'},...
    {'to','0 0 0'},...
    {'up','0 1 -0.2'},...
    {'fov','11'}};
    write_element_with_attributes(docRoot,'PerspectiveCamera',C,[]);
end

function write_light(docRoot)
    D = {{'L','','1.2 1.2 1.2'}};
    write_element_with_attributes(docRoot,'AmbientLight',[],D);
end

function append_child(self,docRoot,type,name,value)
    float_data = docRoot.createElement(type); 
    self.appendChild(float_data);
    float_data.setAttribute('name',name);
    float_data.appendChild(docRoot.createTextNode(value));
end

function write_element_with_attributes(docRoot,name,name_value_strings,children)

    scene = docRoot.getDocumentElement;

    Camera = docRoot.createElement(name); 
    scene.appendChild(Camera);
    
    for i=1:length(name_value_strings)
        Camera.setAttribute(name_value_strings{i}{1},name_value_strings{i}{2});
    end
    
    for i=1:length(children)
        append_child(Camera,docRoot,children{i}{1},children{i}{2},children{i}{3});
    end
end