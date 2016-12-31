function [mesh_frames] = load_mesh_frames(MESH_FILE_TYPE,varargin)
% TODO: have not debugged yet.

    use_my = false;

    ii = 1;
    for ii=1:numel(varargin)
       if(strcmp(varargin{ii},'use_my'))
         use_my = true;
       else
         error('Unknown Switch!\n');
       end
       ii = ii + 1;
    end

    MESH_FILE_PATH = '';
    fnames = dir([MESH_FILE_PATH,['*.',MESH_FILE_TYPE]]);
    numfids = length(fnames);
    mesh_frames = cell(numfids,2);
    if(strcmp(MESH_FILE_TYPE,'mesh')==0)
        mesh_frames = cell(numfids,3);
    end
    
    if(use_my)
        for i = 1:numfids
           switch(MESH_FILE_TYPE)  
               case 'mesh'
                   [mesh_frames{i,1},mesh_frames{i,2},mesh_frames{i,3}] = my_readMESH([MESH_FILE_PATH,fnames(i).name]);  
           end
        end
    else
        for i = 1:numfids
           switch(MESH_FILE_TYPE)
               case 'off'
                   [mesh_frames{i,1},mesh_frames{i,2}] = readOFF([MESH_FILE_PATH,fnames(i).name]);
               case 'obj'
                   [mesh_frames{i,1},mesh_frames{i,2}] = readOBJ([MESH_FILE_PATH,fnames(i).name]);
               case 'ply'
                   [mesh_frames{i,1},mesh_frames{i,2}] = readPLY([MESH_FILE_PATH,fnames(i).name]);   
               case 'mesh'
                   [mesh_frames{i,1},mesh_frames{i,2},mesh_frames{i,3}] = readMESH([MESH_FILE_PATH,fnames(i).name]);    
           end
        end
    end