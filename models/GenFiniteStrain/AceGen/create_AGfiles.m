clc
close all
clear all

folder = 'AGfiles';

FEpos = 5;
EVpos = 9;
SCpos = 13;

% dep       bs  fe  ev  sc
EFdep   = [ 1,  1,  0,  1 ];
EoutFdep= [ 0,  1,  0,  0 ];
RFdep   = [ 1,  1,  1,  1 ];  %Has to be dependent on all as this is used as basis!
dRdXdep = [ 1,  1,  1,  1 ];
dRdFdep = [ 1,  1,  1,  1 ];

% typ = [bs,fe,ev,sc]
% Name of module file
AGname    = @(typ) ['GFS_AG'   '_FE' num2str(typ(2)) '_EV' num2str(typ(3)) '_SC' num2str(typ(4)) '.f90'];
% Name for created acegen subroutine files
RFname    = @(typ) ['RF'   num2str(typ(1)) '_FE' num2str(typ(2)) '_EV' num2str(typ(3)) '_SC' num2str(typ(4)) '.f90'];
dRdXname  = @(typ) ['dRdX' num2str(typ(1)) '_FE' num2str(typ(2)) '_EV' num2str(typ(3)) '_SC' num2str(typ(4)) '.f90'];
dRdFname  = @(typ) ['dRdF' num2str(typ(1)) '_FE' num2str(typ(2)) '_EV' num2str(typ(3)) '_SC' num2str(typ(4)) '.f90'];
EFname    = @(typ) ['EF' num2str(typ(1)) '_FE' num2str(typ(2)) '_SC' num2str(typ(4)) '.f90'];
EoutFname = @(typ) ['EoutF' '_FE' num2str(typ(2)) '.f90'];


%Strategy
% A: Find all relevant files and group their names into a cell array
% 1) Consider all RF's and find which versions we have
% 2) For all of these version, we also need equivalent dRdX and dRdF
% 3) For some of these version, we need to use EF and EoutF multiple times
% B: Create one of each filetype depending on FE, EV and SC

%% Get file names
% Get filenames
cdir = cd(folder);  files = dir;    cd(cdir);
files=files(3:end); files = struct2cell(files); files = files(1, 1:end);

% Setup functions to identify different filetypes
find_RF     = @(str) findprefix(str,'RF');
find_dRdX   = @(str) findprefix(str,'dRdX');
find_dRdF   = @(str) findprefix(str,'dRdF');
find_EF     = @(str) findprefix(str,'EF');
find_EoutF  = @(str) findprefix(str,'EoutF');

% Find different filetypes
RF_files    = files(cellfun(find_RF, files));
dRdX_files  = files(cellfun(find_dRdX, files));
dRdF_files  = files(cellfun(find_dRdF, files));
EF_files    = files(cellfun(find_EF, files));
EoutF_files = files(cellfun(find_EoutF, files));

%% Find all different variants (based on RF)
% Find all variants using RF
Ntypes = length(RF_files);
RF = reshape(cell2mat(RF_files), [], Ntypes)';
BS = str2num(RF(:,3));  Nback = max(BS);    %Back-stresses
FE = str2num(RF(:,2+FEpos));                %Free energies
EV = str2num(RF(:,2+EVpos));                %Evolution equations
SC = str2num(RF(:,2+SCpos));                %Sign convention
types = [BS, FE, EV, SC];

ftypes = unique(types(:,2:end), 'rows');    %Find unique combination excl backstress
Nfiles = length(ftypes(:,1));
%fbacks: Each row gives info on which backstresses are available (1 or 0)
fbacks = false(Nfiles, Nback);  
for i = 1:Nfiles
    cinds = all((ones(Ntypes,1)*ftypes(i,:))==types(:,2:end),2);
    tmp = 1:length(cinds); cinds = tmp(cinds);
    for j = 1:Nback
        fbacks(i,j) = any(BS(cinds)==j);
    end
end

%% Check that all required files are available
missing = false(Nfiles);
for i = 1:Nfiles
    for j = 1:Nback
        if fbacks(i,j)
            % Check residual function (not really nec)
            tmpfun = @(str) strcmp(RFname([j,ftypes(i,:)]), str);
            if ~any(cellfun(tmpfun, RF_files))
                fprintf('Residual function missing\n');
                fprintf('bs=%u, fe=%u, ev=%u, sc=%u\n', [j,ftypes(i,:)]);
                missing(i) = true;
            end
            
            % Check dRdX
            tmpfun = @(str) strcmp(dRdXname([j,ftypes(i,:)]), str);
            if ~any(cellfun(tmpfun, dRdX_files))
                fprintf('dRdX missing: \t');
                fprintf('bs=%u, fe=%u, ev=%u, sc=%u\n', [j,ftypes(i,:)]);
                missing(i) = true;
            end
            
            % Check dRdF
            tmpfun = @(str) strcmp(dRdFname([j,ftypes(i,:)]), str);
            if ~any(cellfun(tmpfun, dRdF_files))
                fprintf('dRdF missing: \t');
                fprintf('bs=%u, fe=%u, ev=%u, sc=%u\n', [j,ftypes(i,:)]);
                missing(i) = true;
            end
            
            % Check EF
            tmpfun = @(str) strcmp(EFname([j,ftypes(i,:)]), str);
            if ~any(cellfun(tmpfun, EF_files))
                fprintf('EF missing: \t');
                fprintf('bs=%u, fe=%u, ev=%u, sc=%u\n', [j,ftypes(i,:)]);
                missing(i) = true;
            end
            
            % Check EoutF
            tmpfun = @(str) strcmp(EoutFname([j,ftypes(i,:)]), str);
            if ~any(cellfun(tmpfun, EoutF_files))
                fprintf('EoutF missing: \t');
                fprintf('bs=%u, fe=%u, ev=%u, sc=%u\n', [j,ftypes(i,:)]);
                missing(i) = true;
            end
        end
    end
end

% if any(missing)
%     if ~strcmp(input('There were files missing, continue? (y/n)', 's'), 'y')
%         return
%     end
% end


%% Create module files
for i = 1:Nfiles
    ctype = ftypes(i,:);    %Current file type
    back = 1:Nback;    back(~fbacks(i,:)) = [];
    fid = fopen(AGname([0,ftypes(i,:)]), 'w');
    
    % Write header
    fprintf(fid, '!AceGen automatic generated file for GFS\n');
    fprintf(fid, '!Free Energy:    %u\n', ftypes(i,1));
    fprintf(fid, '!Evolution Law:  %u\n', ftypes(i,2));
    fprintf(fid, '!SignConvention: %u\n', ftypes(i,3));
    fprintf(fid, '!Allowed backstresses: '); fprintf(fid, '%u, ', back); fprintf(fid,'\n');
    if missing(i)
        fprintf(fid, '!WARNING: There were files missing when this module was created\n');
    end
    
    % Module header
    fprintf(fid, 'module GFS_acegen\n');
    fprintf(fid, 'implicit none\n');
    fprintf(fid, '    contains\n');
    
    
    
    for j = 0:Nback
        typ = [j, ftypes(i,:)];
        if j>0
            fbt = fbacks(i,j);
            files = {RFname(typ), dRdXname(typ), dRdFname(typ), EFname(typ)};
            fprintf(fid, '\n');
        else
            fbt = true;
            files = {EoutFname(typ)};
        end
        if fbt
%             typ = [j, ftypes(i,:)];
%             files = {RFname(typ), dRdXname(typ), dRdFname(typ), EFname(typ)};
            
            for k = 1:length(files)
                fexists = true;
                fname = files{k};
                try    
                    fstr = fileread([folder '/' fname]);
                catch
                    fexists = false;
                end
                if fexists
                    % Remove header
                    ind = regexp(fstr(1:1000), ['SUBROUTINE ' fname(1:2)]);
                    fstr = fstr(ind:end);
                    
                    % Remove internal variable v in call
                    ind = regexp(fstr(1:100), '(v,');
                    fstr((ind+1):(ind+2)) = [];
                    
                    % Add subroutine to 'end' to become 'end subroutine'
                    bpoint = length(fstr)-100;  %Break point
                    tmp = fstr(bpoint:end);
                    ind = regexp(tmp, 'END', 'end');
                    tmp = [tmp(1:ind) ' SUBROUTINE'];
                    fstr = [fstr(1:bpoint-1), tmp];
                    fprintf(fid, '%s\n\n', fstr);
                end
            end
        end
    end
    fprintf(fid, 'end module\n');
    fclose(fid);
end