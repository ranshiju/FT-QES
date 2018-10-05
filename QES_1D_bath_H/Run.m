function Run(varargin) % Version with symmetrical iTEBD
% If you change the folder name, make corresponding changes in
% line-9-Run.m, line-4-Exact_2DIsing.m

if(isunix)
    addpath(genpath('../ClosedFunctions'));
end
addpath(genpath('.'))

% Jbath
var{1}=[0.03,0.001,0.005];
% hbath
var{2}=[0.01,0.025:0.025:1,0.505:0.005:0.52,0.53,0.54]; 
varName={'para.hx2','para.hx'};

if(isunix)
    MainPath='/home/sran/Dropbox/MatlabCode/MPO_Disentangler_Ising/';
    ParaFolder='para/';
else
    MainPath='C:\Users\sran\Dropbox\MatlabCode\MPO_Disentangler_Ising\';
    ParaFolder='para\';
end
ParaFile=['Para_',datestr(datetime,'mm-dd-HH-MM-SS')];

if(exist([MainPath,ParaFolder],'dir')==0)
    mkdir([MainPath,ParaFolder]);
end

Nvar1=numel(var{1});
Nvar2=numel(var{2});
EXP=[MainPath,ParaFolder,ParaFile];

for n1=1:Nvar1
    for n2=1:Nvar2
        para=Parameter;  %#ok<NASGU>
        paraRC=ParameterRC;  %#ok<NASGU>
        paraED=ParameterED;  %#ok<NASGU>
        eval([varName{1},'=',num2str(var{1}(n1)),';']);
        eval([varName{2},'=',num2str(var{2}(n2)),';']);
        save([EXP,'_',num2str(n1),num2str(n2)],'para','paraRC','paraED');
    end
end

for n1=1:Nvar1
    for n2=1:Nvar2
        if(numel(varargin)==0)
            Main([EXP,'_',num2str(n1),num2str(n2)]);
        elseif(numel(varargin)==1)
            if(strcmp(varargin{1},'L'))

            else
                error('Bad input of Run.m');
            end

        else
            error('Bad input of Run.m');
        end
    end
end
Alert
end