function RunTree(varargin) % Version with symmetrical iTEBD
% If you change the folder name, make corresponding changes in
% line-9-Run.m, line-4-Exact_2DIsing.m

var=[2:10];
varName='dc';

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

Nvar=numel(var);
EXP=[MainPath,ParaFolder,ParaFile];

for n=1:Nvar
    para=Parameter;  %#ok<NASGU>
    eval(['para.',varName,'=',num2str(var(n)),';']);
    save([EXP,'_',num2str(n)],'para');
end

for n=1:Nvar
    if(numel(varargin)==0)
        Main([EXP,'_',num2str(n)]);
    elseif(numel(varargin)==1)
        if(strcmp(varargin{1},'L'))
            
        else
            error('Bad input of Run.m');
        end
        
    else
        error('Bad input of Run.m');
    end
end
Alert
end
