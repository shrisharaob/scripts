function [trial, units, pfPars, par] = LoadTestData(varargin)
[filebase, trialName ] = DefaultArgs(varargin, {'jg05-20120315', 'crt1'});
% filebase = 'jg05-20120315';
session = MTASession(filebase,[],'cof');
analysisFolderPath = '~/data/analysis/';
fullpath = ['/data/homes/shrisha/data/nlx/jg05-20120315/' filebase];
par = LoadPar(fullpath);
trialnames = session.list_trialNames;
trial = MTATrial(filebase, {{'CluRes', session.sampleRate}}, trialName);
% units = load([analysisFolderPath  filebase '/' filebase '.SelectedPyrCells.mat']);
% % pfPars = [];
% pfPars = load(['~/data/analysis/' filebase '/' filebase '.FindPFPars.' trialName '.mat']);
end

%ss =tr.Bhv.getState('walk').state;
%ss = ss + 1; % bug xyz index starts from 0
% trajectory = [];
% for  i = 1 : size(ss,1)
% temp = sq(tr.xyz([ss(i,1) : ss(i,2)], 7 ,[1,2])); 
%    trajectory = [ trajectory ; temp];
% end



%         