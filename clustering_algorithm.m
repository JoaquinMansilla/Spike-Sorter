
function[IDX,C]=clustering_algorithm(FeaturePoints,K,value,ListC)
%DO NOT EDIT{
%Put your Clustering Algorithm below.
%K is the number of clusters.
%Feature Points is the matrix of Feature points. Rows of the matrix
%correspond to points, columns correspond to variables.
%IDX is a vector that indicates each point belongs to # cluster.
%C is a  K-by-P matrix that contains k cluster centroids points.  
%The function MUST have this format:
%"function[]=name_of_function()%##name_of_function##%" with no blank spaces, this
%helps the program to recognize the functions 
aux=ListC(value,:);
aux=strtrim(aux);
aux=str2func(aux);
[IDX,C]=aux(FeaturePoints,K);

%Create functions for clustering down HERE:

%Examples

function[IDX,C]=KMEANS(FPoints,K)%##KMEANS##%
[IDX,C]=kmeans(FPoints,K);

function[IDX,C]=Template_Matching(FPoints,K)%##Template_Matching##%
global WFmatrix
IDX=[];
Unitn=0;
WFmatrixaux=double(WFmatrix);
if (exist('Spikes.mat','file')==2)
    auxcd=cd;
    path=which('SpikeSort.m');
    path=path(1:(end-11));
    cd(path);
    load('Spikes.mat');
    sizeS=size(Spikes);
    sizeW=size(WFmatrix);
    for i=1:sizeW(1)
        d=inf;
        for j=1:sizeS(1)
            AutoCorr=Spikes(j,:)*Spikes(j,:)';
            Corr=WFmatrixaux(i,:)*Spikes(j,:)';
            daux=Corr/AutoCorr;
            daux=abs(daux-1);
            if daux<d
                aux=j;
                d=daux;
            end
        end
        AutoCorr=Spikes(aux,:)*Spikes(aux,:)';
        if(((WFmatrixaux(i,:)*Spikes(aux,:)')>0.8*AutoCorr) && ((WFmatrixaux(i,:)*Spikes(aux,:)')<1.2*AutoCorr))
            auxIDX(i)=aux;
        end
    end
    uIDX=unique(auxIDX);
    [N,X]= hist(auxIDX,uIDX);
    if X(1)==0
        X=X(2:end);
        N=N(2:end);
        uIDX=uIDX(2:end);
    end 
    IDX=zeros(sizeW(1),1);
    for i=1:length(uIDX)
        if N(i)>100
            Unitn=Unitn+1;
            aux=find(auxIDX==X(i));
            IDX(aux)=Unitn;
            C(Unitn,:)=Spikes(X(i),:);
        end
    end
else
    errordlg('No Spikes Record Found','Error')
end
cd(auxcd);

function[IDX,C]=algorithm2(FPoints,K)%##algorithm2##%

function[IDX,C]=algorithm3(FPoints,K)%##algorithm3##%
