function[auxFeature]=Features(Feat)
%DO NOT EDIT{
%Put your Featute Algorithm below.
%Feature Points is the matrix of Feature points. Rows of the matrix
%correspond to points, columns correspond to variables.
%The function MUST have this format:
%"function[]=name_of_function()%##name_of_function##%" with no blank spaces, this
%helps the program to recognize the functions}
%Change the function call after including new input and output function
%variables.
global ListF
aux=ListF(Feat,:);
aux=strtrim(aux);
aux=str2func(aux);
auxFeature=aux();
%DO NOT EDIT{
%Create Feature functions down HERE:}
%Examples


function[auxFeature]=PC1()%##PC1##%
global WFmatrix
auxWFmatrix=double(WFmatrix);
[COEFF, SCORE] = princomp(auxWFmatrix,'econ');
%[COEFF, SCORE] = pca(WFmatrix','NumComponents',4); alternative option
%(faster, available on new Matlab versions
auxFeature=[SCORE(:,1)];


function[auxFeature]=PC2()%##PC2##%

global WFmatrix
auxWFmatrix=double(WFmatrix);
[COEFF, SCORE] = princomp(auxWFmatrix,'econ');
%[COEFF, SCORE] = pca(WFmatrix','NumComponents',4); alternative option
%(faster, available on new Matlab versions
auxFeature=[SCORE(:,2)];


function[auxFeature]=PC3()%##PC3##%

global WFmatrix
auxWFmatrix=double(WFmatrix);
[COEFF, SCORE] = princomp(auxWFmatrix,'econ');
%[COEFF, SCORE] = pca(WFmatrix','NumComponents',4); alternative option
%(faster, available on new Matlab versions
auxFeature=[SCORE(:,3)];


function[auxFeature]=PC4()%##PC4##%

global WFmatrix
auxWFmatrix=double(WFmatrix);
[COEFF, SCORE] = princomp(auxWFmatrix,'econ');
%[COEFF, SCORE] = pca(WFmatrix','NumComponents',4); alternative option
%(faster, available on new Matlab versions
auxFeature=[SCORE(:,4)];


function[auxFeature]=Slice1()%##Slice1##%
global WFmatrix
global handles2Spk
global Fs
x=str2double(get(handles2Spk.Slice_1,'String'));
x=round(x*Fs/1000);
if x==0
    errordlg('No Slice Found','Error')
    auxFeature=zeros(size(WFmatrix(:,1)));
else
    auxFeature=(WFmatrix(:,x));
end


function[auxFeature]=Slice2()%##Slice2##%
global WFmatrix
global handles2Spk
global Fs
x=str2double(get(handles2Spk.Slice_2,'String'));
x=round(x*Fs/1000);
if x==0
    errordlg('No Slice Found','Error')
    auxFeature=zeros(size(WFmatrix(:,1)));
else
    auxFeature=(WFmatrix(:,x));
end


function[auxFeature]=CustomFeature1()%##CustomFeature1##% 


function[auxFeature]=CustomFeature2()%##CustomFeature2##%

