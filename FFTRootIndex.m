function Index=FFTRootIndex(f,tolFFt)
if nargin<2
    tolFFt=0;
end
V=CZifft(f);
if min(V(:))<0
    error('~ (f>0)')
end
V=V.^0.5;
g=CZfft(V,tolFFt);
[a,b]=find(g);
[~, sortb_Index]=sort(abs(b),'descend');
  Index= a(sortb_Index,:);
end