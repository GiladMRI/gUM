function Out=gScoreImageSimilarity(A,AbsB,BDx,BDy)
A=A*grmss(AbsB)/grmss(A);
ArecoXA=abs(A);
RecDx=diff(ArecoXA,1,1);
RecDy=diff(ArecoXA,1,2);

RMSScr=grmss(ArecoXA-AbsB);
DScr=grmss(cat(5,RecDx-BDx,(RecDy-BDy).'));

Out=RMSScr+DScr;
