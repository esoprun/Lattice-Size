//This is the implementation of the basis reduction algorithm from 
//[1] Lattice Size and Generalized Basis Reduction in Dimension 3, by A. Harrison and 
//J. Soprunova  Discrete Comput Geom (2021), https://doi.org/10.1007/s00454-020-00269-x
//This algorithm computes the lattice size of a plane lattice polygon P, as defined in 
//[2] Lattice Size of Plane Convex Bodies by Harrison, Soprunova, Tierney, 
//https://arxiv.org/pdf/1709.03451.pdf
//We are grateful to Ivan Soprunov for translating our MATLAB code to Magma.

function LatticeWidth(P,h);
//finds the lattice Width of a lattice polygon P in the direction of h\in\Z^2
maxP:=Maximum([v[1]*h[1]+v[2]*h[2] : v in P]);
minP:=Minimum([v[1]*h[1]+v[2]*h[2] : v in P]);
return maxP-minP;
end function;

function ReductionStep(P);
//This step is based on Proposition 2.7 and Theorem 2.8 from [1]
h1:=<1,0>;
h2:=<0,1>;

wid1:=LatticeWidth(P,h1);
wid2:=LatticeWidth(P,h2);
//Next ensure that |h1|<=|h2|.
if wid2 lt wid1 then 
    h2:=<1,0>;
    h1:=<0,1>;
    wid1:=LatticeWidth(P,h1);
    wid2:=LatticeWidth(P,h2);
end if;

//Find i in Z such that |ih1+h2|-> min. First, check four values of i=\pm floor(|h2|/|h1|)
//and \pm ceiling(|h2/h1|) and find such m with smallest value of |ih1+h2|. 
//See the proof of Theorem 2.8.
M:=[Floor(wid2/wid1),Ceiling(wid2/wid1),-Floor(wid2/wid1),-Ceiling(wid2/wid1)];
minW,j:=Minimum([LatticeWidth(P,<i*h1[1]+h2[1],i*h1[2]+h2[2]>) : i in M]);
if minW lt wid1 then // Now we use Proposition 2.7
   if minW lt (2/3)*wid2 then
       return <M[j]*h1[1]+h2[1],M[j]*h1[2]+h2[2]>,h1,0;
   else // Find two vectors of smallest norm among {h1, h1+h2,h1-h2,2h1+h2, 2h1-h2}. 
       FiveVectors:=[h1,<h1[1]+h2[1],h1[2]+h2[2]>,<h1[1]-h2[1],h1[2]-h2[2]>,<2*h1[1]+h2[1],2*h1[2]+h2[2]>,<2*h1[1]-h2[1],2*h1[2]-h2[2]>];
       MinFive1,j:=Minimum([LatticeWidth(P,f) : f in FiveVectors]);
       FourVectors:=Remove(FiveVectors,j);
       MinFive2,k:=Minimum([LatticeWidth(P,f) : f in FourVectors]);
       return FiveVectors[j],FourVectors[k],1;
    end if;
else 
    m:=Floor(2*wid2/wid1); //This is based on Theorem 2.8: to find min|i*h1+h2| over i in Z
                          //one only needs to consider i with |i|<=Floor(2|h2|/|h1|).
    minW,j:=Minimum([LatticeWidth(P,<i*h1[1]+h2[1],i*h1[2]+h2[2]>) : i in [-m..m]]);
    return h1, <(j-m-1)*h1[1]+h2[1],(j-m-1)*h1[2]+h2[2]>, 1;
end if;
end function;

// P:=[<0,0>,<3,5>,<7,9>,<8,12>]; [h1,h2]=[<-1,1>,<1,0>]
// P:=[<0,0>,<2,3>,<2,7>,<4,8>]; [h1,h2]=[<1,-2>,<0,1>]

function LatticeSize(P);
terminate:=0;
count:=0; //number of iterations
A:=[<1,0>,<0,1>];
while terminate eq 0 do
    h1,h2,terminate:=ReductionStep(P);
    P:=[<h1[1]*p[1]+h1[2]*p[2],h2[1]*p[1]+h2[2]*p[2]> : p in P];//Apply matrix [h1;h2] to P
    A:=[<h1[1]*a[1]+h1[2]*a[2],h2[1]*a[1]+h2[2]*a[2]> : a in A];//Apply matrix [h1;h2] to A
    count:=count+1;
end while;
//find l1,l2,l3,l4. See the definition 2.1 from Lattice Size of Plane Convex Bodies by 
//Harrison, Soprunova, Tierney
L:=[Maximum([p[1]+p[2] : p in P])-Minimum([p[1] : p in P])-Minimum([p[2] : p in P]),
-Minimum([p[1]+p[2] : p in P])+Maximum([p[1] : p in P])+Maximum([p[2] : p in P]),
Maximum([p[1]-p[2] : p in P])-Minimum([p[1] : p in P])+Maximum([p[2] : p in P]),
Maximum([-p[1]+p[2] : p in P])+Maximum([p[1] : p in P])-Minimum([p[2] : p in P])];
l,j:=Minimum(L);
if j eq 1 then A:=[<A[1][1],A[2][1]>,<A[1][2],A[2][2]>]; end if;
if j eq 2 then A:=[<-A[1][1],-A[2][1]>,<-A[1][2],-A[2][2]>]; end if;
if j eq 3 then A:=[<A[1][1],A[2][1]>,<-A[1][2],-A[2][2]>]; end if;    
if j eq 4 then A:=[<-A[1][1],-A[2][1]>,<A[1][2],A[2][2]>]; end if;
return l, A, count;//l=ls(P), A is the matrix that computes ls(P)
end function;


// EXPERIMENTS:

// Using reduction algorithm:

//for i in [1..100] do
//    n:=Random([5..20]);
//    Q:=RandomPolytope(2,n,100);
//    V:=Vertices(Q);
//    P:=[<Vector(v)[1],Vector(v)[2]> : v in V];
//    time LatticeSize(P);
//end for;


// Using Castryck-Cools (onion skins) algorithm available at 
//https://homes.esat.kuleuven.be/~wcastryc/code/basic_commands.m

//for i in [1..100] do
//    n:=Random([5..20]);
//    Q:=RandomPolytope(2,n,100);
//    V:=Vertices(Q);
//    P:=[<Integers()!Vector(v)[1],Integers()!Vector(v)[2]> : v in V];
//    LatticePolytope(P);
//    time LatticeSizeRecursiveSigma(LatticePolytope(P));
//end for;    
