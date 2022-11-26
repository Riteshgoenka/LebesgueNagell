N:=3^3*5^2;
 
G3<a3>:=FullDirichletGroup(3);
G5<a5>:=FullDirichletGroup(5);

G:=DirichletGroup(N, CyclotomicField(4));

eps:=G!(a3*a5);

listk:=[1,2,4,5,3,6,7,10];

/* Define base field. */
R<x>:=PolynomialRing(RationalField());
K1<z3>:=NumberField(x^2-6);
R1<X>:=PolynomialRing(K1);
K<z1>:=NumberField(X^2+1);
zeta_4:=z1;

load "all-675.data";

IsBaseField:=function(qexp)
for i:=1 to #qexp do
  if (qexp[i] notin K) then
    return false,i;
  end if;
end for;
return true;
end function;

/* Build indexed prime list. */
Bpl:=1000;
pl:=[2];
p:=NextPrime(pl[#pl]);
while (p lt Bpl) do
  pl:=Append(pl,p);
  p:=NextPrime(pl[#pl]);
end while;

AddSeq:=function(s)
l:={};
for i:=1 to #s do
  l:=l join {s[i]};
end for;
return l;
end function;


Boundq:=function(qexp,i,F)
 q:=pl[i];
 aq:=qexp[i];
 S:={2,3,q};

 // Bounds from good primes q.
 if (KroneckerSymbol(12,q) eq 1) and (KroneckerSymbol(8,q) eq 1) then
   q4:=Floor(SquareRoot(2*q));
   for u:=-q4 to q4 do
     S:=S join AddSeq(PrimeDivisors(IntegerRing()!AbsoluteNorm((F!(u))-aq)));
   end for;
 elif (KroneckerSymbol(12,q) eq -1) and (KroneckerSymbol(8,q) eq 1) then
   q4:=Floor(SquareRoot(2*q));
   for u:=-q4 to q4 do
     S:=S join AddSeq(PrimeDivisors(IntegerRing()!AbsoluteNorm((F!(u*z1))-aq)));
   end for;
 elif (KroneckerSymbol(12,q) eq 1) and (KroneckerSymbol(8,q) eq -1) then
   q4:=Floor(SquareRoot(2*q));
   for u:=-q4 to q4 do
     S:=S join AddSeq(PrimeDivisors(IntegerRing()!AbsoluteNorm((F!(u*z3))-aq)));
   end for;
 elif (KroneckerSymbol(12,q) eq -1) and (KroneckerSymbol(8,q) eq -1) then
   q4:=Floor(SquareRoot(2*q));
   for u:=-q4 to q4 do
     S:=S join AddSeq(PrimeDivisors(IntegerRing()!AbsoluteNorm(F!(u*z1*z3)-aq)));
   end for;
 end if;

 // Bounds from bad primes q.
 epsmq:=(eps^(-1))(q);
 S:=S join AddSeq(PrimeDivisors(IntegerRing()!AbsoluteNorm(F!epsmq*(q+1)^2-aq^2)));

 return S;
end function;


BM:=10;
Boundf:=function(qexp)
first:=true;
F:=Parent(qexp[1]);
listq:={};

// Use first BM primes q to get a bound.
for i:=1 to BM do
  if (qexp[i] notin K) and (pl[i] notin {2,3}) then
    listq:=listq join {pl[i]};
    J:=Boundq(qexp,i,F);
    if first then
      S:=J;
      first:=false;
    else
      S:=S meet J;
    end if;
  end if;
end for;
return S,listq;
end function;


print "Compute the bounds on p from forms whose coefficient field is not M_beta.";

T:={};
L:={};
BL:={};
for k in listk do
  print(qexplist[k][1][1]);
  mxl:=Factorization(PolynomialRing(K)!qexplist[k][1][1]);

  for j:=1 to #mxl do
  B1:={};
  mx:=mxl[j][1];
  qexpx:=qexplist[k][2];

  W<a>:=NumberField(mx);

  qexp:=[];
  for i:=1 to #qexpx do
    qexp:=Append(qexp,Evaluate(qexpx[i],a));
  end for;

  if not IsBaseField(qexp) then
    print "Form #",k,": absolute degree of coefficient field",Degree(qexplist[k][1][1]),": conjugate #",j;
    B,Q:=Boundf(qexp);
    B1:=B1 join B;
    T:=T join B;
    L:=L join Q;
  end if;
  end for;
  if B1 ne {} then
    BL:=BL join {k};
  end if;
end for;

print "Primes in bound",T;
print "Primes used in bound",L;
print "Forms eliminated",BL;