include("graph.jl")
include("core.jl")
include("alg1.jl")
include("schur.jl")

using LinearAlgebra
using Laplacians

fname = open("filename.txt", "r")
str   = readline(fname);
nn     = parse(Int, str);

for nnnn=1:nn
#######S0 \cup S1

str = readline(fname);
str = split(str);
G   = get_graph(str[1]);
on=G.n;om=G.m;
Gc=findconnect(G)
G=Gc;
n=G.n;m=G.m;

beta=0.05;
eps=1;
p=Int(round(log(n)/eps^2));
p=100;
len=Int(round(1/beta^2));
len=300;
selc_num=100;
kmax=3;

fans1=open("fans2.txt","a");
println(fans1,str[1],' ',on,' ',om,' ',n,' ',m)
println(fans1,"p=",p,",len=",len,",selc_num=",selc_num);



################Init
L=lapsp(G);
a=adjsp(G)

d=zeros(n);
d=[L[i,i] for i=1:n];
S=[];
for i=1:selc_num
    push!(S,argmax(d))
    d[argmax(d)]=0;
end

G=get_graph_union(G,S)
n=G.n;m=G.m;
L=lapsp(G);
a=adjsp(G);
d=zeros(n);
d=[L[i,i] for i=1:n];
S=[];
push!(S,argmax(d))


V=union(1:n)
F=setdiff(V,S)
t1=time()

W,We,rat,lavg=init(G,S,p,len);
t2=time()
println(fans1,"time:",t2-t1,",ratio:",rat,",average length:",lavg)
close(fans1)

f=approxchol_lap(a)
yi=ones(n)
yi[S].=0;
h=f(yi)



######### add one ##########

S0=S;S1=[];
t3=time()
xxx=0;
for i=1:kmax
    fans1=open("fans2.txt","a");
    tt1=time()
    xxx,W,We=addTerminal(S0,S1,G,W,We,p,h)
    push!(S1,xxx)
    tt2=time()
    println(fans1,i,' ',tt2-tt1)
    close(fans1)
end
t4=time()

fans1=open("fans2.txt","a");
for i=1:kmax
    println(fans1,getans(G,S0,S1[1:i]))
end
close(fans1)




S_deg=zeros(Int,kmax)
d[S0].=0
for i=1:kmax
    S_deg[i]=argmax(d)
    d[argmax(d)]=0
end


bet=BetweennessCentrality(G)
clo=ClosenessCentrality(G)
pag=PageRankCentrality(G)

S_bet=[]
S_clo=[]
S_pag=[]
bet[S0].=0
clo[S0].=0
pag[S0].=0
for i=1:kmax
    dd=argmax(bet)
    push!(S_bet,dd)
    bet[dd]=0;
    dd=argmax(clo)
    push!(S_clo,dd)
    clo[dd]=0;
    dd=argmax(pag)
    push!(S_pag,dd)
    pag[dd]=0;
end

texa1=time()
S_acc=[]
FF=union(1:n)
setdiff!(FF,S0)
for i=1:kmax
    sel=0;tmpmax=0;
    for j in FF
        xx=getans(G,S0,union(S_acc,j))
        if xx>tmpmax
            tmpmax=xx
            sel=j
        end
    end
    push!(S_acc,sel)
    setdiff!(FF,sel)
end
texa2=time()

S_hic=[]
hic_sco=zeros(n);
for j in F
    hic_sco[j]=getans(G,S0,union(j))
end
for i=1:kmax
    ddd=argmax(hic_sco)
    push!(S_hic,ddd)
    hic_sco[ddd]=0;
end



fans1=open("fans2.txt","a");
println(fans1,"exacttime:",texa2-texa1,"fasttime:",t4-t3)
for i=1:kmax
    println(fans1,getans(G,S0,S1[1:i]),' ',getans(G,S0,S_deg[1:i]),' ',getans(G,S0,S_bet[1:i]),' ',getans(G,S0,S_clo[1:i]),' ',getans(G,S0,S_pag[1:i]),' ',getans(G,S0,S_hic[1:i]),' ',getans(G,S0,S_acc[1:i]))
end
println(fans1)
close(fans1)


end


close(fname)
