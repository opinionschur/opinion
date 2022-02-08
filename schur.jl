include("graph.jl")


function init(G,T,p,len)
    n=G.n;
    m=G.m;
    t=size(T)[1];
    #beta=t/m;

    # a=spzeros(n,n);
    # a[T,T].+=1;
    wlks=0;
    wlkedge=0;
    wend=0;
    wuse=0;
    wgood=0;
    wleng=0;
    wall=0;
    W=[];
    We=[];
    xyinT=zeros(Int,n)
    #for i in T
        xyinT[T].=1;
    #end
    # I=[];J=[];K=[];
    for j=1:m
        x=G.u[j];
        y=G.v[j];
        # if (xyinT[x]==1) && (xyinT[y]==1)
        #     a[x,y]+=1;
        #     a[y,x]+=1;
        #     continue
        # end
        for i=1:p
            w1=Int64[];
            w2=Int64[];
            xx=x;yy=y;
            lx=1;ly=1;
            push!(w1,xx);push!(w2,yy);
            # while (xyinT[xx]==0) && (lx<len)
            #     # nx=G.nbr[xx][Int(floor(rand()*(size(G.nbr[xx])[1])))+1];
            #     nx = rand(G.nbr[xx])
            #     lx+=1;
            #     push!(w1,nx);
            #     xx=nx;
            #     # if lx==len
            #     #     break
            #     # end
            # end
            for ll in 1:len
                if xyinT[xx]==1
                    lx=ll;
                    break
                end
                nx = rand(G.nbr[xx])
                push!(w1,nx);
                xx=nx;
            end
            # while (xyinT[yy]==0) && (ly<len)
            #     # ny=G.nbr[yy][Int(floor(rand()*(size(G.nbr[yy])[1])))+1];
            #     ny = rand(G.nbr[yy])
            #     ly+=1;
            #     push!(w2,ny);
            #     yy=ny;
            #     # if ly==len
            #     #     break
            #     # end
            # end
            for ll in 1:len
                if xyinT[yy]==1
                    ly=ll;
                    break
                end
                ny = rand(G.nbr[yy])
                push!(w2,ny);
                yy=ny;
            end
            ii=w1[lx];jj=w2[ly];
            # if (xyinT[ii]==1) && (xyinT[jj]==1) && (ii != jj)
            #     # push!(I,ii,jj)
            #     # push!(J,jj,ii)
            #     # push!(K,1/(p*(lx+ly-1)),1/(p*(lx+ly-1)))
            #     a[ii,jj]+=1/(p*(lx+ly-1));
            #     a[jj,ii]+=1/(p*(lx+ly-1));
            #     # wlkedge+=(lx+ly-1)
            # end

            ###########buyong a
            # a[ii,jj]+=xyinT[ii]*xyinT[jj]/(p*(lx+ly-1));
            # a[jj,ii]+=xyinT[ii]*xyinT[jj]/(p*(lx+ly-1));


            # if (xyinT[ii]==1) && (xyinT[jj]==1)
            #     wuse+=1;
            #     wlkedge+=(lx+ly-1)
            # end
            # a[ii,jj]+=((xyinT[ii]==1) && (xyinT[jj]==1) && (ii != jj))*1/(p*(lx+ly-1));
            # a[jj,ii]+=((xyinT[ii]==1) && (xyinT[jj]==1) && (ii != jj))*1/(p*(lx+ly-1));
            # # if (xyinT[ii]==1) && (xyinT[jj]==1)
            #     wuse+=1;
            # end
            # if (xyinT[ii]==1) || (xyinT[jj]==1)
            #     wend+=1;
            #     wleng+=lx+ly-1
            # end
            # wlks+=1;
            # wij=zeros(lx+ly);
            # wij[1:lx]=w1[lx:-1:1];wij[lx+1:lx+ly]=w2[1:ly];
            #wij=vcat(w2[ly:-1:1],w1)
            wlks+=1;
            if xyinT[ii]*xyinT[jj]==1
                wuse+=1;
                reverse!(w2)
                append!(w2,w1)
                push!(W,w2);
                push!(We,ly);
                wall+=lx+ly-1;
            end
        end
    end
    # a=sparse(I,J,K)
    # a[T,T].-=1;
    # println("wlks",' ',wlks,' ',wuse,' ',wuse/wlks)
    # return a[T,T],W,We,wuse/wlks,wall/wuse;
    return W,We,wuse/wlks,wall/wuse;
end


function init_direct(G,T,p,len)
    n=G.n;
    m=G.m;
    t=size(T)[1];
    #beta=t/m;

    a=spzeros(n,n);
    wlks=0;
    wlkedge=0;
    wuse=0;
    wgood=0;
    W=[];
    We=[];
    xyinT=zeros(Int,n)
    for i in T
        xyinT[T]=1;
    end
    for j=1:m
        x=G.u[j];
        y=G.v[j];
        if (xyinT[x]==1) && (xyinT[y]==1)
            a[x,y]+=1;
            #a[y,x]+=1;
            continue
        end
        for i=1:p
            w1=zeros(len);
            w2=zeros(len);
            xx=x;yy=y;
            lx=1;ly=1;
            w1[lx]=xx;w2[ly]=yy;
            while !(xyinT[xx]==1)
                nx=G.nbr_in[xx][Int(floor(rand()*(size(G.nbr_in[xx])[1])))+1];
                lx+=1;
                w1[lx]=nx;
                xx=nx;
                if lx==len
                    break
                end
            end
            while !(xyinT[yy]==1)
                ny=G.nbr_out[yy][Int(floor(rand()*(size(G.nbr_out[yy])[1])))+1];
                ly+=1;
                w2[ly]=ny;
                yy=ny;
                if ly==len
                    break
                end
            end
            ii=Int(w1[lx]);jj=Int(w2[ly]);
            if (xyinT[ii]==1) && (xyinT[jj]==1) && (ii != jj)
                a[ii,jj]+=1/(p*(lx+ly-1));
                #a[jj,ii]+=1/(p*(lx+ly-1));
                wlkedge+=(lx+ly-1)
                wuse+=1;
            end

            wlks+=1;
            # wij=zeros(lx+ly);
            # wij[1:lx]=w1[lx:-1:1];wij[lx+1:lx+ly]=w2[1:ly];
            wij=vcat(w1,w2[ly:-1:1])
            push!(W,wij);
            push!(We,lx);
        end
    end
    #println("wlks",' ',wlks,' ',wlkedge,' ',wuse,' ',wuse/wlks)
    return a[T,T],W,We;
end


function addTerminal(S0,S1,G,W,We,p,h)
    n=G.n;
    m=G.m;
    S=union(S0,S1);
    s0=size(S0)[1];s1=size(S1)[1];s=s0+s1;
    V=union(1:n);
    F=setdiff(V,S);
    orig=zeros(Int,s);  ## orig[1..s]  -> 1..n
    label=zeros(Int,n); ## label[1..n] -> 1..s

    tt1=time();
    for i=1:s
        orig[i]=S[i];
        label[S[i]]=i;
    end
    # change_a=[spzeros(n,n) for i=1:n];
    change_a=zeros(n,s,s);
    change_new=zeros(n,s);
    w=size(W)[1];
    vis_l=zeros(Int,n);
    vis_r=zeros(Int,n);
    first_appear_l=ones(1000);
    first_appear_r=ones(1000);
    tt2=time()
    einS=zeros(Int,n);
    einS[S].=1;
    # pl=1/(p*len);

    for i=1:w
        lenn=size(W[i])[1];
        pl=1/(p*(lenn-1));
        e_s=W[i][1];e_t=W[i][lenn];
        # println(e_s,' ',e_t)
        vis_l[W[i][1:We[i]]]=1:We[i];
        vis_r[W[i][lenn:-1:(We[i]+1)]]=lenn:-1:(We[i]+1);
        # first_appear_l[1:We[i]].=(vis_l[W[i][1:We[i]]].==1:We[i]);
        # first_appear_r[We[i]+1:lenn].=(vis_r[W[i][We[i]+1:lenn]].==(We[i]+1:lenn));
        change_a[W[i][2:We[i]],label[e_s],label[e_t]].-=pl;
        change_a[W[i][2:We[i]],label[e_t],label[e_s]].-=pl;
        change_a[W[i][We[i]+1:lenn-1],label[e_s],label[e_t]].-=pl;
        change_a[W[i][We[i]+1:lenn-1],label[e_t],label[e_s]].-=pl;
        # change_new[W[i][2:We[i]],label[e_t]]-=first_appear_l[2:We[i]].*(vis_r[W[i][2:We[i]]].==0)./(p*(2-lenn:We[i]-lenn));
        # change_new[W[i][We[i]+1:lenn-1],label[e_s]]+=first_appear_r[We[i]+1:lenn-1].*(vis_l[W[i][We[i]+1:lenn-1]].==0)./(p*(We[i]:lenn-2));
        change_new[W[i][2:We[i]],label[e_t]]-=(vis_l[W[i][2:We[i]]].==2:We[i]).*(vis_r[W[i][2:We[i]]].==0)./(p*(2-lenn:We[i]-lenn));
        change_new[W[i][We[i]+1:lenn-1],label[e_s]]+=(vis_r[W[i][We[i]+1:lenn-1]].==(We[i]+1:lenn-1)).*(vis_l[W[i][We[i]+1:lenn-1]].==0)./(p*(We[i]:lenn-2));

        vis_l[W[i][1:We[i]]].=0;
        vis_r[W[i][(We[i]+1):lenn]].=0;
    end
    #
    # for i=1:w
    #      # ttt1=time()
    #     len=size(W[i])[1];
    #     pl=1/(p*(len-1));
    #     e_s=W[i][1];e_t=W[i][len];
    #      ### e_s,e_t : 1..n
    #     # if (einS[e_s]==0) && (einS[e_t]==0)
    #     #     continue
    #     # end
    #     # for j=1:We[i]
    #     #     vis_l[W[i][j]]=j;
    #     # end
    #     # ttt2=time()
    #     vis_l[W[i][1:We[i]]]=1:We[i];
    #     first_appear_l[1:We[i]].=(vis_l[W[i][1:We[i]]].==1:We[i]);
    #     # for j=len:-1:(We[i]+1)
    #     #     vis_r[W[i][j]]=j;
    #     # end
    #     vis_r[W[i][len:-1:(We[i]+1)]]=len:-1:(We[i]+1);
    #     first_appear_r[len:-1:(We[i]+1)].=(vis_r[W[i][len:-1:(We[i]+1)]].==len:-1:(We[i]+1));
    #
    #     # ttt2=time()
    #     # for j=1:We[i]
    #     #     if  (j==vis_l[Int(W[i][j])]) && (vis_r[Int(W[i][j])]==0)
    #     #         change_new[Int(W[i][j]),Int(label[e_t])]+=1/(p*(len-j));
    #     #         # if einS[e_s]==1
    #     #         change_a[Int(W[i][j])][Int(label[e_s]),Int(label[e_t])]-=pl;
    #     #         change_a[Int(W[i][j])][Int(label[e_t]),Int(label[e_s])]-=pl;
    #     #         # end
    #     #     end
    #     # end
    #     # ttt3=time()
    #     ############不写循环
    #     change_a[W[i][2:We[i]],label[e_s],label[e_t]].-=pl;
    #     change_a[W[i][2:We[i]],label[e_t],label[e_s]].-=pl;
    #     change_new[W[i][2:We[i]],label[e_t]]-=first_appear_l[2:We[i]].*(vis_r[W[i][2:We[i]]].==0)./(((2:We[i]).-len).*p)
    #
    #
    #     #############end不写循环
    #
    #
    #     # ttt3=time()
    #     # for j=We[i]+1:len
    #     #     if (j==vis_r[Int(W[i][j])]) && (vis_l[Int(W[i][j])]==0)
    #     #         change_new[Int(W[i][j]),Int(label[e_s])]+=1/(p*(j-1));
    #     #         # if einS[e_t]==1
    #     #         change_a[Int(W[i][j])][Int(label[e_s]),Int(label[e_t])]-=pl;
    #     #         change_a[Int(W[i][j])][Int(label[e_t]),Int(label[e_s])]-=pl;
    #     #         # end
    #     #     end
    #     # end
    #
    #     ##########
    #     change_a[W[i][len-1:-1:We[i]+1],label[e_s],label[e_t]].-=pl;
    #     change_a[W[i][len-1:-1:We[i]+1],label[e_t],label[e_s]].-=pl;
    #     change_new[W[i][len-1:-1:We[i]+1],label[e_s]]-=first_appear_l[len-1:-1:We[i]+1].*(vis_l[W[i][len-1:-1:We[i]+1]].==0)./(((len-1:-1:We[i]+1).-len).*p)
    #     #########
    #
    #
    #      # ttt4=time()
    #     vis_l[W[i][1:We[i]]].=0
    #     vis_r[W[i][(We[i]+1):len]].=0
    #     # for j=1:len
    #     #     vis_l[Int(W[i][j])]=0;
    #     #     vis_r[Int(W[i][j])]=0;
    #     # end
    #     #  ttt5=time()
    #     # if i==1133
    #     #     println("intime:")
    #     #     println("intime:")
    #     #     println(ttt2-ttt1,' ',ttt3-ttt2,' ',ttt4-ttt3,' ',ttt5-ttt4)
    #     # end
    # end


    tt3=time();
    ans=zeros(n);
    x=ones(Int,n);
    x[S0].=0;
    maxx=0;selc=0;
    thetatmp=spzeros(n);
    atmp=zeros(s+1,s+1);
    JJJ=zeros(s+1);
    for i=1:s+1
        JJJ[i]=(i-1)*(s+1)+i;
    end
    for i in F
        # atmp[1:s,1:s]=a+change_a[i,1:s,1:s];
        atmp[1:s,1:s]=change_a[i,1:s,1:s];
        vect=change_new[i,1:s];
        #println()
        #println(vect)
        atmp[s+1,1:s]=vect;
        atmp[1:s,s+1]=vect;
        for j=1:s+1
            atmp[j,j]=0;
            atmp[j,j]=-sum(atmp[j,:])
        end
        atmp=-atmp;
        # S1tmp=union(S1,i);

        Stmp=[S;i]
        thetatmp[Stmp]=atmp*x[Stmp];
        opitmp=h[Stmp]'*thetatmp[Stmp];
        thetatmp[i]=0;
        if opitmp>maxx
            maxx=opitmp
            selc=i
        end
    end
    tt4=time();
    lft=0;rht=0;
    for i=1:w

        len=size(W[i])[1];
        lft=1;rht=len;
        for j=We[i]:-1:1
            if W[i][j]==selc
                lft=j;
                break;
            end
        end
        for j=We[i]+1:len
            if W[i][j]==selc
                rht=j
                break
            end
        end
        # if lft+rht>0
        #     if lft==0
        #         lft=1;
        #     end
        #     if rht==0
        #         rht=len;
        #     end
        if lft!=1 || rht!=len
            We[i]=We[i]-lft+1;
            W[i]=W[i][lft:rht];
        end
    end
    tt5=time()
    # i=selc;
    # atmp=zeros(s+1,s+1);
    # atmp[1:s,1:s]=a+change_a[i,:,:];
    # vect=change_new[i,:];
    # #println()
    # #println(vect)
    # atmp[s+1,1:s]=vect;
    # atmp[1:s,s+1]=vect;
    # for j=1:s+1
    #     atmp[j,j]=0;
    # end
    println("fun time:")
    println(tt2-tt1,' ',tt3-tt2,' ',tt4-tt3,' ',tt5-tt4)
    return selc,W,We;
end
