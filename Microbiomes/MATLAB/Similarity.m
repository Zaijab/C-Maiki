M=csvread('CMAIKI.FAke.DATA.csv',1,1);
%M( ~any(M,2), : ) = [];  %rows
%M( :, ~any(M,1) ) = [];  %columns
b=sum(M.',2);
c=var(b);
[~,B]=size(M);
j=1;
for i=1:B
    if b(i)>c
        a(j)=i;
        j=j+1;
    elseif b(i)<=0.3
        a(j)=i;
        j=j+1;
    end
end
[~,C]=size(a);
for i=1:C
    d=a(C-i+1);
    M(:,d)=[];
    b(d)=[];
end
[~,B]=size(M);
N=(M.'./b).';
Dist1=zeros(B,B);
Dist2=zeros(B,B);
for i=1:B
    for j=i:B
        Dist1(i,j)=norm(N(:,j)-N(:,i),1);
        Dist2(i,j)=norm(N(:,j)-N(:,i));
    end
end

name=1:1:B;
x=1;
y=0.1;
z=1.4;
[s, t]=size(x:y:z);
Weights1=zeros(B,B,t);
cmap=flipud(copper(256));
f=12;
g=0:max(M(f,:))/10:max(M(f,:));
for k=x:y:z
    epsilon=k;
    DGraph1=Dist1;
    for i=1:B
        for j=i+1:B             
            if Dist1(i,j)>=epsilon
                DGraph1(i,j)=0;
            end
        end
    end
    Weights1(:,:,s)=DGraph1;
    G1=graph(DGraph1,'upper');
    f1=figure('Name',['P=1, ',char(949),'=', num2str(k)],...
        'NumberTitle','off');
    %plot(G1,'EdgeLabel',G1.Edges.Weight)  
    p1=plot(G1);
    for i=1:B
        row=round((size(cmap,1)-1)*(M(f,i)/max(M(f,:)))+1);
        color=cmap(row,:);
        highlight(p1,i,'NodeColor',color)
    end
    colormap(cmap)
    p=colorbar;
    set(p,'YTickLabel',g);
    %layout(p1,'force','WeightEffect','direct')
    %view(3)
    %X = kamada_kawai_spring_layout(G1);
    %p1=gplot(G1,X);
    %h(s)=f1;
    %s=s+1;
end

x=0;
y=0.1;
z=.5;
[~, t]=size(x:y:z);
Weights2=zeros(B,B,t);
for k=x:y:z
    epsilon=k;
    DGraph2=Dist2;
    for i=1:B
        for j=i+1:B  
            if Dist2(i,j)>=epsilon
                DGraph2(i,j)=0;
            end
        end
    end
    Weights2(:,:,s)=DGraph2;
    G2=graph(DGraph2,'upper');
    %n2=neighbors(G2,2);
    %n53=neighbors(G2,53);
    %n187=neighbors(G2,187);
    f2=figure('Name',['P=2, ',char(949),'=', num2str(k)],...
        'NumberTitle','off');
    %plot(G2,'EdgeLabel',G2.Edges.Weight)
    p2=plot(G2);
    for i=1:B
        row=round((size(cmap,1)-1)*(M(f,i)/max(M(f,:)))+1);
        color=cmap(row,:);
        highlight(p2,i,'NodeColor',color)
    end
    colormap(cmap)
    p=colorbar;
    set(p,'YTickLabel',g);
    %layout(p2,'force','WeightEffect','Direct')
    %view(3)
    %X = kamada_kawai_spring_layout(G2);
    %p2=gplot(G2,X);
    %h(s)=f2;
    %s=s+1;
end    

%csvwrite('graphweights1.csv',Weights1)
%csvwrite('graphweights2.csv',Weights2)
%savefig(h,'IntensityColorGraphsWater.fig')
%close(h)
