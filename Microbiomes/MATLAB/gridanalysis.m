function gridanalysis(a,g,A,r,w,t,random,P,metric,ep,p,q,edge,vertexc,heat)
% a is the list of points on the grid that were chosen
% g is the list of colors in hexadecimal format for the points chosen
% A is the list of columns to use a vectors
% r is the sites of the road
% P is the P norm to normalize the species across the sites.
% metric is the metric to use when evaluating 'closeness'
% ep is the epsilon to use to evaluate 'closeness'
% p is the power to raise the edge coordinates to.
% q is the root to take of the edge coordinates
% edge is used to tell if edges should be used, 0~No, 1~Yes
% vertexc is used to tell if the vertex should be colored
% heat is used to tell if a heatmap should be generated
%------------------------------------------------------------------------
   M=csvread('grid_coord.csv');                                            %Gets the coordinate data
    g=hex2dec(g)./hex2dec('ffffff');                                       %Converts colors into numbers
   [~, b]=size(r);
   [~, c]=size(a);
   Dist1=zeros(c,c);
   Dist2=zeros(c,b);
   for i=1:c                                                               %Creating a distance matrix between the chosen points
       for j=1:c
           Dist1(i,j)=norm(M(a(i),:)-M(a(j),:));
       end
   end
   for i=1:b                                                               %Creating a distance matrix between the chosen points and the road
       for j=1:c
            Dist2(j,i)=norm(M(a(j),:)-M(r(i),:));
       end
   end
   N=zeros(c,6);
   M(:,1)=M(:,1)-mean(M(:,1));                                             %Recentering the x-values so they do not dominate
   M(:,2)=M(:,2)-mean(M(:,2));                                             %Recentering the y-values so they do not dominate
   for i=1:c
       N(i,1)=M(a(i),1);                                                   %Listing of the x-values
       N(i,2)=M(a(i),2);                                                   %Listing of the y-values
       N(i,3)=g(i);                                                        %Listing of the colors
       N(i,4)=min(Dist2(i,:));                                             %Listing the min distance between the chosen points and the road.
       d=sort(Dist1(i,:));                                                 %Listing the min distance between the chosen points
       N(i,5)=d(2);
       N(i,6)=sqrt(sum(Dist1(i,:).^2));                                    %Listing total distance between each point
   end
   [d , ~]=size(A);
   I=csvread('CMAIKI.FAke.DATA.csv',1,1);
   [~, e]=size(I);
   J=zeros(c,e);
   random=sort(random);
   for i=1:c
       J(i,:)=I(random(i),:);
   end
   if heat~=0
       Heatmap(J)
   end
   for i=1:c
       J(i,:)=I(random(i),:)./(norm(I(random(i),:),P));
   end
   K=zeros(c,c);
   for i=1:c
       for j=i:c
           K(i,j)=norm(I(i,:)-I(j,:),metric);
           if K(i,j)>ep
               K(i,j)=0;
           end
       end
   end
   n=0:0.01:1;
   names=["Longitude, ","Latitude, ","Color, ","Distance To Road, ",...
       "Min Distance To Other Sites, ","Distance To Other Sites, "...
       "Longitude","Latitude","Color","Distance To Road",...
       "Min Distance To Other Sites","Distance To Other Sites"];
   labels=cellstr(num2str((1:c)'));
   cmap=flipud(copper(256));
   g=min(K(K>0)):(max(max(K))-min(K(K>0)))/10:max(max(K));
   vcolor={[.8 0 .8],[0 0 1],[0 1 0],[0.71 0.25 0.05]};    %Road, Water, Tree, Soil
   for i=1:d
      Y=zeros(c,3);
      for j=1:c
          m=sqrt(N(j,A(i,1))^2+N(j,A(i,2))^2+N(j,A(i,3))^2);               %Finding the norm of the 3 columns chosen vectors
          for k=1:3
              Y(j,k)=N(j,(A(i,k)))/m;                                      %Normalizing the 3 columns
          end
      end
      figure('Name',[num2str(names(A(i,1))),num2str(names(A(i,2))),...     %Creating the Figure
          num2str(names(A(i,3)+6))],'NumberTitle','off');
      xlabel(num2str(names(A(i,1)+6)));
      ylabel(num2str(names(A(i,2)+6)));
      zlabel(num2str(names(A(i,3)+6)));
      hold on
      [x, y, z] = sphere(128);                                             %Creating the Sphere
      h = surfl(x, y, z);
      shading interp
      if edge~=0
          for j=1:c
              for k=i:c
                  if K(j,k)~=0
                      B(:,1)=n*Y(j,1)+(1-n)*Y(k,1);
                      B(:,2)=n*Y(j,2)+(1-n)*Y(k,2);
                      B(:,3)=n*Y(j,3)+(1-n)*Y(k,3);
                      if p~=1
                          B=B./((sum((B.').^p).^(1/q)).');
                      end
                      row=round((size(cmap,1)-1)*((K(j,k)-min(K(K>0)...
                          ))/(max(max(K))-min(K(K>0))))+1);
                      color=cmap(row,:);
                      plot3(B(:,1),B(:,2),B(:,3),'Color',color,...
                          'LineWidth',4)
                  end
              end
          end
      end
      if edge~=0
          colormap(cmap)
          cthree=colorbar;
          set(cthree,'YTickLabel',g);
      end
      set(h, 'FaceAlpha', 0.2,'FaceColor',[.95 .95 .95],'FaceLighting',...
          'flat')
      light('Position',[1 1 1],'Style','local')
      s=40;
      if vertexc~=0
          for j=1:c
              if ismember(a(j),r)                                          %Road
                  plot3(Y(j,1),Y(j,2),Y(j,3),'.','MarkerSize',s,...
                      'Color',vcolor{1})
              elseif ismember(a(j),w)                                      %Water
                  plot3(Y(j,1),Y(j,2),Y(j,3),'.','MarkerSize',s,...
                      'Color',vcolor{2})
              elseif ismember(a(j),t)                                      %Tree
                  plot3(Y(j,1),Y(j,2),Y(j,3),'.','MarkerSize',s,...
                      'Color',vcolor{3})
              else
                  plot3(Y(j,1),Y(j,2),Y(j,3),'.','MarkerSize',s,...       %Soil
                      'Color',vcolor{4})
              end
          end
      else
          plot3(Y(:,1),Y(:,2),Y(:,3),'.','MarkerSize',s,'Color','b')
      end
      text(Y(:,1),Y(:,2),Y(:,3),labels,'VerticalAlignment','bottom',...
          'HorizontalAlignment','right')
      if vertexc~=0
          h = zeros(3, 1);
          h(1) = plot(NaN,NaN,'o','MarkerFaceColor',vcolor{1},...
              'MarkerEdgeColor',vcolor{1});
          h(2) = plot(NaN,NaN,'o','MarkerFaceColor',vcolor{2},...
              'MarkerEdgeColor',vcolor{2});
          h(3) = plot(NaN,NaN,'o','MarkerFaceColor',vcolor{3},...
              'MarkerEdgeColor',vcolor{3});
          h(4) = plot(NaN,NaN,'o','MarkerFaceColor',vcolor{4},...
              'MarkerEdgeColor',vcolor{4});
          legend(h,'Road', 'Water', 'Tree', 'Soil')
      end
      hold off
  end
end

%Example ----------------------------------------------------------

%Sites=[5 9 18 19 20 29 41 42 52 60 76 82 87 92 108 110 111 114 115 143];
%Color=["b7b020", "587611", "1b3909", "777b1f", "204c05", "4a8215",...
%    "2c540a","213d04", "4f880a", "333f0a", "41584b", "c6c105",...
%    "283818","2f3608", "b8bc04", "4c5342", "7b7e4f", "8daf09",...
%    "756a1e", "112e0a"];
%Vectors=[1 2 3; 1 2 4; 1 2 5; 1 2 6; 1 3 4; 1 3 5; 1 3 6; 1 4 5; 1 4 6;...
%       1 5 6; 2 3 4; 2 3 5; 2 3 6; 2 4 5; 2 4 6; 2 5 6; 3 4 5; 3 4 6;...
%       3 5 6; 4 5 6];
%Road=[85 86 97 98 110 111 122 123 135];
%Water=[27, 38, 39, 40, 50, 51, 62, 63, 75, 76, 87, 88, 100, 112, 124,...
%       125, 136, 137];
%Tree=[5:1:12, 17:1:24, 40, 41, 42, 47, 48, 49, 52, 53, 54, 59, 60, 61,...
%       64, 65, 66, 71, 72, 73, 74, 77,86, 87,91, 92, 93, 99, 101:1:105,...
%       113:1:117,126:1:129,138:1:144];
%random=[ 24, 25, 38, 8, 22, 6, 21, 23, 36, 1, 12, 35, 16, 15, 11,...
%       37, 28, 32, 2, 30];
%P=2;metric=2;epsilon=0.4;p=2;q=2;
%edge=1;vertexc=1;heat=1;%0~No, 1~Yes
%gridanalysis(Sites,Color,Vectors,Road,Water, Tree, random,P,metric,...
%epsilon,p,q,edge,vertexc,heat)
%-------------------------------------------------------------------
%Sites=[8 16 24 31 32 35 37 38 45 56 67 78 99 100 103 110 116 128 131 136];
%Color=["305565", "4a6519", "2c3210", "6d6326", "44453e", "353c16", "3b5a62",...
% "405b17", "597720", "99a437", "717d25", "82732b", "464338", "5a8d28",...
% "937f34", "384617", "7b7328", "3f551a", "666525" "69942c"];
%random=[38, 19, 31, 17, 12, 30, 36, 4, 20, 28, 6, 26, 22, 13, 34, 24,...
%32, 7, 3, 33];
%Vectors=[1 2 4; 1 3 5; 1 4 5];
%Road=[85 86 97 98 110 111 122 123 135];
%Water=[27, 38, 39, 40, 50, 51, 62, 63, 75, 76, 87, 88, 100, 112, 124,...
%       125, 136, 137];
%Tree=[5:1:12, 17:1:24, 40, 41, 42, 47, 48, 49, 52, 53, 54, 59, 60, 61,...
%       64, 65, 66, 71, 72, 73, 74, 77,86, 87,91, 92, 93, 99, 101:1:105,...
%       113:1:117,126:1:129,138:1:144];
%P=2;metric=2;epsilon=0.4;p=2;q=2; 
%edge=0;vertexc=1;heat=1;%0~No, 1~Yes
%gridanalysis(Sites,Color,Vectors,Road,Water, Tree, random,P,metric,...
%epsilon,p,q,edge,vertexc,heat)