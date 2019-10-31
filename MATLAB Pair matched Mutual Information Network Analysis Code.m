clc
close all 
clear all
tic 
B = importdata ('BDS_3MS_PMMELDNA14.csv'); % Survivors
A = importdata ('BDS_3MNS_PMMELDNA14.csv'); % Non-survivors, labelled as A because the pair-matched matrix arrays will never be bigger than the smallest dataset
BB = importdata ('BDS_6MS_PMMELDNA14.csv'); 
AA = importdata ('BDS_6MNS_PMMELDNA14.csv');
BBB = importdata ('BDS_12MS_PMMELDNA14.csv'); 
AAA = importdata ('BDS_12MNS_PMMELDNA14.csv');
toc
tmpA = size(A);
C = NaN(tmpA); 
D = NaN(tmpA);
for i = 1:size(A,1)
    tmp = B(B(:,6)==A(i,6),:);
     C(i,:) = A(i,:);
    if isempty(tmp) % No matches
        tmp = B(B(:,6)<=A(i,6)+0.5 & B(:,6)>=A(i,6)-0.5,:); % Expand range of check
    end
    if size(tmp,1)==1 & ~isempty(tmp)
        D(i,:) = tmp;
    elseif size(tmp,1)>1
        r = randi([1,size(tmp,1)]);
        D(i,:) = tmp(r,:);
    end
end
del=find(all(isnan(D),2));
for i = size(del,1):-1:1;
    rmv = del(i,:);
    D(rmv,:) = [];
    C(rmv,:) = [];
end 
D(:,6) = [];
C(:,6) = [];
Surv3M = D;
NonSurv3M = C;
tmpAA = size(AA);
CC = NaN(tmpAA); 
DD = NaN(tmpAA);
for i = 1:size(AA,1)
    tmp = BB(BB(:,6)==AA(i,6),:);
     CC(i,:) = AA(i,:);
    if isempty(tmp) % No matches
        tmp = BB(BB(:,6)<=AA(i,6)+0.5 & BB(:,6)>=AA(i,6)-0.5,:); % Expand range of check
    end
    if size(tmp,1)==1 & ~isempty(tmp)
        DD(i,:) = tmp;
    elseif size(tmp,1)>1
        r = randi([1,size(tmp,1)]);
        DD(i,:) = tmp(r,:);
    end
end
del=find(all(isnan(DD),2));
for i = size(del,1):-1:1;
    rmv = del(i,:);
    DD(rmv,:) = [];
    CC(rmv,:) = [];
end 
DD(:,6) = [];
CC(:,6) = [];
Surv6M = DD;
NonSurv6M = CC;
tmpAAA = size(AAA);
CCC = NaN(tmpAAA); 
DDD = NaN(tmpAAA);
for i = 1:size(AAA,1)
    tmp = BBB(BBB(:,6)==AAA(i,6),:);
     CCC(i,:) = AAA(i,:);
    if isempty(tmp) % No matches
        tmp = BBB(BBB(:,6)<=AAA(i,6)+0.5 & BBB(:,6)>=AAA(i,6)-0.5,:); % Expand range of check
    end
    if size(tmp,1)==1 & ~isempty(tmp)
        DDD(i,:) = tmp;
    elseif size(tmp,1)>1
        r = randi([1,size(tmp,1)]);
        DDD(i,:) = tmp(r,:);
    end
end
del=find(all(isnan(DDD),2));
for i = size(del,1):-1:1;
    rmv = del(i,:);
    DDD(rmv,:) = [];
    CCC(rmv,:) = [];
end 
DDD(:,6) = [];
CCC(:,6) = [];
Surv12M = DDD;
NonSurv12M = CCC;
default = size(C);
nn = default(1,2)
A =  Surv3M;
A = A';
B = NonSurv3M;
B = B';
C = Surv6M;
C = C';
D = NonSurv6M;
D = D';
E = Surv12M;
E = E';
F = NonSurv12M;
F = F';
toc
th = 0.75; % threshold significance 
default = size(A,1);
X1 = zeros(default);
X2 = zeros(default);
X3 = zeros(default);
X4 = zeros(default);
P1B = zeros(default);
P2B = zeros(default);
P3B = zeros(default);
P4B = zeros(default);
P4B = zeros(default);
P6B = zeros(default);
for j=1:nn
   for i=j:(nn-1)
      X1(i+1,j)=kernelmi(A(j,:),A(i+1,:)); 
   end
end
for p=1:nn
    for q=p:nn
        X2(p,q)=X1(q,p);
    end
end
X=X1+X2;
for c=1:nn;
    for d=1:nn;
        if X(c,d)>=th;
            P1B(c,d)=X(c,d);
        else P1B(c,d)=0;
        end
    end
end
P1Bgraph=graph(P1B);
P1BAbsgraph=abs(P1B); 
Weight=P1BAbsgraph;
node_names = {'HE','PSShunt','Ascites','Diabetes','Pugh','Alb','Tot Bili','PT_pC','Creatinine','INR','Ammonia','Na','Hb','CRP'};
P1BRgraph = graph(P1BAbsgraph,node_names);
LWidths=5*P1BRgraph.Edges.Weight/max(P1BRgraph.Edges.Weight);
figure(1)
subplot (1,2,1)
plot(P1BRgraph,'LineWidth',LWidths,'Layout','force','UseGravity',true,'NodeColor','green','EdgeColor','black','MarkerSize',4)
% "'EdgeLabel',P2BRgraph.Edges.Weight," - option to label edges with weight in the plot function code
% Degree is the number of edges connecting to each node
% Betweenness measures how often each graph node appears on a shortest path between two nodes in the graph
% Closeness is the inverse sum of the distance from a node to all other nodes in the graph
% Eigenvector is the relative score assigned to the value of a node 
% Shortest paths is the shortest path possible between any 2 nodes
title ('Survivor 3M')
P1BR_degree = centrality(P1BRgraph,'Degree'); 
P1BRgraph.Nodes.degree = P1BR_degree;
P1BR_betweenness = centrality(P1BRgraph,'Betweenness'); 
P1BRgraph.Nodes.betweenness = P1BR_betweenness;
P1BR_closeness = centrality(P1BRgraph,'closeness'); 
P1BRgraph.Nodes.closeness = P1BR_closeness;
P1BR_eigenvector = centrality(P1BRgraph,'eigenvector'); 
P1BRgraph.Nodes.eigenvector = P1BR_eigenvector;
output1 = P1BRgraph.Nodes;
P1BRshortestpaths=distances(P1BRgraph)
for e=1:nn
   for f=e:(nn-1)
      X3(f+1,e)=kernelmi(B(e,:),B(f+1,:));
   end
end
for r=1:nn
    for s=r:nn
        X4(r,s)=X3(s,r);
    end
end
XT=X3+X4;
for m=1:nn;
    for n=1:nn;
        if XT(m,n)>=th;
            P2B(m,n)=XT(m,n);
        else P2B(m,n)=0;
        end
    end
end
P2Bgraph=graph(P2B)
P2BAbsgraph=abs(P2B);
Weight=P2BAbsgraph;
node_names = {'HE','PSShunt','Ascites','Diabetes','Pugh','Alb','Tot Bili','PT_pC','Creatinine','INR','Ammonia','Na','Hb','CRP'};
P2BRgraph = graph(P2BAbsgraph,node_names);
LWidths=5*P2BRgraph.Edges.Weight/max(P2BRgraph.Edges.Weight);
figure(1)
subplot (1,2,2)
plot(P2BRgraph,'LineWidth',LWidths,'Layout','force','UseGravity',true,'NodeColor','red','EdgeColor','black','MarkerSize',4)
title ('Non survivor 3M')
P2BR_degree = centrality(P2BRgraph,'Degree'); 
P2BRgraph.Nodes.degree = P2BR_degree;
P2BR_betweenness = centrality(P2BRgraph,'Betweenness'); 
P2BRgraph.Nodes.betweenness = P2BR_betweenness;
P2BR_closeness = centrality(P2BRgraph,'closeness'); 
P2BRgraph.Nodes.closeness = P2BR_closeness;
P2BR_eigenvector = centrality(P2BRgraph,'eigenvector'); 
P2BRgraph.Nodes.eigenvector = P2BR_eigenvector;
output2 = P2BRgraph.Nodes
P2BRshortestpaths=distances(P2BRgraph) 
for j=1:nn
   for i=j:(nn-1)
      X1(i+1,j)=kernelmi(C(j,:),C(i+1,:)); 
   end
end
for p=1:nn
    for q=p:nn
        X2(p,q)=X1(q,p);
    end
end
X=X1+X2;
for c=1:nn;
    for d=1:nn;
        if X(c,d)>=th;
            P3B(c,d)=X(c,d);
        else P3B(c,d)=0;
        end
    end
end
P3Bgraph=graph(P3B)
P3BAbsgraph=abs(P3B); 
Weight=P3BAbsgraph;
node_names = {'HE','PSShunt','Ascites','Diabetes','Pugh','Alb','Tot Bili','PT_pC','Creatinine','INR','Ammonia','Na','Hb','CRP'};
P3BRgraph = graph(P3BAbsgraph,node_names);
LWidths=5*P3BRgraph.Edges.Weight/max(P3BRgraph.Edges.Weight);
figure(2)
subplot (1,2,1)
plot(P3BRgraph,'LineWidth',LWidths,'Layout','force','UseGravity',true,'NodeColor','green','EdgeColor','black','MarkerSize',4)
title ('Survivor 6M')
P3BR_degree = centrality(P3BRgraph,'Degree'); 
P3BRgraph.Nodes.degree = P3BR_degree;
P3BR_betweenness = centrality(P3BRgraph,'Betweenness'); 
P3BRgraph.Nodes.betweenness = P3BR_betweenness;
P3BR_closeness = centrality(P3BRgraph,'closeness');
P3BRgraph.Nodes.closeness = P3BR_closeness;
P3BR_eigenvector = centrality(P3BRgraph,'eigenvector');
P3BRgraph.Nodes.eigenvector = P3BR_eigenvector;
output3 = P3BRgraph.Nodes;
P3BRshortestpaths=distances(P3BRgraph)
for e=1:nn
   for f=e:(nn-1)
      X3(f+1,e)=kernelmi(D(e,:),D(f+1,:));
   end
end
for r=1:nn
    for s=r:nn
        X4(r,s)=X3(s,r);
    end
end
XT=X3+X4;
for m=1:nn;
    for n=1:nn;
        if XT(m,n)>=th;
            P4B(m,n)=XT(m,n);
        else P4B(m,n)=0;
        end
    end
end
P4Bgraph=graph(P4B)
P4BAbsgraph=abs(P4B);
Weight=P4BAbsgraph;
node_names = {'HE','PSShunt','Ascites','Diabetes','Pugh','Alb','Tot Bili','PT_pC','Creatinine','INR','Ammonia','Na','Hb','CRP'};
P4BRgraph = graph(P4BAbsgraph,node_names);
LWidths=5*P4BRgraph.Edges.Weight/max(P4BRgraph.Edges.Weight);
figure(2)
subplot (1,2,2)
plot(P4BRgraph,'LineWidth',LWidths,'Layout','force','UseGravity',true,'NodeColor','red','EdgeColor','black','MarkerSize',4)
title ('Non survivor 6M')
P4BR_degree = centrality(P4BRgraph,'Degree'); 
P4BRgraph.Nodes.degree = P4BR_degree;
P4BR_betweenness = centrality(P4BRgraph,'Betweenness'); 
P4BRgraph.Nodes.betweenness = P4BR_betweenness;
P4BR_closeness = centrality(P4BRgraph,'closeness');
P4BRgraph.Nodes.closeness = P4BR_closeness;
P4BR_eigenvector = centrality(P4BRgraph,'eigenvector');
P4BRgraph.Nodes.eigenvector = P4BR_eigenvector;
output4 = P4BRgraph.Nodes;
P4BRshortestpaths=distances(P4BRgraph)
for j=1:nn
   for i=j:(nn-1)
      X1(i+1,j)=kernelmi(E(j,:),E(i+1,:)); 
   end
end
for p=1:nn
    for q=p:nn
        X2(p,q)=X1(q,p);
    end
end
X=X1+X2;
for c=1:nn;
    for d=1:nn;
        if X(c,d)>=th;
            P5B(c,d)=X(c,d);
        else P5B(c,d)=0;
        end
    end
end
P5Bgraph=graph(P5B)
P5BAbsgraph=abs(P5B); 
Weight=P5BAbsgraph;
node_names = {'HE','PSShunt','Ascites','Diabetes','Pugh','Alb','Tot Bili','PT_pC','Creatinine','INR','Ammonia','Na','Hb','CRP'};
P5BRgraph = graph(P5BAbsgraph,node_names);
LWidths=5*P5BRgraph.Edges.Weight/max(P5BRgraph.Edges.Weight);
figure(3)
subplot (1,2,1)
plot(P5BRgraph,'LineWidth',LWidths,'Layout','force','UseGravity',true,'NodeColor','green','EdgeColor','black','MarkerSize',4)
title ('Survivor 12M')
P5BRgraph.Nodes;
P5BR_degree = centrality(P5BRgraph,'Degree'); 
P5BRgraph.Nodes.degree = P5BR_degree;
P5BR_betweenness = centrality(P5BRgraph,'Betweenness'); 
P5BRgraph.Nodes.betweenness = P5BR_betweenness;
P5BR_closeness = centrality(P5BRgraph,'closeness');
P5BRgraph.Nodes.closeness = P5BR_closeness;
P5BR_eigenvector = centrality(P5BRgraph,'eigenvector');
P5BRgraph.Nodes.eigenvector = P5BR_eigenvector;
output5 = P5BRgraph.Nodes
P5BRshortestpaths=distances(P5BRgraph)
for e=1:nn
   for f=e:(nn-1)
      X3(f+1,e)=kernelmi(F(e,:),F(f+1,:));
   end
end
for r=1:nn
    for s=r:nn
        X4(r,s)=X3(s,r);
    end
end
XT=X3+X4;
for m=1:nn;
    for n=1:nn;
        if XT(m,n)>=th;
            P6B(m,n)=XT(m,n);
        else P6B(m,n)=0;
        end
    end
end
P6Bgraph=graph(P6B)
P6BAbsgraph=abs(P6B);
Weight=P6BAbsgraph;
node_names = {'HE','PSShunt','Ascites','Diabetes','Pugh','Alb','Tot Bili','PT_pC','Creatinine','INR','Ammonia','Na','Hb','CRP'};
P6BRgraph = graph(P6BAbsgraph,node_names);
LWidths=5*P6BRgraph.Edges.Weight/max(P6BRgraph.Edges.Weight);
figure(3)
subplot (1,2,2)
plot(P6BRgraph,'LineWidth',LWidths,'Layout','force','UseGravity',true,'NodeColor','red','EdgeColor','black','MarkerSize',4)
title ('Non survivor 12M')
P6BR_degree = centrality(P6BRgraph,'Degree'); 
P6BRgraph.Nodes.degree = P6BR_degree;
P6BR_betweenness = centrality(P6BRgraph,'Betweenness'); 
P6BRgraph.Nodes.betweenness = P6BR_betweenness;
P6BR_closeness = centrality(P6BRgraph,'closeness');
P6BRgraph.Nodes.closeness = P6BR_closeness;
P6BR_eigenvector = centrality(P6BRgraph,'eigenvector');
P6BRgraph.Nodes.eigenvector = P6BR_eigenvector;
output6 = P6BRgraph.Nodes;
P6BRshortestpaths=distances(P6BRgraph)



