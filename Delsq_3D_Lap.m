function D = Delsq_3D_Lap(G)
[m,~] = size(G);

p = find(G);
i = G(p);
j = G(p);
s = -6*ones(size(p));

   k = -1; % Up
   Q = G(p+k);
   q = find(Q);
   i = [i; G(p(q))]; 
   j = [j; Q(q)];   
   s = [s; ones(length(q),1)];  

   k = m; % Down
   Q = G(p+k); 
   q = find(Q);
   i = [i; G(p(q))];
   j = [j; Q(q)];  
   s = [s; ones(length(q),1)]; 

   k = 1; % Right
   Q = G(p+k);
   q = find(Q);
   i = [i; G(p(q))]; 
   j = [j; Q(q)];  
   s = [s; ones(length(q),1)]; 

   k = -m; % Left
   Q = G(p+k);
   q = find(Q); 
   i = [i; G(p(q))];  
   j = [j; Q(q)];  
   s = [s; ones(length(q),1)];  

   k = -m*m; % Front
   Q = G(p+k);
   q = find(Q); 
   i = [i; G(p(q))];  
   j = [j; Q(q)];  
   s = [s; ones(length(q),1)];  

   k = m*m; % Back
   Q = G(p+k);
   q = find(Q); 
   i = [i; G(p(q))];  
   j = [j; Q(q)];  
   s = [s; ones(length(q),1)];  
   
D = sparse(i,j,s);

end

