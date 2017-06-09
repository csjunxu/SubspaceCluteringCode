function acc = SubspaceSegmentation( SegmentatiomMethod , X , gnd , Par )

switch SegmentatiomMethod
        
    case 'LSRd0po'
        C = LSRd0po( X , Par ) ;
        
    case 'LSRpo'
        C = LSRpo( X , Par ) ;
        
end

for i = 1 : size(C,2)
   C(:,i) = C(:,i) / max(abs(C(:,i))) ;    
end

nCluster = length( unique( gnd ) ) ;
Z = ( abs(C) + abs(C') ) / 2 ;
idx = clu_ncut(Z,nCluster) ;
acc = compacc(idx,gnd) ;

