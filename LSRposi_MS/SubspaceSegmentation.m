function acc = SubspaceSegmentation( SegmentatiomMethod , X , gnd , Par )

switch SegmentatiomMethod
    
    case 'LSR1'
        C = LSR1( X , Par ) ;
        
    case 'LSR2'
        C = LSR2( X , Par ) ;
        
    case 'LSRd0po'
        C = LSRd0po( X , Par ) ;
        
    case 'LSRpo'
        C = LSRpo( X , Par ) ;
        
end

nCluster = length( unique( gnd ) ) ;
Z = ( abs(C) + abs(C') ) / 2 ;
idx = clu_ncut(Z,nCluster) ;
acc = compacc(idx,gnd) ;
