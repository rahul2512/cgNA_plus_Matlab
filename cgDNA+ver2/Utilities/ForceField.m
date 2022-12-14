ValenceAngleList = [4 1 6 ; 1 6 2 ; 6 2 3 ; 2 3 4 ; 3 4 1 ; 10 3 2 ; ... 
    10 3 4 ; 3 4 5 ; 5 4 1 ; 9 5 4 ; 5 9 11 ; 3 10 12 ; 1 6 7 ; 2 6 7] ;

ValenceAngleData = [ 2 12 9 6 8 10 10 7 3 4 5 11 1 13 ] ;

DihedralAngleList = [ 11 9 5 4 ; 9 5 4 3 ; 5 4 3 10 ; 4 3 10 12 ; ...
    1 6 7 8 ; 4 1 6 2 ; 1 6 2 3 ; 6 2 3 4 ; 2 3 4 1 ; 3 4 1 6 ; ...
    9 5 4 1 ; 2 3 10 12 ; 1 6 7 13 ; 2 6 7 8 ; 2 6 7 13 ] ;

DihedralAngleData = [ 1 6 2 8 4 5 3 9 10 7 12 13 11 6 6 ] ;

NonBondedList = [ 11 4 ; 9 3 ; 9 1 ; 5 10 ; 5 6 ; 5 2 ; 4 12 ; 6 10 ; 1 10 ];
NonBondedData = [ 4 3 ; 2 1 ; 2 2 ; 5 2 ; 5 3 ; 5 3 ;  3 4 ; 3 2 ; 3 4  ] ;

BondedList = [ 1 6 ; 6 2 ; 2 3 ; 3 10 ; 3 4 ; 4 5 ; 9 5 ; 1 4 ] ;
BondedData = [ 6 5 1 2 1 3 4 6  ] ;

BondedFF = {
'CE-CT' , 310.0 , 1.526 ;
'OS-CE' , 320.0 , 1.410 ;
'CI-CT' , 310.0 , 1.526 ;
'OS-CI' , 320.0 , 1.410 ;
'CT-CT' , 310.0 , 1.526 ;           
'CT-OS' , 320.0 , 1.410 ;      
};

NonBondedFF = {
'CE' , 1.9080 , 0.1094 ;
'OS' , 1.6837 , 0.1700 ;             
'CT' , 1.9080 , 0.1094 ;             
'P' , 2.1000 , 0.2000 ;           
'CI' , 1.9080 , 0.1094 ;
};

ChargesBackboneAtoms = {
'O4''' , -0.369100 ;
'C2''' , -0.085400 ; 
'C3''' ,  0.071300 ;
'C4''' ,  0.162900 ;
'C5''' , -0.006900 ;
'C1''' , [0.043100  0.068000 -0.011600 0.035800] ;
'N*'   ,  []       ;
'C*'   ,  []       ;
'O5''' , -0.495400 ;
'O3''' , -0.523200 ;   
'P'''  ,  1.165900 ;
'P'''  ,  1.165900 ;
}; 

ValenceAngleFF ={
'OS-CT-N*', 50.0 , 109.50 ;
'CT-OS-CT' , 60.0 , 109.50 ; 
'CI-CT-OS' , 50.0 , 109.50 ;
'OS-CI-CT' , 50.0 , 109.50 ;
'P -OS-CI' , 100.0 , 120.50 ;
'CT-CE-CT' , 40.0 , 109.50 ;
'CI-CT-CE' , 40.0 , 109.50 ;
'CE-CT-OS' , 50.0 , 109.50 ;
'CE-CT-CT' , 40.0 , 109.50 ; 
'OS-CE-CT' , 50.0 , 109.50 ; 
'P -OS-CE' , 100.0 , 120.50 ; 
'CT-CT-OS' , 50.0 , 109.50 ; 
'CT-CT-N*' , 50.0 , 109.50
} ;

% Data for 'OS-CT-N*-C*' 
% A=1,T=2,C=3,G=4

DihedralAngleFF= {
'X -CI-OS-X' , 3 , 1.150 , 0.0 , 3.0 ;
'X -CE-CT-X' , 9 , 1.400 , 0.0 , 3.0 ;
'X -CT-CT-X' , 9 , 1.40 , 0.0 , 3.0 ;      
'X -CT-N*-X' , 6 , 0.00 , 0.0 , 2.0 ;              
'CT-CT-OS-CT' , [1 1] , [0.383 0.1] , [0.0 180.0] , [3.0 2.0] ;
'CE-CT-CI-OS' , [1 1 1] , [1.178040 0.092102 0.962830] , [ 190.97653 295.63279 348.09535 ] , [1.0 2.0 3.0] ;  
'CE-CT-OS-CT' , [1 1] , [0.3833 0.100] , [0.0  180.0] , [3.0 2.0] ;     
'P -OS-CE-CT' , [1 1 1] , [1.079933 0.319017 0.054836] , [297.76503 200.61045 151.62658] , [1.0 2.0 3.0] ;	
'CT-CT-CE-CT' , [1 1] , [1.2 0.9] , [38.45 -2.23] , [1.0 3.0] ;           	
'CT-CE-CT-OS' , 1 , 1.2 , 322.41 , 1.0 ;
'OS-CT-N*-C*' , ones(4,3), [0.41 1.91 0.91 ; 0.17 1.67 0.99 ; 0.26 1.52 1.03; 0.44 1.73 0.81] , ... 
                [225.4 0.7 120.1 ; 210.9 4.7 188.0 ; 209.6 7.8 184.8 ; 219.8 4.0 87.5  ] , repmat([3.0 2.0 1.0], [4 1] );
'X -CI-CT-X'  , 9 , 1.400 , 0.0 , 3  ;
'P-OS-CE-CT'  , [1 1 1] , [1.079933 0.319017 0.054836] , [297.76503 200.61045 151.62658] , [1 2 3] 
} ;
