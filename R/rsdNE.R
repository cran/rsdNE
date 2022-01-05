#'@name sym
#'@aliases sym
#'@title This generates a class of symmetric rotatable response surface designs with neighbour effects under a polynomial model of order (s1-1)
#'@examples
#' library(rsdNE)
#' sym(2,2,0.5)
#' ##output:
#' ## X matrix
#'#      [,1] [,2] [,3]
#'# [1,]    1   -1   -1
#'# [2,]    1    1    1
#'# [3,]    1    1   -1
#'# [4,]    1   -1    1
#'# [5,]    1   -1   -1
#'# [6,]    1    1    1
#'# [7,]    1   -1    1
#'# [8,]    1    1   -1
#'# [9,]    1   -1   -1
#'#[10,]    1    1    1
#'## Z prime Z matrix

#'#     [,1] [,2] [,3]
#'#[1,]   32    0    0
#'#[2,]    0    4    0
#'#[3,]    0    0    4
#'## Z prime Z inverse matrix
#'#     [,1]   [,2]   [,3]
#'#[1,] 0.03125 0.00 0.00
#'#[2,] 0.00000 0.25 0.00
#'#[3,] 0.00000 0.00 0.25
#'#[1] "total number of runs" "8"
#'#[1] "variance of esitmated response" "0.5312"
#'@references Sarika et al.2009, Communications in Statistics-Theory and Methods; Sarika et al.2013, Ars Combinatoria
#' @param s1 Number of levels of n1 factors, 1<s1<=6
#' @param n1 Number of factors, 1<n1<=4
#' @param c Value of alpha (Coefficient of neighbour effects), 0<=c<=1
#' @description This function generates symmetrical rotatable response surface designs in the presence of neighbour effects for n1 factors each at s1 levels.
#'@returns  his function generates rotatable designs as well as Z_prime_Z matrix,
#'inv(Z_primeZ) matrix and variance estimated response for the (s1^n1) factorial combination.
#'@author
#'Ashutosh Dalal, Division of Design of Experiments,ICAR-IASRI, New Delhi.
#'Seema Jaggi, Education Division, ICAR, Krishi Anusandhan Bhawan - II, Pusa, New Delhi.
#'Eldho Varghese,Fishery Resources Assessment Division,ICAR-CMFRI, Kochi.
#'Subhasish Sarkar, Division of Computer Application,ICAR-IASRI, New Delhi.
#'Arpan Bhowmik, Division of Design of Experiments,ICAR-IASRI, New Delhi.
#'Cini Varghese, Division of Design of Experiments,ICAR-IASRI, New Delhi.
#'Anindita Datta, Division of Design of Experiments,ICAR-IASRI, New Delhi.
#'Soumen Pal, Division of Computer Application,ICAR-IASRI, New Delhi.
#'@export
sym<-function(s1,n1,c){

  # s1 no. of levels corresponding to the n1 no. of factors
  s1
  n1

  # s2 no. of levels corresponding to the n2 no. of factors
  s2=1
  n2=0
  # s3 no. of levels corresponding to the n3 no. of factors
  s3=1
  n3=0
  N=(s1^n1)*(s2^n2)*(s3^n3)

  #alpha value lies between 0 to 1
  c

  k<-max(n1,n2,n3)


  #generation of columns of the design matrix
  if(s1>1){
    l1=as.integer(s1/2)
    v1<-seq(l1,-l1)
    if((s1)%%2==0){

      v1<-seq((s1-1),-(s1-1),by=-2)

    }
    v1<-t(v1)




    y11<-matrix(,nrow = (s1^n1)*(s2^n2)*(s3^n3), ncol = 0)
    j=1
    e=n1
    while(j<=n1 && e>0){
      x1<-matrix(,nrow = 1,ncol = 0)
      for (i in v1){
        if(n1>1){
          x2<-t(rep(c(i),times=(s1^(n1-j))))
        }else{
          x2<-t(rep(c(i),times=N/(s1)))
        }

        x1<-cbind(x1,x2)
      }

      if(n1>1){
        x<-t(rep(c(x1),times=((s1^(n1-e))*(s2^n2)*(s3^n3))))
      }else{ x<-t(rep(c(x1),times=1))}
      y11<-cbind(y11,t(x))
      j<-j+1
      e=e-1
    }
    z1<-y11
    b=1
    while (b<k) {

      y12<-matrix(,nrow= (s1^n1)*(s2^n2)*(s3^n3),ncol=n1)
      a=1

      while (a<=n1) {

        y12[,((a+(b-1))%%(n1))+1]=z1[,a]

        a=a+1
      }
      y11<-rbind(y11,y12)

      b=b+1
    }

  }else{
    y12<-matrix(,nrow = k*(s1^n1)*(s2^n2)*(s3^n3), ncol = 0)
  }




  if(s2>1){

    l2=as.integer(s2/2)
    v2<-seq(l2,-l2)
    if((s2)%%2==0){


      v2<-seq((s2-1),-(s2-1),by=-2)

    }
    v2<-t(v2)


    y21<-matrix(,nrow = (s1^n1)*(s2^n2)*(s3^n3), ncol = 0)
    j=1
    e=n2
    while(j<=n2 && e>0){
      x1<-matrix(,nrow = 1,ncol = 0)

      for (i in v2){

        if(n2>1){
          x2<-t(rep(c(i),times=s2^(n2-j)))
        }else{x2<-t(rep(c(i),times=N/s2))}

        x1<-cbind(x1,x2)

      }
      if(n2>1){
        x<-t(rep(c(x1),times=(s1^(n1))*(s2^(n2-e))*(s3^n3)))
      }else{x<-t(rep(c(x1),times=1))}

      y21<-cbind(y21,t(x))
      j<-j+1
      e=e-1
    }
    z2<-y21



    b=1
    while (b<k) {

      y22<-matrix(,nrow= (s1^n1)*(s2^n2)*(s3^n3),ncol=n2 )
      a=1

      while (a<=n2) {

        y22[,((a+(b-1))%%(n2))+1]=z2[,a]

        a=a+1
      }
      y21<-rbind(y21,y22)

      b=b+1
    }

  }else{
    y21<-matrix(,nrow = k*(s1^n1)*(s2^n2)*(s3^n3), ncol = 0)
  }



  if(s3>1){

    l3=as.integer(s3/2)
    v3<-seq(l3,-l3)
    if((s3)%%2==0){

      v3<-seq((s3-1),-(s3-1),by=-2)

    }
    v3<-t(v3)

    y31<-matrix(,nrow = (s1^n1)*(s2^n2)*(s3^n3), ncol = 0)
    j=1
    e=n3
    while(j<=n3 && e>0){
      x1<-matrix(,nrow = 1,ncol = 0)

      for (i in v3){
        if(n3>1){
          x2<-t(rep(c(i),times=s3^(n3-j)))
        }else{x2<-t(rep(c(i),times=N/s3))}

        x1<-cbind(x1,x2)

      }
      if(n3>1){
        x<-t(rep(c(x1),times=s1^n1*s2^n2*s3^(n3-e)))
      }else{x<-t(rep(c(x1),times=1))}

      y31<-cbind(y31,t(x))
      j<-j+1
      e=e-1
    }
    z3<-y31



    b=1
    while (b<k) {

      y32<-matrix(,nrow= (s1^n1)*(s2^n2)*(s3^n3),ncol=n3 )
      a=1

      while (a<=n3) {

        y32[,((a+(b-1))%%(n3))+1]=z3[,a]

        a=a+1
      }
      y31<-rbind(y31,y32)

      b=b+1
    }


  }else{
    y31<-matrix(,nrow = k*(s1^n1)*(s2^n2)*(s3^n3), ncol = 0)
  }
  X<-cbind(rep(c(1),times=k*(s1^n1)*(s2^n2)*(s3^n3)),y11,y21,y31)

  ##########################################

  if(s1>2 ){
    x1<-matrix(,nrow=k*(s1^n1)*(s2^n2)*(s3^n3),ncol=0)

    p=2
    while(p<s1){
      i=1+1
      while(i<=1+n1){
        x1<-cbind(x1,(X[,i])^p)
        i=i+1


      }
      p=p+1
    }
    X1<-x1
    X<-cbind(X,X1)
  }
  #####################

  if(s2>2 ){
    x1<-matrix(,nrow=k*(s1^n1)*(s2^n2)*(s3^n3),ncol=0)

    p=2
    while(p<s2){
      j=1+n1+1
      while(j<=1+n1+n2){
        x1<-cbind(x1,(X[,j])^p)
        j=j+1
      }
      p=p+1
    }
    X1<-x1
    X<-cbind(X,X1)
  }
  ####################

  if(s3>2 ){
    x1<-matrix(,nrow=k*(s1^n1)*(s2^n2)*(s3^n3),ncol=0)


    p=2
    while(p<s3){
      d=1+n1+n2+1
      while(d<=1+n1+n2+n3 ){
        x1<-cbind(x1,(X[,d])^p)
        d=d+1
      }
      p=p+1
    }
    X1<-x1
    X<-cbind(X,X1)
  }

  # X matrix

  X<-rbind(X[(k*(s1^n1)*(s2^n2)*(s3^n3)),],X,X[1,])
  print('X matrix')
  print(X)


  #####################################################
  # G matrix

  g=k*(s1^n1)*(s2^n2)*(s3^n3)

  mat<-matrix(0,nrow=k*(s1^n1)*(s2^n2)*(s3^n3),ncol=k*(s1^n1)*(s2^n2)*(s3^n3)+2)
  while(g>0){

    mat[c(g),c(g,g+1,g+2)]<-c(c,1,c)
    G<-mat

    g=g-1
  }

  #############


  #matrix multiplication Z=GX

  Z= G%*%X

  Z_prime_Z=t(Z)%*%Z

  #variance of estimated response
  k1=1

  var<-c()

  while(k1<=k*(s1^n1)*(s2^n2)*(s3^n3))
  {
    V=t(X[k1,])
    b<-t(V)
    v_y_hat<-V %*%solve(Z_prime_Z) %*% b

    var<-c(var,v_y_hat)


    k1<-k1+1
  }

  print("Z_prime_Z matrix")
  print(Z_prime_Z)
  print("inv(Z_prime_Z) matrix")
  print(solve(Z_prime_Z))
  variance_esitmated_response<-round(var,digits = 4 )

  vec1<-c(variance_esitmated_response)
  variance_of_estimated_response<-unique(vec1)

  vec2<-c(1:length(vec1))
  total_number_of_runs<-max(vec2)
  vec3<-c('total number of runs',max(vec2))
  vec4<-c('variance of estimated response',unique(vec1) )
  print(vec3)
  print(vec4)



}


















#'@name asym1
#'@aliases asym1
#'@title This generates a class of asymmetric rotatable response surface designs with neighbour effects under a second order model
#'@examples
#'library(rsdNE)
#' asym1(1,1,0.5)
#' ##X matrix
#'#      [,1] [,2] [,3] [,4]
#'#[1,]    1   -1   -1    1
#'#[2,]    1    1    1    1
#'#[3,]    1    1    0    0
#'#[4,]    1    1   -1    1
#'#[5,]    1   -1    1    1
#'#[6,]    1   -1    0    0
#'#[7,]    1   -1   -1    1
#'#[8,]    1    1    1    1
#'##Z prime Z matrix
#'#      [,1] [,2] [,3] [,4]
#'#[1,]   24    0    0   16
#'#[2,]    0   12    0    0
#'#[3,]    0    0    1    0
#'#[4,]   16    0    0   11
#'##Z prime Z imverse matrix
#'#      [,1]       [,2]   [,3]  [,4]
#'#[1,]  1.375 0.00000000    0   -2
#'#[2,]  0.000 0.08333333    0    0
#'#[3,]  0.000 0.00000000    1    0
#'#[4,] -2.000 0.00000000    0    3
#'#[1] "total number of runs" "6"
#'#[1] "variance of esitmated response" "1.4583"
#'@references Verma et al.2021, Communication in Statistics – Simulation and Computation
#' @description This function generates asymmetrical rotatable response surface designs in the presence of neighbour effects for 2n factors, n factors at 2 levels and another n factors at 3 levels.
#' @param n1 n1 factors having 2 levels, 1<=n1<=5
#' @param n2 n2 factors having 3 levels, 1<=n2<=5
#' @param c Value of alpha (Coefficient of neighbour effects), 0<=c<=1
#' @note Here 3 types of cases have been considered:
#' (2^n1*3^n2), where, n1=n2=n;
#' (2^n1*3), where, n1=n and n2=1;
#' (2*3^n2), where, n1=1 and n2=n.
#'@returns  This function generates rotatable designs as well as Z_prime_Z matrix,
#'inv(Z_primeZ) matrix and variance estimated response for the (2^n1 * 3^n2) factorial combination.
#'@author
#'Ashutosh Dalal, Division of Design of Experiments,ICAR-IASRI, New Delhi.
#'Seema Jaggi, Education Division, ICAR, Krishi Anusandhan Bhawan - II, Pusa, New Delhi.
#'Eldho Varghese,Fishery Resources Assessment Division,ICAR-CMFRI, Kochi.
#'Subhasish Sarkar, Division of Computer Application,ICAR-IASRI, New Delhi.
#'Arpan Bhowmik, Division of Design of Experiments,ICAR-IASRI, New Delhi.
#'Cini Varghese, Division of Design of Experiments,ICAR-IASRI, New Delhi.
#'Anindita Datta, Division of Design of Experiments,ICAR-IASRI, New Delhi.
#'Soumen Pal, Division of Computer Application,ICAR-IASRI, New Delhi.
#'@export
asym1<-function(n1,n2,c){


  # s1 no. of levels corresponding to the n1 no. of factors
  s1=2
  n1
  # s2 no. of levels corresponding to the n2 no. of factors
  s2=3
  n2

  N=(s1^n1)*(s2^n2)
  c #alpha value lies between 0 to 1
  k<-max(n1,n2)

  #generation of columns of the design matrix

  if(s1>1){
    l1=as.integer(s1/2)
    v1<-seq(l1,-l1)
    if((s1)%%2==0){


      v1<-seq((s1-1),-(s1-1),by=-2)

    }
    v1<-t(v1)


    y11<-matrix(,nrow = (s1^n1)*(s2^n2), ncol = 0)
    j=1
    e=n1
    while(j<=n1 && e>0){
      x1<-matrix(,nrow = 1,ncol = 0)
      for (i in v1){
        if(n1>1){
          x2<-t(rep(c(i),times=(s1^(n1-j))*(s2^n2)))
        }else{
          x2<-t(rep(c(i),times=N/(s1)))
        }

        x1<-cbind(x1,x2)

      }

      if(n1>1){
        x<-t(rep(c(x1),times=(N/((s1^e)*(s2^n2)))))
      }else{ x<-t(rep(c(x1),times=1))}

      y11<-cbind(y11,t(x))
      j<-j+1
      e=e-1
    }

    z1<-y11
    b=1
    while (b<k) {

      y12<-matrix(,nrow= (s1^n1)*(s2^n2),ncol=n1)
      a=1

      while (a<=n1) {

        y12[,((a+(b-1))%%(n1))+1]=z1[,a]

        a=a+1
      }
      y11<-rbind(y11,y12)

      b=b+1
    }

  }else{
    y12<-matrix(,nrow = k*(s1^n1)*(s2^n2), ncol = 0)
  }



  ######################################

  if(s2>1){

    l2=as.integer(s2/2)
    v2<-seq(l2,-l2)
    if((s2)%%2==0){


      v2<-seq((s2-1),-(s2-1),by=-2)

    }
    v2<-t(v2)


    y21<-matrix(,nrow = (s1^n1)*(s2^n2), ncol = 0)
    j=1
    e=n2
    while(j<=n2 && e>0){
      x1<-matrix(,nrow = 1,ncol = 0)

      for (i in v2){

        if(n2>1){
          x2<-t(rep(c(i),times=s2^(n2-j)))
        }else{x2<-t(rep(c(i),times=1))}

        x1<-cbind(x1,x2)

      }
      if(n2>1){
        x<-t(rep(c(x1),times=(s1^(n1))*(s2^(n2-e))))
      }else{x<-t(rep(c(x1),times=N/s2))}

      y21<-cbind(y21,t(x))
      j<-j+1
      e=e-1
    }

    z2<-y21



    b=1
    while (b<k) {

      y22<-matrix(,nrow= (s1^n1)*(s2^n2),ncol=n2 )
      a=1

      while (a<=n2) {

        y22[,((a+(b-1))%%(n2))+1]=z2[,a]

        a=a+1
      }
      y21<-rbind(y21,y22)

      b=b+1
    }


  }else{
    y21<-matrix(,nrow = k*(s1^n1)*(s2^n2), ncol = 0)
  }



  ############################

  X<-cbind(rep(c(1),times=k*(s1^n1)*(s2^n2)),y11,y21)

  ##########################################


  if(s1>2 ){
    x1<-matrix(,nrow=k*(s1^n1)*(s2^n2),ncol=0)

    p=2
    while(p<s1){
      i=1+1
      while(i<=1+n1){
        x1<-cbind(x1,(X[,i])^p)
        i=i+1


      }
      p=p+1
    }
    X1<-x1
    X<-cbind(X,X1)
  }
  #####################

  if(s2>2 ){
    x1<-matrix(,nrow=k*(s1^n1)*(s2^n2),ncol=0)

    p=2
    while(p<s2){
      j=1+n1+1
      while(j<=1+n1+n2){
        x1<-cbind(x1,(X[,j])^p)
        j=j+1
      }
      p=p+1
    }
    X1<-x1
    X<-cbind(X,X1)
  }
  ####################

  # generated X matrix

  X<-rbind(X[(k*(s1^n1)*(s2^n2)),],X,X[1,])
  print('X matrix')
  print(X)


  #####################################################
  # G matrix

  g=k*(s1^n1)*(s2^n2)

  mat<-matrix(0,nrow=k*(s1^n1)*(s2^n2),ncol=k*(s1^n1)*(s2^n2)+2)
  while(g>0){

    mat[c(g),c(g,g+1,g+2)]<-c(c,1,c)
    G<-mat

    g=g-1
  }
  #print(G)

  #matrix multiplication Z=GX

  Z= G%*%X

  Z_prime_Z=t(Z)%*%Z

  #variance of estimated response
  k1=1

  var<-c()

  while(k1<=k*(s1^n1)*(s2^n2))
  {
    V=t(X[k1,])
    b<-t(V)
    v_y_hat<-V %*%solve(Z_prime_Z) %*% b

    var<-c(var,v_y_hat)


    k1<-k1+1
  }
  print("Z_prime_Z matrix")
  print(Z_prime_Z)
  print("inv(Z_prime_Z) matrix")
  print(solve(Z_prime_Z))
  variance_esitmated_response<-round(var,digits = 4 )

  vec1<-c(variance_esitmated_response)
  variance_of_estimated_response<-unique(vec1)

  vec2<-c(1:length(vec1))
  total_number_of_runs<-max(vec2)
  vec3<-c('total number of runs',max(vec2))
  vec4<-c('variance of estimated response',unique(vec1) )
  print(vec3)
  print(vec4)


}


#' @param s1 Number of levels of n1 factors, 1<s1<=8
#' @param n1 Number of factors, 1<=n1<=4
#' @param s2 Number of levels of n2 factors, 1<s2<=8
#' @param n2 Number of factors, 1<=n2<=4
#' @param c Value of alpha (Coefficient of neighbour effects), 0<=c<=1
#'@name asym2
#'@aliases asym2
#'@title This generates a class of asymmetric rotatable response surface designs with neighbour effects under a polynomial model of order max(s1,s2)-1
#'@examples
#'library(rsdNE)
#'asym2(2,2,5,2,0.5)
#'@references Dalal, 2021, Unpublished M.Sc. Thesis, IARI, New Delhi
#' @description This function generates asymmetrical rotatable response surface designs in the presence of neighbour effects for (n1 + n2) factors, n1 factors at s1 levels and another n2 factors at s2 levels.
#'@note  Here s1 and s2 both not even at the same time and s1 not equal to s2.
#'@returns his function generates rotatable designs as well as Z_prime_Z matrix,
#'inv(Z_primeZ) matrix and variance estimated response for the (s1^n1 * s2^n2) factorial combination.
#'@author
#'Ashutosh Dalal, Division of Design of Experiments,ICAR-IASRI, New Delhi.
#'Seema Jaggi, Education Division, ICAR, Krishi Anusandhan Bhawan - II, Pusa, New Delhi.
#'Eldho Varghese,Fishery Resources Assessment Division,ICAR-CMFRI, Kochi.
#'Subhasish Sarkar, Division of Computer Application,ICAR-IASRI, New Delhi.
#'Arpan Bhowmik, Division of Design of Experiments,ICAR-IASRI, New Delhi.
#'Cini Varghese, Division of Design of Experiments,ICAR-IASRI, New Delhi.
#'Anindita Datta, Division of Design of Experiments,ICAR-IASRI, New Delhi.
#'Soumen Pal, Division of Computer Application,ICAR-IASRI, New Delhi.
#'@export

asym2<-function(s1,n1,s2,n2,c){

  # s1 no. of levels corresponding to the n1 no. of factors
  s1
  n1

  # s2 no. of levels corresponding to the n2 no. of factors
  s2
  n2
  # s3 no. of levels corresponding to the n3 no. of factors
  s3=1
  n3=0
  N=(s1^n1)*(s2^n2)*(s3^n3)

  #alpha value lies between 0 to 1
  c

  k<-max(n1,n2,n3)


  #generation of columns of the design matrix
  if(s1>1){
    l1=as.integer(s1/2)
    v1<-seq(l1,-l1)
    if((s1)%%2==0){

      v1<-seq((s1-1),-(s1-1),by=-2)

    }
    v1<-t(v1)




    y11<-matrix(,nrow = (s1^n1)*(s2^n2)*(s3^n3), ncol = 0)
    j=1
    e=n1
    while(j<=n1 && e>0){
      x1<-matrix(,nrow = 1,ncol = 0)
      for (i in v1){
        if(n1>1){
          x2<-t(rep(c(i),times=(s1^(n1-j))))
        }else{
          x2<-t(rep(c(i),times=N/(s1)))
        }

        x1<-cbind(x1,x2)
      }

      if(n1>1){
        x<-t(rep(c(x1),times=((s1^(n1-e))*(s2^n2)*(s3^n3))))
      }else{ x<-t(rep(c(x1),times=1))}
      y11<-cbind(y11,t(x))
      j<-j+1
      e=e-1
    }
    z1<-y11
    b=1
    while (b<k) {

      y12<-matrix(,nrow= (s1^n1)*(s2^n2)*(s3^n3),ncol=n1)
      a=1

      while (a<=n1) {

        y12[,((a+(b-1))%%(n1))+1]=z1[,a]

        a=a+1
      }
      y11<-rbind(y11,y12)

      b=b+1
    }

  }else{
    y12<-matrix(,nrow = k*(s1^n1)*(s2^n2)*(s3^n3), ncol = 0)
  }




  ######################################


  if(s2>1){

    l2=as.integer(s2/2)
    v2<-seq(l2,-l2)
    if((s2)%%2==0){


      v2<-seq((s2-1),-(s2-1),by=-2)

    }
    v2<-t(v2)


    y21<-matrix(,nrow = (s1^n1)*(s2^n2)*(s3^n3), ncol = 0)
    j=1
    e=n2
    while(j<=n2 && e>0){
      x1<-matrix(,nrow = 1,ncol = 0)

      for (i in v2){

        if(n2>1){
          x2<-t(rep(c(i),times=s2^(n2-j)))
        }else{x2<-t(rep(c(i),times=N/s2))}

        x1<-cbind(x1,x2)

      }
      if(n2>1){
        x<-t(rep(c(x1),times=(s1^(n1))*(s2^(n2-e))*(s3^n3)))
      }else{x<-t(rep(c(x1),times=1))}

      y21<-cbind(y21,t(x))
      j<-j+1
      e=e-1
    }
    z2<-y21



    b=1
    while (b<k) {

      y22<-matrix(,nrow= (s1^n1)*(s2^n2)*(s3^n3),ncol=n2 )
      a=1

      while (a<=n2) {

        y22[,((a+(b-1))%%(n2))+1]=z2[,a]

        a=a+1
      }
      y21<-rbind(y21,y22)

      b=b+1
    }

  }else{
    y21<-matrix(,nrow = k*(s1^n1)*(s2^n2)*(s3^n3), ncol = 0)
  }



  ############################


  if(s3>1){

    l3=as.integer(s3/2)
    v3<-seq(l3,-l3)
    if((s3)%%2==0){

      v3<-seq((s3-1),-(s3-1),by=-2)

    }
    v3<-t(v3)

    y31<-matrix(,nrow = (s1^n1)*(s2^n2)*(s3^n3), ncol = 0)
    j=1
    e=n3
    while(j<=n3 && e>0){
      x1<-matrix(,nrow = 1,ncol = 0)

      for (i in v3){
        if(n3>1){
          x2<-t(rep(c(i),times=s3^(n3-j)))
        }else{x2<-t(rep(c(i),times=N/s3))}

        x1<-cbind(x1,x2)

      }
      if(n3>1){
        x<-t(rep(c(x1),times=s1^n1*s2^n2*s3^(n3-e)))
      }else{x<-t(rep(c(x1),times=1))}

      y31<-cbind(y31,t(x))
      j<-j+1
      e=e-1
    }
    z3<-y31



    b=1
    while (b<k) {

      y32<-matrix(,nrow= (s1^n1)*(s2^n2)*(s3^n3),ncol=n3 )
      a=1

      while (a<=n3) {

        y32[,((a+(b-1))%%(n3))+1]=z3[,a]

        a=a+1
      }
      y31<-rbind(y31,y32)

      b=b+1
    }


  }else{
    y31<-matrix(,nrow = k*(s1^n1)*(s2^n2)*(s3^n3), ncol = 0)
  }
  X<-cbind(rep(c(1),times=k*(s1^n1)*(s2^n2)*(s3^n3)),y11,y21,y31)

  ##########################################

  if(s1>2 ){
    x1<-matrix(,nrow=k*(s1^n1)*(s2^n2)*(s3^n3),ncol=0)

    p=2
    while(p<s1){
      i=1+1
      while(i<=1+n1){
        x1<-cbind(x1,(X[,i])^p)
        i=i+1


      }
      p=p+1
    }
    X1<-x1
    X<-cbind(X,X1)
  }
  #####################

  if(s2>2 ){
    x1<-matrix(,nrow=k*(s1^n1)*(s2^n2)*(s3^n3),ncol=0)

    p=2
    while(p<s2){
      j=1+n1+1
      while(j<=1+n1+n2){
        x1<-cbind(x1,(X[,j])^p)
        j=j+1
      }
      p=p+1
    }
    X1<-x1
    X<-cbind(X,X1)
  }
  ####################

  if(s3>2 ){
    x1<-matrix(,nrow=k*(s1^n1)*(s2^n2)*(s3^n3),ncol=0)


    p=2
    while(p<s3){
      d=1+n1+n2+1
      while(d<=1+n1+n2+n3 ){
        x1<-cbind(x1,(X[,d])^p)
        d=d+1
      }
      p=p+1
    }
    X1<-x1
    X<-cbind(X,X1)
  }

  # X matrix

  X<-rbind(X[(k*(s1^n1)*(s2^n2)*(s3^n3)),],X,X[1,])
  print(X)


  #####################################################
  # G matrix

  g=k*(s1^n1)*(s2^n2)*(s3^n3)

  mat<-matrix(0,nrow=k*(s1^n1)*(s2^n2)*(s3^n3),ncol=k*(s1^n1)*(s2^n2)*(s3^n3)+2)
  while(g>0){

    mat[c(g),c(g,g+1,g+2)]<-c(c,1,c)
    G<-mat

    g=g-1
  }

  #############


  #matrix multiplication Z=GX

  Z= G%*%X

  Z_prime_Z=t(Z)%*%Z

  #variance of estimated response
  k1=1

  var<-c()

  while(k1<=k*(s1^n1)*(s2^n2)*(s3^n3))
  {
    V=t(X[k1,])
    b<-t(V)
    v_y_hat<-V %*%solve(Z_prime_Z) %*% b

    var<-c(var,v_y_hat)


    k1<-k1+1
  }

print("Z_prime_matrix")
  print(Z_prime_Z)
  print("inv(Z_prime_Z) matrix")
  print(solve(Z_prime_Z))
  variance_esitmated_response<-round(var,digits = 4 )

  vec1<-c(variance_esitmated_response)
  variance_of_estimated_response<-unique(vec1)

  vec2<-c(1:length(vec1))
  total_number_of_runs<-max(vec2)
  vec3<-c('total number of runs',max(vec2))
  vec4<-c('variance of estimated response',unique(vec1))
  print(vec3)
  print(vec4)
}

