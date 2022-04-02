
# Data normalization into the range [0-1]
NormData <- function( data ){

  # Adjustment for dividing by 0
  epi <- 0.000000001

  # Normalization
  out <- apply( data, 2, function( x, e ){ ( x - min( x ) ) / ( max( x ) + e - min( x ) )}, e = epi )

  return( out )

}

DataPreprocess <- function( data ){

  factorTOnumeric <- function( factor, numbers ) {
    num <- NA
    for ( i in 1:length( levels(factor) ) ) {
      num <- ifelse( factor == levels(factor)[i], numbers[i], num )
    }
    return( num )
  }
  characterTOnumeric <- function( char, numbers ) {
    num <- NA
    val <- unique( char )
    for(i in 1:length( val ) ){
      num <- ifelse( char == val[i], numbers[i], num )
    }
    return( num )
  }

  idx <- c()

  for( i in 1:ncol( data ) ){

    if( is.factor( data[,i] ) ){

      data[, i ] <- factorTOnumeric( data[, i], 1:length( levels( data[, i] ) ) )
      idx <- c( idx, i )

    }else if( is.character( data[,i] ) ){

      data[, i ] <- characterTOnumeric( data[, i], 1:length( unique( data[, i] ) ) )
      idx <- c( idx, i )

    }else{

      max_ <- max( data[,i] )
      min_ <- min( data[,i] )
      if( ( max_- min_ ) >= 0.1 ){

        idx <- c( idx, i )

      }

    }

  }

  # Take only non constant attributes
  data <- data[ , idx, drop = F ]

  # Data normalization into the range [0-1]
  out <- NormData( data )

  return( list( data = out, idx = idx ) )

}

# formula 4: Mutual Information (MI) I(X:Y)
get_I <- function( x, y, cons = NULL ){

  if( is.null( cons ) ){

    data1 <- x
    data2 <- y
    I <- get_H_mix( x, ncol( x ), nrow( x ) ) - get_H_cond( data1, data2 )

  }else{

    I <- get_H_cond( x, cons ) - get_H_cond( x, cbind( y, cons ) )

  }

  return( I )

}

# formula 3: Conditional Entropy H(X|Y)
get_H_cond <- function( x, y ){

  data1 <- cbind( x, y )
  data2 <- y

  H <- get_H_mix( data1, ncol( data1 ), nrow( data1 ) ) -
    get_H_mix( data2, ncol( data2 ), nrow( data2 ) )

  return( H )

}

normaliseArray <- function( inp ){

  minval <- min( inp )
  maxval <- max( inp )
  normvec <- floor( ( inp - minval ) / 0.2 )
  # normvec <- floor( ( inp - minval ) / 0.01 )

  num_state <- ( maxval - minval ) / 0.2 + 1
  # num_state <- ( maxval - minval ) / 0.01

  num_state <- floor( num_state )

  return( list( num_state = num_state, normvec = normvec ) )

}

# formula 1: Entropy H(X)
get_H_mix <- function( x, n_features, n_sample ){

  joint_states <- 1
  num_states <- integer()
  norm_vec <- list()

  for( i in 1:n_features ){

    temp <- normaliseArray( x[, i, drop = F ] )

    num_states[ i ] <- temp$num_state
    norm_vec[[ i ]] <- temp$normvec

    joint_states <- joint_states * num_states[ i ]

  }

  bias <- 1
  if( n_features > 1 ){

    for( i in 2:n_features ){

      bias[ i ] <- bias[ i-1 ] * num_states[ i-1 ]

    }

  }

  jointStateCounts <- integer( joint_states )
  for( i in 1:n_sample ){

    idx <- 1

    for( j in 1:n_features ){

      idx <- idx + norm_vec[[ j ]][ i ] * bias[ j ]

    }

    jointStateCounts[ idx ] <- jointStateCounts[ idx ] + 1

  }

  jointStateProbs <- jointStateCounts / n_sample
  epi <- 0.0000001
  val <- 0
  H <- 0

  for( i in 1:joint_states ){

    val <- jointStateProbs[ i ]
    if( val > 0 ){

      H <- H - val * log( val )

    }

  }

  return( H )

}

get_I_selected_ICA2 <- function( data, label, selected = NULL, samplePct = 0.1 ){

  Nstructure <- length( selected )
  I <- matrix( -Inf, ncol( data ), 1 )

  h <- 0
  W <- S1 <- W1 <- list()

  for( index in 1:ncol( data ) ){

    h[ index ] <- 0; tmp <- 0; flag <- 0

    # Filter assesing whether process current variable or skip
    if( Nstructure > 0 ){

      for( i in 1:Nstructure ){

        if( any( selected[[ i ]] == index ) ){

          flag <- 1
          break

        }

      }

    }

    if( flag == 0 ){

      h1 <- MIL[ index ]

      if( Nstructure > 0 ){

        for( i in 1:Nstructure ) { #formula 20: iteration over K independent subsets

          tmp <- tmp + length( selected[[ i ]] )

          curdata1 <- cbind( data[, selected[[ i ]], drop = F ], data[, index, drop = F ] )

          temp <- myICA_lila( t( curdata1 ), ncol( curdata1 ), 1, index, i, samplePct )
          S1[[ index ]] <- temp$Zica
          A <- temp$A
          TT <- temp$TT
          mu <- temp$mu

          curdata2 <- cbind( curdata1, label )
          W1[[ index ]] <- A %*% TT

          temp <- myICA_lila( t( curdata2 ), ncol( curdata2 ), 2, index, i, samplePct )
          S2 <- temp$Zica
          A <- temp$A
          TT <- temp$TT
          mu <- temp$mu

          W[[ index ]] <- A %*% TT
          S2 <- t( S2 )

          h[ index ] <- h[ index ] -
            MIH[ index ] +
            MIHXY[ index ] +
            get_H_mix( S2[, ncol( S2 ), drop = F ], 1, nrow( data ) ) +
            log( abs( W[[ index ]][ nrow( W[[ index ]] ), ncol( W[[ index ]] ) ] ) )

        }

      }

      if( h[ index ] == 0 ){

        I[ index ] <- h1

      }else{

        # formula 20: h1 = I(x_i :y) = MIL
        I[ index ] <- h1 + h[ index ]

      }

    }

  }

  idxtmp <- which( I == max( I ) )

  idx <- idxtmp[ 1 ]
  assign( "apI", c( apI, max( I ) ), envir = .GlobalEnv )

  if( length( W1 ) == 0 ){

    temp <- incre_ICA( t( data[ , idx, drop = F ] ), W = NULL, samplePct = samplePct )
    Zica[[ 1 ]] <<- temp$Zica
    Wout <- temp$A

  }else{

    Wout <- W1[[ idx ]]
    Zica[[ length( Zica ) ]] <<- S1[[ idx ]]

  }

  return( list( idx = idx, Wout = Wout ) )

}

incremental_ICA <- function( W_ , X ){

  n <- nrow( W_ )
  W <- matrix( runif( n+1 ), n+1, n+1 )

  W[ 1:n, 1:n ] <- W_
  wtmp <- runif( n+1 )

  detW <- det( W_ )
  lr <- 0.0001
  m <- ncol( X )
  maxIter <- 100
  batchsize <- max( floor( m / 10 ) )

  for( i in 1:maxIter ){

    dia <- apply( W_ * diag( n ), 2, sum )
    dia <- c( dia, 1 )
    dia[ seq(2, length( dia ), by = 2 ) ] <- -dia[ seq(2, length( dia ), by = 2 ) ]

    rd <- sample( c( 1, m ), batchsize, replace = T )

    x <- X[, rd, drop = F ]
    tmp <- ( 1 - 2 * g( ( wtmp * x ) ) )

    step <- tmp %*% t( x ) + repmat( detW %*% W[ n+1, n+1 ], 1, n+1 ) / dia
    wtmp <- wtmp + lr * step
  }

  W[ nrow( W ) ,] <- wtmp
  S <- W %*% X

  return( list( S = S, W = W ) )

}

incre_ICA <- function( Z, W = NULL, samplePct = 0.1 ){

  maxSamples <- 1000

  temp <- myCenter( Z )
  Zc <- temp$Zc
  mu <- temp$mu

  d <- dim( Zc )[ 1 ]
  n <- dim( Zc )[ 2 ]

  if ( n > maxSamples ){

    Zct <- Zc[, sample( 1:n, n * samplePct ) ]#

  }else{

    Zct <- Zc

  }

  if( is.null( W ) ){

    temp <- compute_W( Zc, samplePct )
    W <- temp$W

  }else{

    W <- mod_W( W, Zc, samplePct )

  }

  A <- W
  Zica <- A %*% Zc

  return( list( Zica = Zica, A = A, mu = mu ) )

}

compute_W <- function( X, samplePct ){

  W <- matrix( runif( nrow( X ) * nrow( X ) ), nrow( X ), nrow( X ) )

  lr <- 0.003
  m <- ncol( X )
  maxIter <- 30
  batchsize <- m * samplePct

  for( i in 1:maxIter ){

    rd <- matrix( sample( 1:m, batchsize ), ncol = 1 )
    x <- X[, rd, drop = F ]

    tmp <- ( 1 - 2 * g( ( W %*% x ) ) )

    step <- tmp %*% t( x ) + solve( t( W ) )

    W <- W + lr * step

  }

  S <- W %*% X

  return( list( S = S, W = W ) )

}

mod_W <- function( w, X, samplePct ){

  W <- matrix( 0, nrow( X ), nrow( X ) )
  W[1:nrow( w ), 1:ncol( w ) ] <- w

  W[ nrow( W ), ] <- runif( nrow( X ) )
  W[ nrow( W ), ncol( W ) ] <- 1

  lr <- 0.001
  m <- ncol( X )
  maxIter <- 100
  batchsize <- m * samplePct

  for( i in 1:maxIter ){

    rd <- matrix( sample( 1:m, batchsize ), ncol = 1 )
    x <- X[, rd, drop = F ]

    tmp <- ( 1 - 2 * g( ( W %*% x ) ) )

    step <- tmp * t( x ) + inv( t( W ) )

    W <- W + lr * step

    W[ nrow( W ) - 1, ncol( W ) ] <- 0
    W[ 1:nrow( w ), 1:ncol( w ) ] <- w

  }

  S <- W %*% X

  return( list( S = S, W = W ) )

}

g <- function( x ){
  1 / ( 1 + exp( -x ) )
}

myICA_lila <- function( Z, r, o, m, l, samplePct ){

  assign( "coun" , 1, envir = .GlobalEnv )

  eps <- 1e-6
  maxSamples <- 1000
  maxIters <- 20

  temp <- myCenterAndWhiten( Z )
  Zcw <- temp$Zcw
  TT <- temp$TT
  mu <- temp$mu

  d <- dim( Zcw )[ 1 ]
  n <- dim( Zcw )[ 2 ]

  if( n > maxSamples ){

    Zcwt_indx <- sample( 1:n, n * samplePct )
    Zcwt <- Zcw[, Zcwt_indx ]

  }else{

    Zcwt <- Zcw

  }

  normRows <- function( X ){
    out <- X * 1 / sqrt( apply( X^2, 1, sum ) )
    out[is.na(out)] <- 0
    out
  }

  W <- normRows( matrix( runif( r * d ), r, d ) )

  k <- 0
  err <- Inf

  if( coun == 54 ){

    f <- 0

  }

  while( err > eps && k < maxIters ){

    k <- k + 1

    Wlast <- W

    Sk <- Wlast %*% Zcwt
    Sk <- array( Sk, c( nrow( Sk ), 1, ncol( Sk ) ) )

    G <- Sk * exp( -0.5 * Sk^2 )
    Gp <- Sk * G

    G_perm <- G[,,]
    Zcwt_perm <- t( Zcwt )
    G_Zcwt_perm <- array_mult( G_perm, Zcwt_perm )
    # W <- Reduce( "+", G_Zcwt_perm ) / length( G_Zcwt_perm ) + apply( Gp, 1 , mean ) * Wlast
    W <- Reduce_cpp( G_Zcwt_perm ) / length( G_Zcwt_perm ) + apply( Gp, 1 , mean ) * Wlast

    W <- normRows( W )


    b1 <- which( is.nan( W ) )
    b2 <- which( W == Inf )
    b <- length( b1 ) + length( b2 )

    if( b > 0 ){

      h <- 0

    }

    temp <- try( svd( W ), T )

    U <- temp$u
    S <- diag( temp$d )

    W <- U %*% diag( 1 / diag( S ) ) %*% t( U ) %*% W
    W[is.na(W)] <- 0

    err <- max( 1 - sapply( 1:nrow(W), function( i, x, y ){ sum( x[i,] * y[i,] ) }, x = W, y = Wlast ) )

  }

  A <- W
  Zica <- A %*% Zcw

  return( list( Zica = Zica, A = A, TT = TT, mu = mu ) )

}

myCenterAndWhiten <- function( Z ){

  temp <- myCenter( Z )
  Zc <- temp$Zc
  mu <- temp$mu

  temp <- myWhiten( Zc )
  Zcw <- temp$Zw
  TT <- temp$TT

  return( list( Zcw = Zcw, TT = TT, mu = mu ) )

}

myCenter <- function( Z ){

  mu <- apply( Z, 1, mean  )

  Zc <- Z - mu

  return( list( Zc = Zc, mu = mu ) )

}

myWhiten <- function( Z ){

  R <- cov( t( Z ) )

  temp <- svd( R )
  U <- temp$u
  S <- diag( temp$d )

  temp <- diag( 1 / sqrt( diag( S ) ) )
  if( any( !is.finite( temp ) ) ){

    temp[ !is.finite( temp ) ] <- 0

  }
  TT <- U %*% temp %*% t( U )
  Zw <- TT %*% Z

  return( list( Zw = Zw, TT = TT ) )

}
