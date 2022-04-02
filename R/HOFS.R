HOFS <- function( data, label, C = 0.2, Nfeatures = 10, samplePct = 0.1, verbose = T, seed = 666 ){

  set.seed( seed )

  # Assing global variables
  sapply( c( "MIL", "MIH", "MIHXY", "COV", "apI" ), function( x ){ assign( x, NULL, envir = .GlobalEnv ) } )
  sapply( c( "Zica" ), function( x ){ assign( x, list(), envir = .GlobalEnv ) } )

  # Attribute preparation, i.e. normalziation
  temp <- DataPreprocess( data )
  data <- temp$data
  idx <- temp$idx

  temp <- DataPreprocess( label )
  label <- temp$data

  Nfeature <- ncol( data )
  Nsamples <- nrow( data )
  N <- min( 50, length( idx ) )

  Stmp <- Atmp <- rank <- NULL; tmp <- 1; s <- 0
  selected <- list()

  if( verbose ) { print( "Preparing single feature entropy" ) }

  for( i in 1:Nfeature ){
    MIL[ i ] <<- get_I( data[, i, drop = F ], label )
    MIH[ i ] <<- get_H_mix( data[, i, drop = F ], 1, nrow( data ) )
    MIHXY[ i ] <<- get_H_mix( cbind( data[, i, drop = F ], label ), 2, nrow( label ) )
    MIY <- get_H_mix( label, 1, nrow( data ) )
  }

  combo <- 1
  W <- list()

  temp <- get_I_selected_ICA2( data, label, selected, samplePct )
  index <- temp$idx

  Wout <- temp$Wout

  if( verbose ) { print( paste( 1,' feature selected' ) ) }
  s <- s + 1
  selected[[ combo ]] <- index
  rank <- c( rank, index )

  while( s < N & s < Nfeatures ){

    if( verbose ) { print( paste( s + 1,' feature selected' ) ) }

    temp <- get_I_selected_ICA2( data, label, selected, samplePct )

    index <- temp$idx
    W[[ combo ]] <- temp$Wout

    rank <- c( rank, index )
    argcov <- NULL

    for( i in 1:combo ){

      argcov[ i ] <- 0

      Nfea <- length( selected[[ i ]] )

      for( j in 1:Nfea ){

        C_ <- cor( cbind( data[ , index, drop = F ], data[ , selected[[ i ]][ j, drop = F ] ] ) )
        argcov[ i ] <- argcov[ i ] + C_[ 1, 2 ]# abs()

      }

      argcov[ i ] <- argcov[ i ] / Nfea

    }

    MaxC <- max( argcov )

    if( MaxC <= C ){

      Atmp <- NULL
      n <- length( Zica )

      temp <- Zica[[ n ]]
      Zica[[ n + 1 ]] <<- temp[ nrow( temp ), , drop = F ]

      temp <- Zica[[ n ]]
      Zica[[ n ]] <<- temp[ -nrow( temp ), , drop = F ]

      temp <- W[[ combo ]]
      W[[ combo + 1 ]] <- temp[ nrow( temp ), ncol( temp ), drop = F ]

      temp <- W[[ combo ]]
      W[[ combo ]] <- temp[ 1:( nrow( temp ) - 1 ), 1:( ncol( temp ) - 1 ), drop = F ]

      combo <- combo + 1
      selected[[ combo ]] <- index

      s <- s + 1

    }else{

      subidx <- which( argcov == MaxC )
      selected[[ subidx ]] <- c( selected[[ subidx ]], index )

      s <- s + 1

    }

  }

  MIL <- MIL; MIH <- MIH; MIHXY <- MIHXY; MIY <- MIY
  if( !is.null( colnames( data ) ) ){

    names(MIL) <- names(MIH) <- names(MIHXY) <- colnames( data )
    names(MIY) <- colnames( label )

    for( i in 1:length(selected) ){
      temp <- selected[[i]]
      names(temp) <- colnames(data)[temp]
      selected[[i]] <- temp
    }

  }

  # Remove global variables
  rm( list = c( "MIL", "MIH", "MIHXY", "COV", "Zica", "apI" ), envir = .GlobalEnv )

  names( selected ) <- paste0( "Subset_", 1:length(selected) )
  return( list( Selected = selected, Rank = rank, EntropyY = MIY, EntropyX = MIH[rank],
                Joint_EntropyXY = MIHXY[rank], Mutual_InformationXY = MIL[rank] ) )
}

