/***********************************************************************************
 * quad.cpp
 * 
 * Code to generate Gaussian quadrature nodes and weights
 * 
 * 30jan2016
 * Philip Barrett << Chicago
 * 
 ***********************************************************************************/

#include "quad.hpp"

// [[Rcpp::export]]
arma::vec quad_nodes_1d( int n_nodes, double sig=1, double mu=0 ){
// Computes the nodes for one-dimensional Guassian quadrature
  vec nodes = zeros( n_nodes );
  switch ( n_nodes ) 
  {
  case 1:
    nodes << 0;
    break;
  case 2:
    nodes << 0.7071067811865475 << -0.7071067811865475;
    break;     
  case 3:
    nodes << 1.224744871391589 << 0 << -1.224744871391589;
    break;
  case 4:
    nodes << 1.650680123885785 << 0.5246476232752903 << -0.5246476232752903 << -1.650680123885785;
    break;
  case 5:
    nodes << 2.020182870456086 << 0.9585724646138185 << 0 << -0.9585724646138185 << -2.020182870456086;
    break;
  case 6:
    nodes << 2.350604973674492 << 1.335849074013697 << 0.4360774119276165 << -0.4360774119276165 << -1.335849074013697 << -2.350604973674492;
    break;
  case 7:
    nodes << 2.651961356835233 << 1.673551628767471 << 0.8162878828589647 << 0 << -0.8162878828589647 << -1.673551628767471 << -2.651961356835233;
    break;
  case 8:
    nodes << 2.930637420257244 << 1.981656756695843 << 1.157193712446780 << 0.3811869902073221 << -0.3811869902073221 << -1.157193712446780 << -1.981656756695843 << -2.930637420257244;
    break;
  case 9:
    nodes << 3.190993201781528 << 2.266580584531843 << 1.468553289216668 << 0.7235510187528376 << 0 << -0.7235510187528376 << -1.468553289216668 << -2.266580584531843 << -3.190993201781528;
    break;
  case 10:
    nodes << 3.436159118837738 << 2.532731674232790 << 1.756683649299882 << 1.036610829789514 << 0.3429013272237046 << -0.3429013272237046 << -1.036610829789514 << -1.756683649299882 << -2.532731674232790 << -3.436159118837738;
    break;
  }
  return mu + sqrt(2.0) * sig * nodes ;
      // Scale nodes so that they can be used directly in integration
}


// [[Rcpp::export]]
arma::vec quad_weights_1d( int n_nodes ){
// Computes the weights for one-dimensional Guassian quadrature
  vec weights = zeros( n_nodes );
  switch ( n_nodes ) 
  {
  case 1:
    weights << 1.77245385;   /* Sqrt(pi) */
    break;
  case 2:
    weights << 0.8862269254527580 <<  0.8862269254527580;
    break;     
  case 3:
    weights << 0.2954089751509193 << 1.181635900603677 << 0.2954089751509193;
    break;
  case 4:
    weights << 0.08131283544724518 << 0.8049140900055128 << 0.8049140900055128 << 0.08131283544724518;
    break;
  case 5:
    weights << 0.01995324205904591 << 0.3936193231522412 << 0.9453087204829419 << 0.3936193231522412 << 0.01995324205904591;
    break;
  case 6:
    weights << 0.004530009905508846 << 0.1570673203228566 << 0.7246295952243925 << 0.7246295952243925 << 0.1570673203228566 << 0.004530009905508846;
    break;
  case 7:
    weights << 0.0009717812450995192 << 0.05451558281912703 << 0.4256072526101278 << 0.8102646175568073 << 0.4256072526101278 << 0.05451558281912703 << 0.0009717812450995192;
    break;
  case 8:
    weights << 0.0001996040722113676 << 0.01707798300741348 << 0.2078023258148919 << 0.6611470125582413 << 0.6611470125582413 << 0.2078023258148919 << 0.01707798300741348 << 0.0001996040722113676;
    break;
  case 9:
    weights << 0.00003960697726326438 << 0.004943624275536947 << 0.08847452739437657 << 0.4326515590025558 << 0.7202352156060510 << 0.4326515590025558 << 0.08847452739437657 << 0.004943624275536947 << 0.00003960697726326438;    
    break;
  case 10:
    weights << 7.640432855232621e-06 << 0.001343645746781233 << 0.03387439445548106 << 0.2401386110823147 << 0.6108626337353258 << 0.6108626337353258 << 0.2401386110823147 << 0.03387439445548106 << 0.001343645746781233 << 7.640432855232621e-06;
    break;
  }
  weights = weights / 1.77245385;
      // Divide weights by sqrt(pi) (also rescales them to sum to 1)
  return weights;
}

// [[Rcpp::export]]
List quad_nodes_weights( int n_nodes, int n_dim, 
                          arma::vec sig, arma::vec mu ){
// Computes the integration nodes in n_dim dimensions
  mat nodes = zeros( pow( n_nodes, n_dim ), n_dim ) ;
  vec weights = ones( pow( n_nodes, n_dim ) ) ;
      // Initialize the nodes and weights
  vec nodes_1d = quad_nodes_1d( n_nodes, 1, 0 ) ;
  vec weights_1d = quad_weights_1d( n_nodes ) ;
      // The one-dimensional nodes and weights.  Do not scale the nodes yet -
      // needs doing on a dimension-by-dimension basis.
  umat idx = zeros<umat>( pow( n_nodes, n_dim ), n_dim ) ;
      // The matrix of indices
  for( int i = 1 ; i < pow( n_nodes, n_dim ) ; i++ ){
  // Loop to compute the indices of the quadrature nodes
    idx.row(i) = idx.row(i-1) ;
        // Start out by copying the preceiding row down
    if( idx( i - 1, n_dim - 1 ) < n_nodes - 1 ){
      idx( i, n_dim - 1 ) = idx( i - 1, n_dim - 1 ) + 1 ;
    }   // Increment the last index where possible
    else
    {
      int j = n_dim - 2 ;
      bool complete = false ;
          // Initiate iteration over row elements
      while( j >=0 & !complete ){
        if( idx( i, j ) < n_nodes - 1 ){
        // If the current trial case is something that can be incremented
          idx( i, j ) ++ ;
          idx.row(i).tail( n_dim - 1 - j ) = zeros<urowvec>(n_dim - 1 - j) ;
              // Increment and set subsequent indices to zero
          complete = true ;
              // Break the loop
        }
        else
        {
          j-- ;
              // Keep looking
        }
      }
    }
  }
  
  for( int i = 0 ; i < idx.n_rows ; i++ ){
    for( int j = 0 ; j < idx.n_cols ; j++ ){
      nodes(i,j) = mu(j) + sig(j) * nodes_1d(idx(i,j)) ;
          // Fill in the nodes matrix
      weights(i) = weights(i) * weights_1d(idx(i,j)) ;
    }
  }
  
  List out ;
  out["nodes"] = nodes ;
  out["weights"] = weights ;
  out["idx"] = idx ;
      // Create output
  return out ;
}

// [[Rcpp::export]]
arma::mat quad_nodes_weights_mat( int n_nodes, int n_dim, 
                          arma::vec sig, arma::vec mu ){
// Matrix form of the weights and nodes
  List l_out = quad_nodes_weights( n_nodes, n_dim, sig, mu ) ;
  vec weights = l_out["weights"] ;
  mat nodes = l_out["nodes"] ;
      // Extract the list elements
  mat out( nodes.n_rows, n_dim + 1 ) ;
      // Initialize the output matrix
  out.head_cols( n_dim ) = nodes ;
  out.col( n_dim ) = weights ;
      // Assign the weights and nodes to the different columns
  return out ;
}