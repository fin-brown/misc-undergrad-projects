
polyopt <- function(a, x0, tol=1e-6, MaxIter=100) {
  ####################################################################################
  ##
  ##    An implementation of the Newton-Raphson minimisation method for use on 
  ##    polynomial functions. 
  ##    
  ##    The arguments are:
  ##
  ##    a :        A vector containing the coefficients of the powers of x contained 
  ##               in the polynomial, where a[1] is the coefficient of x^0 (a constant)
  ##               and a[length(a)] is the coefficient of the highest power of x.
  ##               The terms of a must be non-zero and a must be at least 3 terms long.
  ##    x0 :       The initial estimate of the optimising x-value.
  ##    tol :      The tolerance for the stopping conditions of the iteration; when 
  ##               this exceeds either the gradient or the relative change between two
  ##               consecutive estimates the iteration halts. Defaults to 1e-6.
  ##    MaxIter :  The maximum number of iterations allowed before the process halts
  ##
  ##    The function returns a named list containing:
  ##    
  ##    x :        The final estimate of the optimising x-value
  ##    hx :       The value of the function at this optimising x-value. 
  ##    gradient : The gradient of the function at this optimising x-value.
  ##    hessian :  The value of the second derivative at this optimising x-value.
  ##    N.iter :   The number of iterations of the process before termination. 
  ##
  ####################################################################################
  #
  #  Before performing any computation check the arguments are valid. If we are given 
  #  fewer than 3 coefficients throw an error with an appropriate message, if any of 
  #  the coefficients are equal to 0 give a warning but continue the computation.
  #
  if(length(a)<3) stop('The polynomial to be minimised must be at least degree 2.')
  if(any(a==0)) warning('The input polynomial coefficients cannot be equal to 0.')
  #
  #  Start by setting x as the initial estimate x0 given in the arguments and defining 
  #  the variable p as the degree of the polynomial to estimate.
  #
  x <- x0; 
  p <- length(a)-1
  # 
  #  Define h as a function so we can evaluate it with different x-values in the same
  #  scope later on. We use vector multiplication to give the sum of the powers of x
  #  multiplied by their respective coefficients stored in the initial argument 
  #  vector a. 
  #
  h <- function(x) {x^(0:p) %*% a}
  #
  #  Define the vector a1 as the coefficients of h', the coefficient for the highest 
  #  power of x needs to be set to 0 as this power is longer present in h'. Again we 
  #  we use an expression when defining dh.dx, although here it is also allows more  
  #  succinct code further on.
  #
  a1    <- c((a*0:p)[-1], 0) 
  dh.dx <- expression(x^(0:p) %*% a1)
  #
  #  As above define the vector a2 as the coefficients of h" use this to construct
  #  the expression for d2h.dx2.
  #
  a2      <- c((a1 * c(0:(p-1), 0))[-1], 0)
  d2h.dx2 <- expression(x^(0:p) %*% a2) 
  #
  #  Now we need to initialise Iter, the iteration number, to 0 as we have not started
  #  iterating and RelChange to Inf (a place-holder for RelChange which is meaningless
  #  until we have our first estimate). RelChange is the difference between two 
  #  successive estimates and setting it to a large number means it will exceed the
  #  tolerance.
  #
  Iter      <- 0
  RelChange <- Inf
  #
  #  We use a while loop to perform the iteration with conditions that mean it will 
  #  continue to iterate until either dh.dx or RelChange become less than tol, the 
  #  tolerance supplied in the initial arguments (defaulting to 1e-6), or Iter 
  #  reaches MaxIter, the total number of allowed iterations supplied in the initial
  #  arguments (defaulting to 100). 
  #
  while (abs(eval(dh.dx)) > tol & abs(RelChange) > tol & Iter < MaxIter) {
    x.previous <- x                                     # Store the previous estimate. 
    x <- x - (eval(dh.dx) / eval(d2h.dx2))              # Calculate the new estimate.      
    RelChange <- (h(x) - h(x.previous)) / h(x.previous) # Calculate RelChange.
    Iter <- Iter + 1                                    # Move onto next iteration.
  }
  #
  #  Finally, return a list of the results, converting each to type numeric. 
  #
  list(x=as.numeric(x), gradient=as.numeric(eval(dh.dx)), 
       hessian=as.numeric(eval(d2h.dx2)), N.iter=as.numeric(Iter))
}