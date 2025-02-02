#' Drug Repurposing
#' Diffuse the input vector using an enhancement matrix
#' 
#' @param input.vector A numeric vector to be diffused
#' @param input.matrix The enhancement matrix
#' @param beta A parameter controlling the diffusion intensity
#' @param iter.max Maximum number of iterations (default: 10)
#' @param tol Stopping criterion for convergence (default: 10^(-4))
#' @return A numeric vector after diffusion
#' 
#' @examples 
#' input.vec <- c(1,3)
#' enhance.mat <- matrix(c(1,0,0,1), nrow = 2)
#' diffused.vec <- my_diffusion_function(input.vec, enhance.mat, beta = 0.5, iter.max = 10, tol = 1e-4)
#' 
#' @export
diffus_vec=function(input.vector,input.enhance.matrix,beta=0.75,iter.max=10,tol=10^(-4)){
  
        require(SMUT)
      
        
        vector.new=vector.old=input.vector
        
        iter=1
        
        repeat{
                matrix.old=vector.new
                vector.new <- beta*eigenMapMatMult(input.enhance.matrix,matrix.old)+(1-beta)*(input.vector)
                diff.vector=sum(abs((vector.new-matrix.old))) #diff.vector2=max(abs((vector.new-matrix.old)))
               
                if (iter>iter.max|diff.vector < tol){break}
                iter=iter+1
        }
        
        
        
        return(vector.new)
}







