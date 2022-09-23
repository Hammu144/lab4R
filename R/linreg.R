#' Linreq
#'
#' @field x matrix.
#' @field y numeric.
#' @field hat_B matrix.
#' @field hat_y matrix.
#' @field hat_e matrix.
#' @field df integer.
#' @field hat_var_2 matrix.
#' @field hat_var_hat_B matrix.
#' @field data data.frame.
#' @field formula formula.
#' @field t_B numeric.
#'
#' @return Linear Regresion
#' @export linreg
#' @exportClass linreg
#' @importFrom methods new setRefClass
#'
#' @examples
#' data(iris)
#' mod_object <- lm(Petal.Length~Species, data = iris)
#' print(mod_object)


linreg <- setRefClass("linreg",
                      fields = list(x = "matrix", y = "numeric",
                                    hat_B = "matrix", hat_y = "matrix",
                                    hat_e = "matrix", df = "integer",
                                    hat_var_2 = "matrix", hat_var_hat_B = "matrix",
                                    data = "data.frame", formula = "formula",
                                    t_B = "numeric"
                      ),
                      methods = list(
                        initialize = function(formula, data){
                          .self$data <- data
                          .self$formula <-  formula
                          .self$x <-  model.matrix(.self$formula, .self$data)
                          .self$y <-  .self$data[[.self$formula[[2]]]]
                          n <- nrow(.self$x)
                          p <- ncol(.self$x)
                          .self$df <- n-p
                          .self$hat_B <- solve((t( .self$x)%*% .self$x))%*%t( .self$x)%*% .self$y
                          .self$hat_y <- ( .self$x) %*% ( .self$hat_B)
                          .self$hat_e <-  .self$y -  .self$hat_y
                          .self$hat_var_2 <- (t(.self$hat_e)%*%.self$hat_e) / .self$df

                          .self$hat_var_hat_B <- ((.self$hat_var_2)^(2))[1,1] * (solve(t(.self$x)%*%.self$x))

                          for (i in 1:length(.self$hat_B))
                          {
                            tem_vec <- .self$hat_B[i] /(sqrt(.self$hat_var_hat_B[i,i]))
                            .self$t_B = append(.self$t_B,tem_vec)

                          }


                        },
                        plot = function(){
                          base::print(ggplot2::ggplot(data = .self$data, mapping= ggplot2::aes(x= .self$hat_y, y=.self$hat_e))+
                                                  ggplot2::geom_point()+
                                                  ggplot2::stat_summary(fun=median, geom = "line"))
                          y_axis = sqrt(abs( (.self$hat_e  -  mean(.self$hat_e)) / sd(.self$hat_e)))

                          base::print(ggplot2::ggplot(data = .self$data, mapping= ggplot2::aes(x= .self$hat_y, y=y_axis ))+
                                        ggplot2::geom_point()+
                                        ggplot2::stat_summary(fun=median, geom = "line"))


                        },
                        resid = function(){
                          return(.self$hat_e)
                        },

                        pred = function(){
                          return(.self$hat_y)
                        },

                        coef = function(){
                          temp_vec <-    as.vector(.self$hat_B)
                          names(temp_vec)<-
                            c("Intercept", "Speciesversicolor", "Speciesvirginica")
                          return(temp_vec)
                        },

                        summary = function(){
                          v1 <- .self$hat_B
                          v2 <-sqrt(abs(diag(.self$hat_var_hat_B)))
                          v3 <- .self$t_B
                          v4 <- c()
                          for(i in 1:length(.self$t_B)){
                            p = pt(.self$t_B[i], .self$df, lower.tail = FALSE, log.p = FALSE)
                            if((p < 0.01) & (p < 0.05)){v4 <- append(v4,"***")}
                            else if((p > 0.01) & (p < 0.05)){v4 <-append(v4,"**")}
                            else if((p > 0.05) & (p < 0.05)){v4 <-append(v4,"*")}
                            else {v4 <-append(v4," ")}

                          }
                           m <- matrix(c(v1,v2,v3,v4 ),nrow = 3, ncol = 4, byrow = FALSE)
                           rownames(m) <-
                             c("Intercept", "Speciesversicolor", "Speciesvirginica")
                           colnames(m) <-
                             c("Estimate", "Std.Error", "t value", "pr(>|t|)")
                           return(m)



                          cat(paste0("Residual standard error: ", sqrt(.self$hat_var_2), " on ", .self$df , " degrees of freedom"))

                        },
                        print = function(){
                        n_formula <- as.character(.self$formula)
                        cat(paste0("lm(formula =" ,n_formula[2], " " , n_formula[1], " " ,n_formula[3], " data = iris)"))

                        cat(paste0(rownames(.self$hat_B), "/n"))
                        cat(paste0(hat_B))



                        }



                      )
)



# data(iris)
# mod_object <- lm(Petal.Length~Species, data = iris)
# print(mod_object)
# summary(mod_object)
#
#

