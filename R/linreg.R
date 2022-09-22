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
                          .self$x <-  model.matrix(formula, data)
                          .self$y <-  data[[formula[[2]]]]
                          .self$data <- data
                          .self$formula <-  formula
                          n <- nrow(.self$x)
                          p <- ncol(.self$x)
                          .self$df <- n-p
                          .self$hat_B <- solve((t( .self$x)%*% .self$x))%*%t( .self$x)%*% .self$y
                          .self$hat_y <- ( .self$x)%*%( .self$hat_B)
                          .self$hat_e <-  .self$y -  .self$hat_y 
                          .self$hat_var_2 <- (t(.self$hat_e)%*%.self$hat_e)/.self$df
                          
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
                          mod_object <- lm(.self$formula, data = .self$data)
                          print(mod_object)
                        }
                        
                        
                        
                      )
)


