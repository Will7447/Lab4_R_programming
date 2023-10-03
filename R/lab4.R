#' A linear regression function
#' @description
#' A RC class function that can handle linear regression problems like function lm(), it also handles
#' the special functions print(), plot(), pred(), resid(), coef() and summary(). 
#' @field beta_hat Regressions coefficients 
#' @field X independent variables
#' @field y dependent variable
#' @field y_hat itted values
#' @field e_hat residuals
#' @field d_f degrees of freedom
#' @field sigma_hat_2 residual variance
#' @field var_beta_hat variance of the regression coefficients
#' @field t_beta t-values for each coefficient
#' @field p_value p_values for each coefficient
#' @field for_mula formula to linear regression
#' @field data data frame to analyse
#' @field df_name name of the analysed data frame
#' 
#' @importFrom methods new
#' @importFrom plyr is.formula
#' @import ggplot2
#' 
#' @export linreg
#' @exportClass linreg

linreg <- setRefClass('linreg',
                      fields = list(beta_hat='matrix',
                                    X='matrix',
                                    y='numeric',
                                    y_hat='matrix',
                                    e_hat='matrix',
                                    d_f='numeric',
                                    sigma_hat_2='numeric',
                                    var_beta_hat='matrix',
                                    t_beta='matrix',
                                    p_value='matrix',
                                    for_mula='formula',
                                    data='data.frame',
                                    df_name='character'),
                      methods = list(
                        initialize = function(for_mula,data){
                          stopifnot(plyr::is.formula(for_mula),is.data.frame(data))
                          for_mula <<- for_mula
                          data <<- data
                          df_name <<- deparse(substitute(data))
                          X <<- model.matrix(for_mula,data)
                          y <<- data[,all.vars(for_mula)[1]]
                          beta_hat <<- solve(t(X)%*%X)%*%t(X)%*%y
                          y_hat <<- X%*%beta_hat
                          e_hat <<- y-y_hat
                          d_f <<- length(y)-length(all.vars(for_mula))
                          sigma_hat_2 <<- (t(e_hat)%*%e_hat)[1,1]/d_f
                          var_beta_hat <<- sigma_hat_2*solve(t(X)%*%X)
                          t_beta <<- beta_hat/sqrt(diag(var_beta_hat))
                          p_value <<- 2*pt(abs(t_beta),df = d_f,lower.tail = FALSE)
                          },
                        print = function(){
                          cat('Call:\n')
                          cat(paste0('linreg(formula = ',deparse(for_mula),', data = ',df_name,')\n'))
                          cat('\nCoefficients:\n')
                          base::print(beta_hat[,1])
                        },
                        plot = function(){
                          library(ggplot2)
                          index <- 1:length(y)
                          plot_data <- data.frame(index, y_hat, e_hat, sqrt(abs(e_hat/sd(e_hat))))
                          plot_data$label <- ifelse(plot_data$index %in% c(99, 118, 119), as.character(plot_data$index), "")
                          p1 <- ggplot(plot_data, aes(x = y_hat, y = e_hat)) +
                            geom_point(size = 3.5, shape = 1, stroke = 1) +
                            geom_text(aes(label = label), hjust = 1.5, vjust = 0) +
                            theme_test()+
                            labs(x = paste0("Fitted values\nlinreg(",deparse(for_mula),')'), y = "Residuals") +
                            ggtitle("Residuals vs Fitted") +
                            theme(plot.title = element_text(hjust = 0.5))+
                            geom_hline(yintercept = 0, color = "grey", linetype = "dotted")+
                            stat_summary(geom = "line", color = "red", fun = median)
                          
                          p2 <- ggplot(plot_data, aes(x = y_hat, y = sqrt(abs(e_hat/sd(e_hat))))) +
                            geom_point(size = 3.5, shape = 1, stroke = 1) +
                            geom_text(aes(label = label), hjust = 1.5, vjust = 0) +
                            theme_test()+
                            labs(x = paste0("Fitted values\nlinreg(",deparse(for_mula),')'), y = expression(sqrt(abs("Standardized residuals")))) +
                            ggtitle("Scale-Location") +
                            theme(plot.title = element_text(hjust = 0.5))+
                            stat_summary(geom = "line", color = "red", fun = median)
                          base::print(p1)
                          base::print(p2)
                        },
                        resid = function(){
                          return(e_hat)
                        },
                        pred = function(){
                          return(y_hat)
                        },
                        coef = function(){
                          vec_coef = as.vector(beta_hat)
                          names(vec_coef) = row.names(beta_hat)
                          return(vec_coef)
                        },
                        summary = function(){
                          last_col <- matrix(rep("***", ncol(X)), ncol = 1)
                          coe_mat = cbind(beta_hat,sqrt(diag(var_beta_hat)),t_beta,p_value,last_col)
                          coe_mat <- as.data.frame(coe_mat)
                          colnames(coe_mat) <- c('Estimate','Std. Error','t value','Pr(>|t|)','')
                          cat('Coefficients:\n')
                          base::print(coe_mat)
                          cat(paste0('\nResidual standard error: ',round(sqrt(sigma_hat_2),3),' on ',d_f,' degrees of freedom'))
                        }
                                     )
                      )

