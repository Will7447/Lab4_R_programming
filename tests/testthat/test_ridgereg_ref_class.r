data("iris")

test_that("ridgereg rejects errounous input", {
  expect_error(ridgereg_mod <- ridgereg$new(formula = Petal.Length~Sepdsal.Width+Sepal.Length, data=iris, lambda=0.1))
  expect_error(ridgereg_mod <- ridgereg$new(formula = Petal.Length~Sepal.Width+Sepal.Length, data=irfsfdis, lambda=0.1))
  expect_error(ridgereg_mod <- ridgereg$new(formula = Petal.Length~Sepal.Width+Sepal.Length, data=iris, lambda='a'))
  
})

test_that("coef() method gets similar coefficients as lm.ridge()", {
  x_lmridge <- scale(as.matrix(cbind(iris$Sepal.Length, iris$Sepal.Width)))
  y_lmridge <- iris$Petal.Length
  ridge_model <- MASS::lm.ridge(y_lmridge~x_lmridge, lambda = 0.1)
  coefficients <- coef(ridge_model) 
  
  ridgereg_mod <- ridgereg$new(Petal.Length~Sepal.Length+Sepal.Width, data=iris, lambda=0.1)
  
  error <- abs(ridgereg_mod$coef()-coefficients) 
 
  expect_true(all(error < 1e-02))
})


