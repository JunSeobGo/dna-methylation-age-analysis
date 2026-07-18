# Ridge 시각화
lambda_values <- c(0.1, 0.5, 1)


for (l in lambda_values) {
  fit <- glmnet(x_train, y_train, alpha = 0, lambda = l)
  
  pred <- predict(fit, newx = x_test, s = l)
  
  # Ridge 회귀 모델의 계수 확인
  coef(fit, s = 0.1)  # s는 lambda 값
  
  mse <- mean((y_test - pred)^2)
  print(paste("Mean Squared Error:", mse))
  
  result[[as.character(l)]] <- data.frame(pred = as.vector(pred), y_test = y_test, lambda = l)
  
  p <- ggplot(result[[as.character(l)]], aes(x = pred, y = y_test)) +
    geom_point(colour = 'red') +
    stat_smooth(method = lm) +
    labs(x = "Predicted age", y = "Chronological age", title = paste0("Ridge Regression (lambda = ", l, ")")) +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  
  print(p)
}

# Elastic net 시각화
lambda_values <- c(0.1, 0.5, 1)

result <- list()
for (l in lambda_values) {
  fit <- glmnet(x_train, y_train, alpha = 0.5, lambda = l)
  
  pred <- predict(fit, newx = x_test, s = l)
  
  # Elastic net 회귀 모델의 계수 확인
  coef(fit, s = 0.1)  # s는 lambda 값
  
  mse <- mean((y_test - pred)^2)
  print(paste("Mean Squared Error:", mse))
  
  result[[as.character(l)]] <- data.frame(pred = as.vector(pred), y_test = y_test, lambda = l)
  
  p <- ggplot(result[[as.character(l)]], aes(x = pred, y = y_test)) +
    geom_point(colour = 'red') +
    stat_smooth(method = lm) +
    labs(x = "Predicted age", y = "Chronological age", title = paste0("Elastic Net (lambda = ", l, ")")) +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  
  print(p)
}
# Lasso 시각화
lambda_values <- c(0.1, 0.5, 1)


for (l in lambda_values) {
  fit <- glmnet(x_train, y_train, alpha = 1, lambda = l)
  
  pred <- predict(fit, newx = x_test, s = l)
  
  # Lasso 회귀 모델의 계수 확인
  coef(fit, s = 0.1)  # s는 lambda 값
  
  mse <- mean((y_test - pred)^2)
  print(paste("Mean Squared Error:", mse))
  
  result[[as.character(l)]] <- data.frame(pred = as.vector(pred), y_test = y_test, lambda = l)
  
  p <- ggplot(result[[as.character(l)]], aes(x = pred, y = y_test)) +
    geom_point(colour = 'red') +
    stat_smooth(method = lm) +
    labs(x = "Predicted age", y = "Chronological age", title = paste0("Lasso (lambda = ", l, ")")) +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  
  print(p)
}
#---------------------------
p <- ggplot(result, aes(x = pred, y = y_test)) +
  geom_point() +
  stat_smooth(method = lm)+
  labs(x = "Predicted age", y = "Chronological age", main=paste("Ridge Regression,MSE=",round(mse ,2))) 
