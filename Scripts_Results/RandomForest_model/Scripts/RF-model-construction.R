# Iris-RandomForest

library(randomForest)
library(missForest)
library(datasets)
library(caret)

data <- iris
str(data)  # str()显示对象的结构

data$Species <- as.factor(data$Species)  # 转换为因子类型（分类变量）
table(data$Species)

set.seed(123)

ind <- sample(2,nrow(data),replace = TRUE,prob = c(0.7,0.3))
# smple(2) 从1：2中抽样，nrow()抽样次数等于数据行数，replace = TRUE有放回抽样
#prob=c(0.7,0.3)抽样概率，1的概率70%，2的概率30%
train <- data[ind == 1,]
test <- data[ind == 2,]

rf <- randomForest(Species~., data = train, proximity = TRUE)
# Species~. ：公式，Species为因变量，. 表示所有其他变量为自变量。
# proximity = TRUE 计算样本间的邻近度矩阵。

tuned_rf_model <- tuneRF(train[,-5],train[,5],ntreeTry = 500,
                         stepFactor = 1.5, improve = 0.01,trace = TRUE,
                         plot = TRUE)

# train[,-5]训练集特征(排除第五列)，train[,5]训练集标签(第五列)，
# ntreeTry = 500 尝试的树数量，stepFactor mtry调整步长因子，improve = 0.01，最小改进阈值
# trace = TRUE 显示调参过程；plot = TRUE 绘制OOB误差图

# k折交叉验证，将训练集分为k个子集，每次使用其中k-1个子集作为训练数据，
# 剩余一个子集作为验证数据，重复进行k次并计算平均性能指标。

#R中的caret包用于调整机器学习算法。并非所有的机器学习算法参数都可以在caret中进行调整，
#只有影响较大的算法参数才可以在caret中调整。RandomForest算法只有mtry参数可调节。

rf_grid <- expand.grid(mtry = seq(2, 6, by = 1))
# expand.grid() 创建参数网络，测试mtry值2,3,4,5,6
control <- trainControl(method = "cv",number = 10)
# trainCotrol控制训练过程，10cv，10折交叉验证。
rf_model <- train(Species ~ ., data = train, method = "rf", trControl = control,
                  tuneGrid = rf_grid,importance = TRUE)
# tuneGrid调参网格；importance = TRUE 计算变量重要性。

print(rf_model$bestTune)

print(rf_model$results)

best_model <- rf_model$finalModel


#预测测试集内容
predictions <- predict(best_model,newdata = test)

#绘制混淆矩阵
result_matrix <- confusionMatrix(table(predictions,test$Species))
print(result_matrix)
result_matrix$byClass

#预测各个类别的概率
rf_pred <- predict(best_model,test,type = 'prob')
rf_pred <- data.frame(rf_pred)
colnames(rf_pred) <- paste0(colnames(rf_pred),"_pred_RF")

#
