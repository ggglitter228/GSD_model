setwd("G:/Bioinformation/Scripts_Results/RandomForest_model/Scripts/")

# 加载必要的包
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ROCR)

# 设置随机种子保证结果可重复
set.seed(123)

# 1. 数据准备 - 创建二分类问题
# 选择两个类别进行二分类（例如：setosa vs non-setosa）
data(iris)

# 方法1: 创建二分类标签（setosa vs 其他）
iris_binary <- iris
iris_binary$Binary_Class <- ifelse(iris$Species == "setosa", "setosa", "other")
iris_binary$Binary_Class <- as.factor(iris_binary$Binary_Class)
iris_binary$Species <- NULL  # 移除原始多分类标签

# 或者方法2: 只保留两个类别（setosa vs versicolor）
# iris_binary <- iris[iris$Species %in% c("setosa", "versicolor"), ]
# iris_binary$Binary_Class <- as.factor(iris_binary$Species)
# iris_binary$Species <- NULL

#数据探索
str(iris_binary)
table(iris_binary$Binary_Class)
sum(is.na(iris_binary))

# 2. 数据分割
train_index <- createDataPartition(iris_binary$Binary_Class, p = 0.7, list = FALSE)
# p = 0.7训练集中样本所占的比例；还可以设置times，想要分割的次数，默认为1
train_data <- iris_binary[train_index, ]
test_data <- iris_binary[-train_index, ]



# 检查训练集和测试集的类别分布
cat("训练集类别分布:\n")
print(table(train_data$Binary_Class))
cat("测试集类别分布:\n")
print(table(test_data$Binary_Class))

# 3. 基础随机森林模型 - 输出概率
rf_base <- randomForest(Binary_Class ~ .,
                        data = train_data,
                        importance = TRUE,
                        proximity = TRUE,
                        ntree = 500,
                        keep.forest = TRUE)
#importance 是否计算变量重要性，默认FALSE；proximity 是否计算临近矩阵，默认FALSE
#keep.forest是否在结果中保存森林，默认TRUE


print(rf_base)

# 4. 获取预测概率分数
# 训练集预测概率
train_prob <- predict(rf_base, train_data, type = "prob")
#type指定预测的类型，“response”返回类别预测，“prob”返回每个类别的概率，“vote”返回每棵树的投票情况
train_scores <- train_prob[, "setosa"]  # 正类的概率分数

# 测试集预测概率
test_prob <- predict(rf_base, test_data, type = "prob")
test_scores <- test_prob[, "setosa"]  # 正类的概率分数

# 5. 使用不同阈值进行预测
# 定义阈值函数
predict_with_threshold <- function(scores, threshold, positive_class = "setosa") {
  ifelse(scores >= threshold, positive_class, "other")
}

# 尝试不同的阈值
thresholds <- c(0.3, 0.5, 0.7)
results <- list()

for (thresh in thresholds) {
  test_pred <- predict_with_threshold(test_scores, thresh)
  cm <- confusionMatrix(as.factor(test_pred), test_data$Binary_Class, positive = "setosa")

  results[[paste0("threshold_", thresh)]] <- list(
    threshold = thresh,
    predictions = test_pred,
    confusion_matrix = cm,
    accuracy = cm$overall["Accuracy"],
    sensitivity = cm$byClass["Sensitivity"],
    specificity = cm$byClass["Specificity"]
  )
}

# 6. 寻找最佳阈值
# 使用ROC曲线找到最佳阈值
roc_obj <- roc(test_data$Binary_Class, test_scores, levels = c("other", "setosa"))
best_threshold <- coords(roc_obj, "best", ret = "threshold")$threshold

cat("基于ROC曲线的最佳阈值:", best_threshold, "\n")

# 使用最佳阈值进行预测
best_test_pred <- predict_with_threshold(test_scores, best_threshold)
best_cm <- confusionMatrix(as.factor(best_test_pred), test_data$Binary_Class, positive = "setosa")

# 7. 参数调优
# 使用tuneRF进行mtry参数调优
mtry_tune <- tuneRF(
  x = train_data[, -5],  # 排除Binary_Class列
  y = train_data$Binary_Class,
  ntreeTry = 500,
  stepFactor = 1.5,
  improve = 0.01,
  trace = TRUE,
  plot = TRUE
)

best_mtry <- mtry_tune[which.min(mtry_tune[,2]), 1]

# 8. 使用最优参数重新训练模型
rf_tuned <- randomForest(Binary_Class ~ .,
                         data = train_data,
                         mtry = best_mtry,
                         ntree = 1000,
                         importance = TRUE,
                         proximity = TRUE)

print(rf_tuned)

# 获取调优后模型的概率分数
test_prob_tuned <- predict(rf_tuned, test_data, type = "prob")
test_scores_tuned <- test_prob_tuned[, "setosa"]

# 9. 模型评估和比较
# 创建结果比较数据框
model_comparison <- data.frame(
  Model = c("Base RF", "Tuned RF"),
  Best_Threshold = c(best_threshold, coords(roc(test_data$Binary_Class, test_scores_tuned), "best", ret = "threshold")$threshold),
  Test_Accuracy = c(
    best_cm$overall["Accuracy"],
    confusionMatrix(
      as.factor(predict_with_threshold(test_scores_tuned,
                                       coords(roc(test_data$Binary_Class, test_scores_tuned), "best", ret = "threshold")$threshold)),
      test_data$Binary_Class, positive = "setosa")$overall["Accuracy"]
  ),
  Sensitivity = c(
    best_cm$byClass["Sensitivity"],
    confusionMatrix(
      as.factor(predict_with_threshold(test_scores_tuned,
                                       coords(roc(test_data$Binary_Class, test_scores_tuned), "best", ret = "threshold")$threshold)),
      test_data$Binary_Class, positive = "setosa")$byClass["Sensitivity"]
  ),
  Specificity = c(
    best_cm$byClass["Specificity"],
    confusionMatrix(
      as.factor(predict_with_threshold(test_scores_tuned,
                                       coords(roc(test_data$Binary_Class, test_scores_tuned), "best", ret = "threshold")$threshold)),
      test_data$Binary_Class, positive = "setosa")$byClass["Specificity"]
  )
)

print(model_comparison)

# 10. 变量重要性分析
importance_scores <- importance(rf_tuned)
print(importance_scores)

# 绘制变量重要性图
varImpPlot(rf_tuned, main = "变量重要性排序")

# 使用ggplot2绘制重要性图
importance_df <- as.data.frame(importance_scores)
importance_df$Variable <- rownames(importance_df)

ggplot(importance_df, aes(x = reorder(Variable, MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "随机森林变量重要性 (MeanDecreaseAccuracy)",
       x = "变量",
       y = "重要性得分") +
  theme_minimal()

# 11. 详细的ROC分析和绘图
# 基础模型的ROC
roc_base <- roc(test_data$Binary_Class, test_scores, levels = c("other", "setosa"))
roc_tuned <- roc(test_data$Binary_Class, test_scores_tuned, levels = c("other", "setosa"))

# 绘制ROC曲线
plot(roc_base, col = "blue", main = "ROC曲线比较")
plot(roc_tuned, col = "red", add = TRUE)
legend("bottomright",
       legend = c(paste0("Base RF (AUC = ", round(auc(roc_base), 3), ")"),
                  paste0("Tuned RF (AUC = ", round(auc(roc_tuned), 3), ")")),
       col = c("blue", "red"), lwd = 2)

# 12. 概率分数分布可视化
# 创建分数分布数据框
score_df <- data.frame(
  Score = c(test_scores, test_scores_tuned),
  Model = rep(c("Base RF", "Tuned RF"), each = length(test_scores)),
  True_Class = rep(test_data$Binary_Class, 2)
)

# 绘制分数分布密度图
ggplot(score_df, aes(x = Score, fill = True_Class)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Model) +
  labs(title = "模型预测分数分布",
       x = "预测分数 (setosa概率)",
       y = "密度") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1")

# 13. 阈值性能分析
# 分析不同阈值下的性能
threshold_range <- seq(0.1, 0.9, by = 0.05)
performance_df <- data.frame()

for (thresh in threshold_range) {
  pred <- predict_with_threshold(test_scores_tuned, thresh)
  cm <- confusionMatrix(as.factor(pred), test_data$Binary_Class, positive = "setosa")

  performance_df <- rbind(performance_df, data.frame(
    Threshold = thresh,
    Accuracy = cm$overall["Accuracy"],
    Sensitivity = cm$byClass["Sensitivity"],
    Specificity = cm$byClass["Specificity"],
    F1_Score = cm$byClass["F1"]
  ))
}

# 绘制阈值性能曲线
performance_long <- performance_df %>%
  pivot_longer(cols = -Threshold, names_to = "Metric", values_to = "Value")

ggplot(performance_long, aes(x = Threshold, y = Value, color = Metric)) +
  geom_line(size = 1) +
  geom_vline(xintercept = model_comparison$Best_Threshold[2], linetype = "dashed", color = "red") +
  labs(title = "不同阈值下的模型性能",
       x = "分类阈值",
       y = "性能得分") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")

# 14. 输出最终预测结果（包含分数）
final_predictions <- data.frame(
  Actual_Class = test_data$Binary_Class,
  Predicted_Score = test_scores_tuned,
  Predicted_Class = predict_with_threshold(test_scores_tuned, model_comparison$Best_Threshold[2])
)

# 添加预测置信度
final_predictions$Confidence <- ifelse(
  final_predictions$Actual_Class == final_predictions$Predicted_Class,
  "Correct",
  "Incorrect"
)

print(head(final_predictions, 10))

# 15. 保存模型和结果
results <- list(
  base_model = rf_base,
  tuned_model = rf_tuned,
  best_threshold = model_comparison$Best_Threshold[2],
  test_predictions = final_predictions,
  performance_analysis = performance_df,
  variable_importance = importance_df,
  roc_curve = roc_tuned
)

# 输出总结报告
cat("\n随机森林二分类模型总结报告:\n")
cat("================================\n")
cat("最佳mtry参数:", best_mtry, "\n")
cat("最佳分类阈值:", round(model_comparison$Best_Threshold[2], 3), "\n")
cat("测试集准确率:", round(model_comparison$Test_Accuracy[2], 4), "\n")
cat("敏感度:", round(model_comparison$Sensitivity[2], 4), "\n")
cat("特异度:", round(model_comparison$Specificity[2], 4), "\n")
cat("AUC值:", round(auc(roc_tuned), 4), "\n")
cat("\n最重要的变量:", importance_df$Variable[which.max(importance_df$MeanDecreaseAccuracy)], "\n")

# 16. 使用PR曲线评估（可选）
library(PRROC)
pr_curve <- pr.curve(scores.class0 = test_scores_tuned[test_data$Binary_Class == "setosa"],
                     scores.class1 = test_scores_tuned[test_data$Binary_Class == "other"],
                     curve = TRUE)
plot(pr_curve, main = "Precision-Recall曲线")
