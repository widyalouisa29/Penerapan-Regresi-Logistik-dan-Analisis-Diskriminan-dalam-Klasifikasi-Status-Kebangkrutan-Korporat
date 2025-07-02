# 1. Instalasi (cukup sekali)
install.packages("FactoMineR")
install.packages("factoextra")

# 2. Load libraries
library(MVN)
library(biotools) 
library(candisc)
library(MASS) 
library(caret)
library(readxl)
library(dplyr)
library(FactoMineR)
library(factoextra)
library(psych)
library(candisc)

# 3. Baca data
bankruptcy_italian_companies_2023 <- read_excel("bankruptcy_italian_companies_2023.xlsx")
View(bankruptcy_italian_companies_2023)
str(bankruptcy_italian_companies_2023)
summary(bankruptcy_italian_companies_2023)

# 4. Hapus kolom non-numerik ('Company')
data_numeric <- bankruptcy_italian_companies_2023[, -1]

# 5. Ubah karakter ke numerik
data_numeric <- as.data.frame(lapply(data_numeric, function(x) as.numeric(as.character(x))))

# 6. Hapus kolom dengan NA & SD = 0
data_clean <- data_numeric[, colSums(is.na(data_numeric)) == 0]
data_clean <- data_clean[, sapply(data_clean, sd) != 0]

# 7. Simpan target 'Bankrupt'
data_target <- data_clean$Bankrupt

# 8. Buang kolom target dari PCA input
data_pca_input <- data_clean[, !colnames(data_clean) %in% "Bankrupt"]

#Standarisasi
data_scaled <- scale(data_pca_input)


#8. Check KMO
r <- cor(data_scaled)
KMO(r)

#9. Bartlett Test
bartlett.test(data_scaled)

#---------Principal Component Analysis-----------
##function in R----
#prcomp
PCA.mod <- prcomp(data_scaled)
summary(PCA.mod)       #vs t(cumvar)
load <-PCA.mod$rotation       #vs pc$vectors
head(PCA.mod$x)        #vs head(scores)


# Hitung proporsi variansi yang dijelaskan oleh setiap PC
explained_var <- PCA.mod$sdev^2
pve <- explained_var / sum(explained_var)  # proportion of variance explained

# Plot akumulasi variansi (curve style)
plot(cumsum(pve), 
     xlab = "Principal Component", 
     ylab = "Accumulative Prop. of Variation Explained",
     ylim = c(0, 1), 
     type = 'b',        # garis dan titik
     pch = 19,          # bentuk titik bulat
     col = "blue",
     main = "Cumulative Variance Explained (PCA)")
abline(h = 0.8, col = "red", lty = 2)  

#using library FactoMineR
# https://rpubs.com/cahyaalkahfi/pca-with-r
library(FactoMineR)
pca_result <- PCA(data_scaled, 
                  scale.unit = TRUE, 
                  graph = FALSE, 
                  ncp=10)     # ncp kita set agar menghasilkan output semua dimensi (default 5)

fviz_pca_biplot(pca_result,
                geom.ind = "point",              # Tampilkan titik untuk observasi
                col.ind = "steelblue",           # Warna observasi
                fill.ind = "lightblue",          # Warna isi titik observasi (jika shape = 21)
                col.var = "red",                 # Warna panah variabel
                repel = TRUE,                    # Hindari label saling tumpang tindih
                labelsize = 2,                   # Ukuran label variabel
                arrowsize = 0.8,                 # Ukuran panah
                addEllipses = FALSE,             # Jangan tampilkan ellipse (opsional)
                title = "Biplot PCA: 15 Variabel Terpenting",
                select.var = list(contrib = 15)  # Hanya tampilkan 15 variabel dengan kontribusi terbesar
)
# menampilkan ringkasan hasil pca
eig_values <- pca_result$eig          # vs print(cumvar)
pca_result$svd$V        # vs pc$vectors
pca_result$ind['coord'] # vs head(scores)


#Visualisasi
library(factoextra)
# membuat scree plot
fviz_eig(pca_result, 
         addlabels = TRUE, 
         ncp = 15, 
         barfill = "skyblue", 
         barcolor = "darkblue", 
         linecolor = "red")

# 10. Ambil 35 komponen utama
data_reduced <- as.data.frame(pca_result$ind$coord[, 1:10])
data_reduced$Bankrupt <- data_target


#=========================Uji Asumsi========================================
# Mahalanobis distance pada hasil PCA
mahal <- mahalanobis(data_reduced[, 1:10],
                     colMeans(data_reduced[, 1:10]),
                     cov(data_reduced[, 1:10]))

# Tentukan cutoff (misalnya alpha = 0.001 untuk outlier ekstrim)
cutoff <- qchisq(0.999, df = 10)

# Identifikasi indeks outlier
outlier_index <- which(mahal > cutoff)

# Buang outlier
data_reduced_no_outlier <- data_reduced[-outlier_index, ]
data_target_no_outlier <- data_reduced_no_outlier$Bankrupt
data_reduced_no_outlier$Bankrupt <- NULL

# Tambahkan kembali target ke data PCA yang sudah dibersihkan dari outlier
data_reduced_no_outlier$Bankrupt <- as.factor(data_target_no_outlier
str(data_reduced_no_outlier)

#UJI ASUMSI MULTIVARIATE NORMALITY
# Hanya 10 komponen utama tanpa kolom target
mardia_result <- mvn(data_reduced_no_outlier[, 1:10], mvnTest = "mardia", alpha = 0.05)
print(mardia_result)

#UJI ASUMSI HOMOGENITAS RAGAM PERAGAM 
uji_bart <- function(x){ 
  method <- "Bartlett's test of sphericity"
  data.name <- deparse(substitute(x))
  x <- subset(x, complete.cases(x))
  n <- nrow(x)
  p <- ncol(x)
  chisq <- (1 - n + (2 * p + 5) / 6) * log(det(cor(x)))
  df <- p * (p - 1) / 2
  p.value <- pchisq(chisq, df, lower.tail = FALSE)
  names(chisq) <- "Chi-squared"
  names(df) <- "df"
  return(structure(list(statistic = chisq, parameter = df, p.value = p.value,
                        method = method, data.name = data.name), class = "htest"))
}

uji_bart(data_reduced_no_outlier[, 1:10])

#PENGUJIAN PEREBEDAAN RATA-RATA VARIABLE DEPENDEN 
manova_model <- manova(as.matrix(data_reduced_no_outlier[, 1:10]) ~ Bankrupt, data = data_reduced_no_outlier)
summary(manova_model, test = "Wilks")

#KONTRIBUSI VARIABEL PREDIKTOR TERHADAP VARIABEL RESPON
cc <- candisc(manova_model)
cc

#=======================#Analisis Diskriminan#=============================
# 19. Split data untuk LDA 
set.seed(555)
ind <- sample(2, nrow(data_reduced_no_outlier), replace = TRUE, prob = c(0.6, 0.4))
ind
train_data <- data_reduced_no_outlier[ind == 1, ]
summary(train_data)
test_data <- data_reduced_no_outlier[ind == 2, ]
View(test_data)

# 20. Linear Discriminant Analysis (LDA)
lda_model <- lda(Bankrupt ~ ., data = train_data)
lda_model
summary(lda_model)

table(train_data$Bankrupt)
table(test_data$Bankrupt)

# 23. Partimat
library(klaR)
partimat(Bankrupt ~ ., data = train_data, method = "lda")

# 21. Prediksi dan Histogram
p <- predict(lda_model, train_data)


lda_df <- data.frame(LD1 = p$x,  # hanya LD1 karena hanya ada 1 variabel diskiriminan
                     Class = as.factor(train_data$Bankrupt))

# Biplot LD1 (karena hanya ada 1 variabel diskriminan)
ggplot(lda_df, aes(x = LD1, fill = Class)) +
  geom_density(alpha = 0.6) +
  labs(title = "LDA Biplot (LD1)",
       x = "Linear Discriminant 1 (LD1)",
       y = "Density",
       fill = "Status Bangkrut") +
  scale_fill_manual(values = c("0" = "#00AFBB", "1" = "#FC4E07"),
                    labels = c("0" = "Tidak Bangkrut", "1" = "Bangkrut")) +
  theme_minimal()

# Scatter plot hanya LD1 (gunakan jitter agar tidak tumpuk)
ggplot(lda_df, aes(x = LD1, y = 0, color = Class)) +
  geom_jitter(height = 0.05, size = 2, alpha = 0.7) +
  labs(title = "LDA Scatter Plot (LD1)",
       x = "Linear Discriminant 1 (LD1)",
       y = "",
       color = "Bankrupt") +
  theme_minimal()


# 24. Confusion Matrix & Accuracy (Train)
p1 <- predict(lda_model, train_data)$class
conf_train <- table(Predicted = p1, Actual = train_data$Bankrupt)
conf_train
cat("Accuracy (Train):", sum(diag(conf_train)) / sum(conf_train), "\n")

# 25. Confusion Matrix & Accuracy (Test)
p2 <- predict(lda_model, test_data)$class
conf_test <- table(Predicted = p2, Actual = test_data$Bankrupt)
conf_test
cat("Accuracy (Test):", sum(diag(conf_test)) / sum(conf_test), "\n")
#======================#regresi logistik#================================
library(readr)
library(generalhoslem)
library(pscl)
library(car)

# 7. Simpan target 'Bankrupt'
Y <- as.factor(data_clean$Bankrupt)

#Membentuk data frame
data_logistik<-data.frame(data_reduced,Y)
str(data_logistik)

#Memeriksa Asumsi Nonmultikolinieritas
model_vif <- lm(as.numeric(Y) ~ ., data = data_logistik)
vif(model_vif)

Y <- as.factor(data_clean$Bankrupt)
# Set seed untuk reproducibility
set.seed(123)

# Buat indeks pembagian 70% training
n <- nrow(data_logistik)
index_train <- sample(1:n, size = 0.7 * n)

# Pisahkan data training dan testing
train_data <- data_logistik[index_train, ]
test_data <- data_logistik[-index_train, ]

# 6. Fit model regresi logistik
model_train <- glm(Y ~ ., data = train_data, family = binomial())
summary(model_train)

#Uji Signifikansi Keseluruhan Model
pR2(model_train)
qchisq(0.95,9)

#R Square
Rsq<-1-(994.67/1414.57)

#Odds Ratio
beta<-(coef(model_train))
beta

OR_beta<-exp(beta)
OR_beta

cbind(beta,OR_beta)

# 5. Prediksi probabilitas di data test
yp_hat_test <- predict(model_train, newdata = test_data, type = "response")
# 6. Tambahkan prediksi ke test data
test_data$yp_hat <- yp_hat_test
# 7. Klasifikasi (cutoff 0.5)
pred_class <- ifelse(yp_hat_test > 0.5, 1, 0)
pred_class <- as.factor(pred_class)
# Pastikan label asli berupa faktor juga
test_data$Y <- as.factor(test_data$Y)
# 8. Confusion Matrix
table(Predicted = pred_class, Actual = test_data$Y)
# Akurasi
accuracy <- mean(pred_class == as.numeric(test_data$Y) - 1)  # karena Y adalah faktor
print(paste("Akurasi data test:", round(accuracy, 3)))
#  Uji goodness of fit jika perlu
library(generalhoslem)
logitgof(test_data$Y, yp_hat_test) 

# Sensitivitas dan spesifisitas
true_pos <- sum(test_data$Y == 1 & pred_class == 1)
false_neg <- sum(test_data$Y == 1 & pred_class == 0)
sensitivitas_logit <- true_pos / (true_pos + false_neg)

true_neg <- sum(test_data$Y == 0 & pred_class == 0)
false_pos <- sum(test_data$Y == 0 & pred_class == 1)
spesifisitas_logit <- true_neg / (true_neg + false_pos)

# Tampilkan hasil
cat("Akurasi Regresi Logistik:", round(accuracy, 4), "\n")
cat("Sensitivitas:", round(sensitivitas_logit, 4), "\n")
cat("Spesifisitas:", round(spesifisitas_logit, 4), "\n")

#perbandingan
# Buat data frame perbandingan metrik
perbandingan_metrik <- data.frame(
  Model = c("QDA", "Regresi Logistik"),
  Akurasi = c(round(akurasi, 4), round(akurasi_logit, 4)),
  Sensitivitas = c(round(sensitivitas, 4), round(sensitivitas_logit, 4)),
  Spesifisitas = c(round(spesifisitas, 4), round(spesifisitas_logit, 4))
)

# Pastikan package ROCR terpasang
install.packages("ROCR")
library(ROCR)

# Buat objek prediksi dan performa
pred_roc <- prediction(yp_hat_test, test_data$Y)
perf_roc <- performance(pred_roc, "tpr", "fpr")

# Plot ROC curve
plot(perf_roc, col = "blue", lwd = 2, main = "ROC Curve - Logistic Regression")
abline(a = 0, b = 1, lty = 2, col = "gray")  # Garis diagonal random classifier

# Hitung AUC
auc <- performance(pred_roc, "auc")
auc_value <- auc@y.values[[1]]
cat("AUC:", round(auc_value, 4), "\n")

