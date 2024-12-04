# OmicSignatures: Generate omic-based signatures using LASSO regression
OmicSignatures is an R package designed to create omic-based signatures for various exposures, such as dietary and lifestyle factors, using cross-validated LASSO regression. The package includes robust methods for nested cross-validation to evaluate model performance and ensure the reliability of generated signatures. Additionally, it provides tools to correct for covariates in features or exposures either before running the LASSO or within the LASSO model itself.

## Installation
```R
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}

# Install OmicSignatures
devtools::install_github("lazyotter/OmicSignatures")
```
## Definitions
Exposure: Refers to any factor that influences an individual’s biological system, such as dietary intake, physical activity, environmental pollutants, or lifestyle habits. Exposures can be continuous variables (e.g., daily calcium intake) or categorical variables (e.g., smoking status: smoker/non-smoker).
Omics: Represents high-dimensional biological data generated through technologies such as genomics, transcriptomics, proteomics, or metabolomics. These datasets typically include measurements of thousands of biological features (e.g., genes, proteins, metabolites) across samples.

## Documentation
Detailled documentation is available throughout the package. You can access it by adding ? before the function name. For example,
```R
?create_signature
```

## Usage
```R
library(OmicSignatures)

# Load dataset (participants/samples in rows, features/covariates/exposure in columns)
all_data <- read.csv("omics_covars.csv")

# Get names of all features
features <- colnames(all_data)[3:500]

# Covariates
covars <- c("Age", "Sex", "Smoke_Stat", "BMI")

# Ensure covariates format is correct (categorical as factor, numeric as numeric)
all_data <- covar_to_factor(all_data, c("Sex", "Smoke_Stat")
all_data <- covar_to_numeric(all_data, c("Age", "BMI")


## Generate omic signature of exposure
signature <- create_signature(data = all_data, train_idx = c(1:nrow(all_data)), features = features, exposure = "Coffee", covars = covars)

# Get model performance of LASSO
model_perf <- nested_cv_signature(data = data, features = features, exposure = "Coffee", covars = covars)
```
Other methods can be used to correct for covariates. See ```?residuals``` for more information.
