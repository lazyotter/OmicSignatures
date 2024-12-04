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
