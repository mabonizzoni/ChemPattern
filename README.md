# ChemPattern
Wolfram Language routines for chemical pattern recognition. `ChemPattern` is a Mathematica library designed to assist with chemometric analysis, particularly focusing on Linear Discriminant Analysis (LDA) and Principal Component Analysis (PCA). The library contains a set of functions for data visualization, cleaning, and statistical analysis, packaged in a single file (`ChemPattern.m`).

---

### Functions

The `ChemPattern.m` file provides the following functions:

* **`overview`**: Provides a quick quality examination of each measurement in a dataset through sparklines.
* **`outlierPCA`**: Automatically identifies outliers in a dataset using PCA scores.
* **`removeOutliers`**: A utility function to remove specific outlier points from a dataset.
* **`selectVarSubsets`**: Generates a new dataset by selecting measurements whose names match specified criteria.
* **`addLabels`**: Adds row and/or column labels to a matrix.
* **`filterVars`**: An interactive tool to explore and remove variables from an LDA analysis based on their contribution.
* **`groupcontribs`**: Generates a bar chart showing the contributions of variable groups to the overall discrimination in an LDA.
* **`heatmap`**: Produces a heat map of the entire dataset to quickly visualize data quality and distribution.
* **`pairwiseScatterPlot`**: Creates scatter plots for each pair of variables in a dataset.
* **`projectorLDA`**: Projects data from a new dataset onto the transformation rules obtained from a standard LDA on a different dataset.
* **`retainedInfo`**: Calculates the percentage of information retained as a function of filtering threshold.
* **`lda`**: Performs Linear Discriminant Analysis (LDA) on a dataset, with options for outputting scores, eigenvectors, or various plots.
* **`pca`**: Performs Principal Component Analysis (PCA) on a dataset.

### Usage

To use this library, you can load the `chempattern.m` file into a Mathematica session using the `Get` command:

```mathematica
Get["path/to/your/ChemPattern.m"]
